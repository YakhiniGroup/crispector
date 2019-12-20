from constants_and_types import IsEdit, IndelType, AlgResult, Pr, Path, CIGAR, FREQ, IS_EDIT, C_TX, C_MOCK, TX_READ_NUM, \
    MOCK_READ_NUM, TX_EDIT, EDIT_PERCENT, CI_LOW, CI_HIGH, READ_LEN_SIDE, ALIGN_CUT_SITE, ALIGNMENT_W_INS, \
    ALIGNMENT_W_DEL, POS_IDX_E, POS_IDX_S, INDEL_TYPE, ExpType
from exceptions import ClassificationFailed
from utils import Configurator, Logger, color_edit_background, get_read_around_cut_site
from modification_tables import ModificationTables
from modification_types import ModificationTypes
from typing import List, Tuple
from scipy.stats import norm, binom, hypergeom
import numpy as np
import os
import matplotlib as mpl
from matplotlib import pyplot as plt
import warnings
import seaborn as sns #
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import pandas as pd


class CrispectorAlgorithm:
    """
    Core algorithm for CRISPECTOR
    """

    def __init__(self, site_name: str, experiment_name: str, cut_site: int, modification: ModificationTypes,
                 binom_p_l: List[Pr], confidence: Pr, on_target: bool, output: Path):
        """
        :param site_name:
        :param experiment_name:
        :param cut_site:
        :param modification:
        :param binom_p_l: binom_p list
        :param confidence: confidence interval
        :param on_target: is on target flag
        :param output:  output directory path
        """
        self._name = site_name
        self._experiment_name = "{} - {}".format(experiment_name, site_name)
        self._cut_site = cut_site
        self._modifications = modification
        self._binom_p_l = binom_p_l
        self._confidence = confidence
        self._cfg = Configurator.get_cfg()
        self._win_size = self._cfg["window_size"]
        self._tables_offset = self._cut_site - self._win_size  # position offset
        self._n_reads_tx = 0
        self._n_reads_mock = 0
        self._edit: IsEdit = dict()
        self._output = output
        self._tx_df = None
        self._mock_df = None
        self._logger = Logger.get_logger()
        self._on_target = on_target

    def evaluate(self, tables: ModificationTables) -> AlgResult:
        """
        Evaluate editing activity from ModificationTables and create a folder with results graphs and tables
        :param tables: ModificationTables
        :return: Dict with results
        """
        # set algorithm properties
        self._tx_df = tables.tx_reads
        self._mock_df = tables.mock_reads
        self._n_reads_tx = tables.n_reads_tx
        self._n_reads_mock = tables.n_reads_mock
        edited_indexes = []

        # Run evaluation on each modification type
        for table_idx, (indel_type, prior) in enumerate(zip(self._modifications.types, tables.priors)):
            table = tables.tables[table_idx]
            pointers = tables.pointers[table_idx]
            binom_p = self._binom_p_l[table_idx]
            eval_size = (2*self._win_size+1) if indel_type == IndelType.INS else 2*self._win_size
            self._edit[table_idx] = np.array(eval_size*[False])

            # Run evaluation on every position
            for pos_idx in range(eval_size):
                tx_indels, mock_indels = table[[C_TX, C_MOCK], self._tables_offset + pos_idx]
                is_edit = self._classify_position(tx_indels, mock_indels, prior[pos_idx], binom_p)
                self._edit[table_idx][pos_idx] = is_edit
                if is_edit:
                    edited_indexes += pointers[self._tables_offset + pos_idx]

        # Compute editing activity.
        result_dict = self._compute_editing_activity(edited_indexes)

        self._logger.debug("Site {} - Start creating plots and tables".format(self._name))

        # Site output - Only plot for a valid output
        if self._output is not None:
            # Mark all edited modifications as an edit in modification distribution
            tables.set_tx_dist_is_edit_col(self._edit, self._tables_offset, self._win_size)

            # plot modification table
            # TODO - decide which plots are needed
            all_table_idx = list(range(self._modifications.size))
            del_table_idx = [i for i, x in enumerate(self._modifications.types) if x == IndelType.DEL]
            ins_table_idx = [i for i, x in enumerate(self._modifications.types) if x == IndelType.INS]
            mix_table_idx = [i for i, x in enumerate(self._modifications.types) if x in [IndelType.MIXED, IndelType.SUB]]
            tables.plot_tables(self._edit, self._tables_offset, self._output, self._experiment_name, all_table_idx, "all", figsize=(14,20))
            tables.plot_tables(self._edit, self._tables_offset, self._output, self._experiment_name, del_table_idx, "del")
            tables.plot_tables(self._edit, self._tables_offset, self._output, self._experiment_name, ins_table_idx, "ins")
            tables.plot_tables(self._edit, self._tables_offset, self._output, self._experiment_name, mix_table_idx, "mix_and_sub")
            # TODO - needed? tables.dump_tables(self._edit, self._tables_offset, self._output)

            # Dump .csv file with all reads
            self._tx_df.to_csv(os.path.join(self._output, "treatment_aligned_reads.csv.gz"), index=False,
                               compression='gzip')
            self._mock_df.to_csv(os.path.join(self._output, "mock_aligned_reads.csv.gz"), index=False,
                                 compression='gzip')

            # Create edited read table
            self._plot_edited_reads_to_table(tables)

            # Plot mutation distribution
            self._plot_distribution_of_edit_events(tables)
            self._plot_distribution_of_all_modifications(tables)
            self._distribution_of_edit_event_sizes(tables)
            self._distribution_of_edit_event_sizes_3_plots(tables) # TODO - which one is better

            # Plot editing activity
            self._plot_editing_activity(result_dict)


        return result_dict

    def _classify_position(self, tx_indels: int, mock_indels: int, edit_prior: Pr, binom_p: Pr) -> bool:
        """
        Classify for each position if it Edit or not.
        :param tx_indels:
        :param mock_indels:
        :param edit_prior:
        :param binom_p:
        :return: Edit or not bool
        """

        total_reads = self._n_reads_tx + self._n_reads_mock
        no_edit_prior = 1 - edit_prior

        # if no indels in treatment, return True (same results as if run the classifier)
        if tx_indels == 0:
            return False

        # In extreme cases both posteriors are smaller than python minimum number (1e-325).
        # The heuristic to solve this rare event, is to divide by 10 both the treatment and the mock indels.
        for idx in range(5):
            total_edit = tx_indels + mock_indels

            # Likelihood function computation
            no_edit_likelihood = hypergeom.pmf(tx_indels, total_reads, total_edit, self._n_reads_tx)
            edit_likelihood = binom.pmf(k=tx_indels, n=total_edit, p=binom_p)

            # Posterior probability computation
            no_edit_post = no_edit_prior * no_edit_likelihood
            edit_post = edit_prior * edit_likelihood

            # Regular case - Both different from zero
            if (edit_post != 0) or (no_edit_post != 0):
                edit = edit_post > no_edit_post
                return edit
            else:
                tx_indels = tx_indels // 10
                mock_indels = mock_indels // 10

        raise ClassificationFailed()

    def _compute_editing_activity(self, edited_indexes: List) -> AlgResult:
        """
        - Compute editing activity.
        - Compute confidence interval
        - Fill inplace is_edit column in tx_read
        :param edited_indexes: List
        :return: AlgResult
        """
        # Compute editing activity & Fill inplace is_edit column in tx_read
        edited_indexes = list(set(edited_indexes))  # remove duplicates
        self._tx_df[IS_EDIT] = False
        self._tx_df.loc[self._tx_df.index.isin(edited_indexes), IS_EDIT] = True
        edited_reads = self._tx_df.loc[self._tx_df[IS_EDIT], FREQ].sum()
        editing_activity = edited_reads / self._n_reads_tx

        # Compute CI
        CI_low, CI_high = self._compute_confidence_interval(editing_activity)

        # Prepare return dict
        result_d: AlgResult = dict()
        result_d[TX_READ_NUM] = self._n_reads_tx
        result_d[MOCK_READ_NUM] = self._n_reads_mock
        result_d[TX_EDIT] = edited_reads
        result_d[EDIT_PERCENT] = 100*editing_activity
        result_d[CI_LOW] = 100*CI_low
        result_d[CI_HIGH] = 100*CI_high
        return result_d

    def _compute_confidence_interval(self, editing_activity: Pr) -> Tuple[Pr, Pr]:
        """
        Compute confidence interval and returns low & high CI boundary
        :param editing_activity:
        :return: Tuple of low & high CI boundary
        """
        confidence_inv = norm.ppf(self._confidence + (1-self._confidence)/2)
        half_len_CI = confidence_inv * np.sqrt((editing_activity*(1-editing_activity))/self._n_reads_tx)
        return max(0, editing_activity - half_len_CI), editing_activity + half_len_CI


    def _plot_distribution_of_all_modifications(self, tables: ModificationTables):
        tx_dist_d = dict()
        mock_dist_d = dict()
        amplicon_length = len(tables.amplicon)
        for indel_type in [IndelType.DEL, IndelType.SUB, IndelType.INS]:
            tx_dist_d[indel_type] = np.zeros(amplicon_length + 1, dtype=np.int)
            mock_dist_d[indel_type] = np.zeros(amplicon_length + 1, dtype=np.int)

        # Tx - Aggregate modifications
        for row_idx, row in tables.tx_dist.iterrows():
            mod_range = range(row[POS_IDX_S], row[POS_IDX_E])
            tx_dist_d[row[INDEL_TYPE]][mod_range] += row[FREQ]

        # Mock - Aggregate modifications
        for row_idx, row in tables.mock_dist.iterrows():
            mod_range = range(row[POS_IDX_S], row[POS_IDX_E])
            mock_dist_d[row[INDEL_TYPE]][mod_range] += row[FREQ]

        # Find plots heights
        tx_max_indels = np.array([tx_dist_d[IndelType.INS], tx_dist_d[IndelType.DEL], tx_dist_d[IndelType.SUB]]).max()
        mock_max_indels = np.array([mock_dist_d[IndelType.INS], mock_dist_d[IndelType.DEL], mock_dist_d[IndelType.SUB]]).max()
        max_indels = max(tx_max_indels, mock_max_indels)

        # Set font
        mpl.rcParams.update(mpl.rcParamsDefault)
        mpl.rcParams['font.size'] = 18
        mpl.rcParams['axes.labelsize'] = 22
        mpl.rcParams['axes.titlesize'] = 24

        # Create axes
        fig, axes = plt.subplots(nrows=2, ncols=1, figsize=(16, 9), sharex=True)
        plt.subplots_adjust(wspace=0, hspace=0)

        # Plot distributions
        for exp_type in ExpType:
            dist_d = tx_dist_d if exp_type == ExpType.TX else mock_dist_d
            ax = axes[0] if exp_type == ExpType.TX else axes[1]
            positions = list(range(amplicon_length + 1))
            ax.axvline(x=self._cut_site, linestyle='-', color="red", label='Expected cut-site', linewidth=2)
            ax.axvline(x=self._cut_site - self._win_size, linestyle='--', color='k',
                       label='Quantification window', linewidth=1)
            ax.axvline(x=self._cut_site + self._win_size, linestyle='--', color='k', linewidth=1)
            ax.plot(positions, dist_d[IndelType.INS], color='#32b165', label="Insertions", linewidth=3, alpha=0.9)
            ax.plot(positions, dist_d[IndelType.SUB], color='b', label="Substitutions", linewidth=3, alpha=0.9)
            ax.plot(positions, dist_d[IndelType.DEL], color='darkviolet', label="Deletions", linewidth=3, alpha=0.9)
            ax.set_ylim(bottom=0, top=max(int(1.1 * max_indels), 10))
            ax.set_xlim(left=0, right=amplicon_length)
            ax.set_ylabel("{}\nIndel count".format(exp_type.name()))

        # Create legend, title and x-label
        axes[0].legend()
        axes[0].set_title("Distribution Of All Modifications\n{}".format(self._experiment_name),
                     weight='bold', family='serif')

        # Remove Tx first y axe label (overlap with mock
        plt.setp(axes[0].get_yticklabels()[0], visible=False)

        with warnings.catch_warnings():
            warnings.simplefilter("ignore", UserWarning)
            fig.savefig(os.path.join(self._output, 'distribution_of_all_modifications.png'),
                        bbox_inches='tight', dpi=100)
            plt.close(fig)

    def _plot_distribution_of_edit_events(self, tables: ModificationTables):
        dist_d = dict()
        amplicon_length = len(tables.amplicon)
        for indel_type in [IndelType.DEL, IndelType.SUB, IndelType.INS]:
            dist_d[indel_type] = np.zeros(amplicon_length + 1, dtype=np.int)

        # Aggregate modifications
        edit_mod = tables.tx_dist.loc[tables.tx_dist[IS_EDIT]]
        for row_idx, row in edit_mod.iterrows():
            mod_range = range(row[POS_IDX_S], row[POS_IDX_E])
            dist_d[row[INDEL_TYPE]][mod_range] += row[FREQ]

        # Set font
        mpl.rcParams.update(mpl.rcParamsDefault)
        mpl.rcParams['font.size'] = 22
        mpl.rcParams['axes.labelsize'] = 30
        mpl.rcParams['axes.titlesize'] = 26

        # Create axes
        fig = plt.figure(figsize=(16, 9))
        ax = fig.add_axes([0, 0, 1, 1])

        positions = list(range(amplicon_length + 1))
        ax.axvline(x=self._cut_site, linestyle='-', color="red", label='Expected cut-site', linewidth=2)
        ax.axvline(x=self._cut_site - self._win_size, linestyle='--', color='k',
                   label='Quantification window', linewidth=1)
        ax.axvline(x=self._cut_site + self._win_size, linestyle='--', color='k', linewidth=1)

        for indel_type in [IndelType.INS, IndelType.SUB, IndelType.DEL]:
            ax.plot(positions, dist_d[indel_type], color=indel_type.color, label=indel_type.name, linewidth=3, alpha=0.9)

        # Create legend, axes limit and labels
        ax.legend()
        ax.set_title("Distribution Of Edited Events\n{}".format(self._experiment_name),
                     weight='bold', family='serif')
        max_indels = np.array([dist_d[IndelType.INS], dist_d[IndelType.DEL], dist_d[IndelType.SUB]]).max()
        ax.set_ylim(bottom=0, top=max(int(1.1 * max_indels), 10))
        ax.set_xlim(left=0, right=amplicon_length)
        ax.set_xlabel("Position")
        ax.set_ylabel("Indel count")

        with warnings.catch_warnings():
            warnings.simplefilter("ignore", UserWarning)
            fig.savefig(os.path.join(self._output, 'distribution_of_edit_events.png'),
                        bbox_inches='tight', dpi=100)
            plt.close(fig)

    def _plot_edited_reads_to_table(self, tables: ModificationTables):
        """
        Create an HTML page with all the edited reads (aligned to reference)
        :param tables:
        :return:
        """
        length = READ_LEN_SIDE

        # Columns variables
        column_order = list(range(1,length+1)) + list(range(length+2,2*length+4)) + [length+1]
        column_width = (2 * length + 3) * [40]
        column_width[-1] = 5  # cut-site marking
        column_width[-2] = 160  # Frequency
        column_width[-3] = 160  # number of reads

        # Header
        reference_l = [base for base in tables.amplicon[self._cut_site-length:self._cut_site+length]]
        header = []
        header += ["<b>{}<b>".format(bp) for bp in reference_l]
        header += ["<b>Number Of Reads<b>", "<b>Frequency (%)<b>", ""]
        header_colors = color_edit_background(reference_l, reference_l) + ['#c8d4e3', '#c8d4e3', 'black']

        # Prepare edit reads - Cut the relevant part and aggregate
        edit_reads = self._tx_df.loc[self._tx_df[IS_EDIT]]

        partial_reads = pd.DataFrame(columns=[FREQ, ALIGNMENT_W_INS, ALIGNMENT_W_DEL])
        if edit_reads.shape[0] > 0:
            partial_reads[FREQ] = edit_reads[FREQ].copy()
            partial_reads[ALIGNMENT_W_INS] = edit_reads.apply(lambda row : get_read_around_cut_site(row[ALIGNMENT_W_INS],
                                                              row[ALIGN_CUT_SITE], length), axis=1)
            partial_reads[ALIGNMENT_W_DEL] = edit_reads.apply(lambda row: get_read_around_cut_site(row[ALIGNMENT_W_DEL],
                                                              row[ALIGN_CUT_SITE], length), axis=1)
            partial_reads = partial_reads.groupby([ALIGNMENT_W_INS, ALIGNMENT_W_DEL])[FREQ].sum().to_frame(FREQ).reset_index()
            partial_reads = partial_reads.sort_values(by=[FREQ], ascending=False).reset_index(drop=True)

        # Cells values & colors
        values = np.zeros(shape=(partial_reads.shape[0], 2*length+3), dtype='<U20')
        colors = np.zeros(shape=(partial_reads.shape[0], 2*length+3), dtype='<U20')
        for row_idx, (_, row) in enumerate(partial_reads.iterrows()):
            read_w_ins = list(row[ALIGNMENT_W_INS])
            read_w_del = list(row[ALIGNMENT_W_DEL])
            # cell fill color
            colors[row_idx] = np.array(color_edit_background(read_w_ins, read_w_del) + ['lightgrey', 'lightgrey', 'black'])
            values[row_idx] = read_w_del + ["<b>{:,}<b>".format(row[FREQ]),
                                            "<b>{:.4f}<b>".format(100*row[FREQ]/self._n_reads_tx),""]

        # convert to list
        values = [list(values[:, idx]) for idx in range(values.shape[1])]
        colors = [list(colors[:, idx]) for idx in range(colors.shape[1])]

        # Make the plot
        fig = make_subplots(
            rows=2, cols=1,
            row_heights=[0.99, 0.01],
            vertical_spacing=0.03,
            specs=[[{"type": "table"}],
                   [{"type": "scatter"}]]
        )
        # Add the table
        fig.add_trace(
            go.Table(
                columnorder=column_order,
                columnwidth=column_width,
                header=dict(
                    values=header,
                    font=dict(size=14),
                    fill=dict(color=header_colors),
                    align="center"
                ),
                cells=dict(
                    values=values,
                    align="center",
                    fill=dict(color=colors),
                    font=dict(size=14, color='black'))
            ),
            row=1, col=1
        )

        # Add Scatter plots for legend
        fig.add_trace(
            go.Scatter(
                x=[0],
                y=[0],
                mode="lines",
                name="Insertions",
                line=dict(color='indianred')
            ),
            row=2, col=1
        )
        fig.add_trace(
            go.Scatter(
                x=[0],
                y=[0],
                mode="lines",
                name="Deletions",
                line=dict(color='lightgrey')
            ),
            row=2, col=1
        )
        fig.add_trace(
            go.Scatter(
                x=[0],
                y=[0],
                mode="lines",
                name="Substitutions",
                line=dict(color='royalblue')
            ),
            row=2, col=1
        )
        fig.add_trace(
            go.Scatter(
                x=[0],
                y=[0],
                mode="lines",
                name="Cut Site",
                line=dict(color='Black')
            ),
            row=2, col=1
        )

        fig.update_layout(
            title_text="Edited reads around the cut-site",
            title_x=0.5,
            titlefont=dict(size=30, color='black', family='Arial, sans-serif'),
            showlegend = True,
            paper_bgcolor = 'rgba(0,0,0,0)',
            plot_bgcolor = 'rgba(0,0,0,0)',
            legend = go.layout.Legend(
                traceorder="normal",
                font=dict(
                    family="sans-serif",
                    size=20,
                    color="black"
                ),
                bordercolor="Black",
                borderwidth=2
            )
        )
        fig.update_xaxes(color='rgba(0,0,0,0)')
        fig.update_yaxes(color='rgba(0,0,0,0)')

        # Dump to static HTML
        fig.write_html(os.path.join(self._output, 'edited_reads_table.html'))

    def _plot_editing_activity(self, result_d: AlgResult):
        # Set font
        mpl.rcParams.update(mpl.rcParamsDefault)
        sns.set(style="whitegrid")
        mpl.rcParams['font.size'] = 14
        mpl.rcParams['xtick.labelsize'] = 14
        mpl.rcParams['ytick.labelsize'] = 14
        mpl.rcParams['axes.labelsize'] = 16
        mpl.rcParams['axes.titlesize'] = 16

        bar_width = 0.8
        if self._on_target:
            bar_color = "#39ad48"
        else:
            bar_color = "#db5856"

        # Define fix and axes
        fig = plt.figure(figsize=(4, 4))
        ax = fig.add_axes([0, 0, 1, 1])

        # Get bar data
        editing = result_d[EDIT_PERCENT]
        CI_high = result_d[CI_HIGH] - editing
        CI_low = editing - result_d[CI_LOW]

        # plot bar
        ax.bar([0], [result_d[EDIT_PERCENT]], color=bar_color, width=bar_width,
               yerr=[[CI_low], [CI_high]], align='center', ecolor='black', capsize=15)

        # Set labels
        ax.set_ylabel("Editing Activity (%)")

        # Set scale and lim
        y_lim = max(min(1.2 * (editing + CI_high), 100), 0.1)
        ax.set_ylim(0, y_lim)
        ax.set_xlim(-0.8, 0.8)

        # Text below each bar plot + y ticks
        ax.set_xticks([0])
        ax.set_xticklabels([self._name])

        ax.set_title("Editing Activity with {} % CI\n{}".format(100 * self._confidence, self._experiment_name),
                     weight='bold', family='serif')

        ax.text(x=-0.1, y=-0.3, s="Number of edited reads\nEditing activity",
                ha='left', va='bottom', transform=fig.transFigure, family='serif')
        ax.text(x=0.55, y=-0.3, s="- {:,} (out of {:,} reads).\n"
                                  "- {:.2f}%, CI=({:.2f}%$-${:.2f}%).".format(result_d[TX_EDIT], result_d[TX_READ_NUM],
                                                                              editing, result_d[CI_LOW],
                                                                              result_d[CI_HIGH]),
                ha='left', va='bottom', transform=fig.transFigure, family='serif')

        with warnings.catch_warnings():
            warnings.simplefilter("ignore", UserWarning)
            fig.savefig(os.path.join(self._output, 'site_editing_activity.png'),
                        bbox_inches='tight', dpi=200)
            plt.close(fig)

    def _distribution_of_edit_event_sizes(self, tables: ModificationTables):
        # Set font
        mpl.rcParams.update(mpl.rcParamsDefault)
        mpl.rcParams['font.size'] = 20
        mpl.rcParams['xtick.labelsize'] = 20
        mpl.rcParams['ytick.labelsize'] = 20
        mpl.rcParams['axes.labelsize'] = 24
        mpl.rcParams['axes.titlesize'] = 26
        mpl.rcParams['legend.fontsize'] = 24
        bar_width = 1

        amplicon_length = len(tables.amplicon)

        # Create bin borders & names:
        bins = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 20, 30, 40, 50, 100, amplicon_length]
        bins_names = []
        for l, r in zip(bins[:-1], bins[1:]):
            if r - l == 1:
                bins_names.append(str(r))
            else:
                bins_names.append("{}-{}".format(l + 1, r))
        bin_num = len(bins) - 1

        # Get indel size values
        values = {}
        edited_mod = tables.tx_dist.loc[tables.tx_dist[IS_EDIT]]
        for indel in [IndelType.DEL, IndelType.INS, IndelType.SUB]:
            indel_dist = edited_mod.loc[edited_mod['IndelType'] == indel]
            df_bins = pd.cut(indel_dist['indel_length'], bins)
            values[indel] = indel_dist.groupby(df_bins)[FREQ].agg(['sum'])['sum'].values

        # Create the X position of bars
        r = {}
        for idx, indel in enumerate([IndelType.DEL, IndelType.INS, IndelType.SUB]):
            r[indel] = list(range(0 + idx, 4 * bin_num + idx, 4))

        # Create figure
        fig = plt.figure(figsize=(16, 9))
        ax = fig.add_axes([0, 0, 1, 1])

        # Create bar plot
        for indel in [IndelType.DEL, IndelType.INS, IndelType.SUB]:
            ax.bar(r[indel], values[indel], label=indel.name, color=indel.color, edgecolor='black', alpha=0.9, width=bar_width)

        # Labels, ticks and title
        plt.xlabel('Edit event length')
        plt.ylabel('Number of modifications')
        plt.legend(loc='upper right')
        ax.set_xticks(r[IndelType.INS])
        ax.set_xticklabels(bins_names, rotation='vertical')
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.set_title("Edit Event Size Distribution\n{}".format(self._experiment_name), weight='bold', family='serif')

        with warnings.catch_warnings():
            warnings.simplefilter("ignore", UserWarning)
            fig.savefig(os.path.join(self._output, 'distribution_of_edit_events_size.png'),
                    bbox_inches='tight', dpi=100)
            plt.close(fig)

    def _distribution_of_edit_event_sizes_3_plots(self, tables: ModificationTables):
        # Set font
        mpl.rcParams.update(mpl.rcParamsDefault)
        mpl.rcParams['font.size'] = 16
        mpl.rcParams['xtick.labelsize'] = 16
        mpl.rcParams['ytick.labelsize'] = 16
        mpl.rcParams['axes.labelsize'] = 18
        mpl.rcParams['axes.titlesize'] = 26
        mpl.rcParams['legend.fontsize'] = 18
        bar_width = 1

        amplicon_length = len(tables.amplicon)

        # Create bin borders & names:
        bins = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 20, 30, 40, 50, 100, amplicon_length]
        bins_names = []
        for l, r in zip(bins[:-1], bins[1:]):
            if r - l == 1:
                bins_names.append(str(r))
            else:
                bins_names.append("{}-{}".format(l + 1, r))
        bin_num = len(bins) - 1

        # Get indel size values
        values = {}
        edited_mod = tables.tx_dist.loc[tables.tx_dist[IS_EDIT]]
        for indel in [IndelType.DEL, IndelType.INS, IndelType.SUB]:
            indel_dist = edited_mod.loc[edited_mod['IndelType'] == indel]
            df_bins = pd.cut(indel_dist['indel_length'], bins)
            values[indel] = indel_dist.groupby(df_bins)[FREQ].agg(['sum'])['sum'].values

        # Create the X position of bars
        r = list(range(0 , bin_num))

        # Create figure
        fig, axes = plt.subplots(nrows=1, ncols=3, figsize=(16, 6), constrained_layout=True)

        # Create bar plot
        for idx, indel in enumerate([IndelType.DEL, IndelType.INS, IndelType.SUB]):
            axes[idx].bar(r, values[indel], label=indel.name, color=indel.color, edgecolor='black', alpha=0.9, width=bar_width)
            axes[idx].set_xticks(r)
            axes[idx].set_xticklabels(bins_names, rotation='vertical')
            axes[idx].spines['top'].set_visible(False)
            axes[idx].spines['right'].set_visible(False)

            # Labels, ticks and title
            axes[idx].set_xlabel('Edit event length')
            axes[idx].set_ylabel('Number of modifications')
            axes[idx].legend(loc='upper right')

        fig.suptitle("Edit Event Size Distribution\n{}".format(self._experiment_name), weight='bold', family='serif')

        with warnings.catch_warnings():
            warnings.simplefilter("ignore", UserWarning)
            fig.savefig(os.path.join(self._output, 'distribution_of_edit_events_size_3_plots.png'),
                    bbox_inches='tight', dpi=100)
            plt.close(fig)
