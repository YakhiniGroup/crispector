from constants_and_types import IsEdit, IndelType, AlgResult, Pr, Path, CIGAR, FREQ, IS_EDIT, C_TX, C_MOCK, TX_READ_NUM, \
    MOCK_READ_NUM, TX_EDIT, EDIT_PERCENT, CI_LOW, CI_HIGH
from input_processing import InputProcessing
from utils import Configurator, Logger
from modification_tables import ModificationTables
from modification_types import ModificationTypes
from typing import List, Tuple
from scipy.stats import norm, binom, hypergeom  # TODO - add to package requirements
import numpy as np
import os
import matplotlib as mpl
from matplotlib import pyplot as plt
import warnings
import seaborn as sns # TODO - add to package requirements

class CrispectorAlgorithm:
    """
    Core algorithm for CRISPECTOR
    """

    def __init__(self, site_name: str, cut_site: int, modification: ModificationTypes, binom_p_l: List[Pr],
                 confidence: Pr, on_target: bool, output: Path):
        """
        :param site_name:
        :param cut_site:
        :param modification:
        :param binom_p_l: binom_p list
        :param confidence: confidence interval
        :param on_target: is on target flag
        :param output:  output directory path
        """
        self._name = site_name
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
        for table_idx, (indel_type, prior) in enumerate(zip(self._modifications.types, self._modifications.priors)):
            table = tables.tables[table_idx]
            pointers = tables.pointers[table_idx]
            binom_p = self._binom_p_l[table_idx]
            eval_size = (2*self._win_size+1) if indel_type == IndelType.INS else 2*self._win_size
            self._edit[table_idx] = eval_size*[False]

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

        # Site output
        # plot modification table
        tables.plot_tables(self._edit, self._tables_offset, self._output)
        tables.dump_tables(self._edit, self._tables_offset, self._output)
        # Dump .csv file with all reads
        # TODO - add gzip, remove site_name & cigar_path?
        self._tx_df.to_csv(os.path.join(self._output, "treatment_aligned_reads.csv"), index=False)
        self._mock_df.to_csv(os.path.join(self._output, "mock_aligned_reads.csv"), index=False)
        # Plot mutation distribution
        self._plot_modification_distribution(tables)
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

        total_edit = tx_indels + mock_indels
        total_reads = self._n_reads_tx + self._n_reads_mock
        no_edit_prior = 1 - edit_prior

        # if no indels in treatment, return True (same results as if run the classifier)
        if tx_indels == 0:
            return False

        # Likelihood function computation
        no_edit_likelihood = hypergeom.pmf(tx_indels, total_reads, total_edit, self._n_reads_tx)
        edit_likelihood = binom.pmf(k=tx_indels, n=total_edit, p=binom_p)

        # Posterior probability computation
        no_edit_post = no_edit_prior * no_edit_likelihood
        edit_post = edit_prior * edit_likelihood

        # In extreme cases both posteriors are smaller than python minimum number (1e-325).
        # The only observed case is when a positions is both with high number of edits and relative high experiment
        # contamination. In this case mark all treatment edits as edit.
        if (edit_post == 0) and (no_edit_post == 0):
            return True

        edit = edit_post > no_edit_post
        return edit

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

    def _plot_modification_distribution(self, tables):
        dist_d = dict()
        amplicon_length = len(tables.amplicon)
        for indel_type in [IndelType.DEL, IndelType.SUB, IndelType.INS]:
            dist_d[indel_type] = np.zeros(amplicon_length + 1, dtype=np.int)

        # Aggregate modifications
        edit_reads = self._tx_df.loc[self._tx_df[IS_EDIT]]
        for row_idx, row in edit_reads.iterrows():
            pos_idx = 0  # position index
            for length, indel_type in InputProcessing.parse_cigar(row[CIGAR]):
                # For a match - continue
                if indel_type == IndelType.MATCH:
                    pos_idx += length
                # Mismatch or deletions
                elif indel_type in [IndelType.DEL, IndelType.SUB]:
                    dist_d[indel_type][pos_idx:pos_idx + length] += row[FREQ]
                    pos_idx += length
                # Insertions
                elif indel_type == IndelType.INS:
                    dist_d[indel_type][pos_idx] += row[FREQ]

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
        ax.plot(positions, dist_d[IndelType.INS], color='#32b165', label="Insertions", linewidth=3, alpha=0.9)
        ax.plot(positions, dist_d[IndelType.SUB], color='b', label="Substitutions", linewidth=3, alpha=0.9)
        ax.plot(positions, dist_d[IndelType.DEL], color='darkviolet', label="Deletions", linewidth=3, alpha=0.9)

        # Create legend, axes limit and labels
        ax.legend()
        ax.set_title("Edited Reads Modification Distribution", weight='bold')
        max_indels = np.array([dist_d[IndelType.INS], dist_d[IndelType.DEL], dist_d[IndelType.SUB]]).max()
        ax.set_ylim(bottom=0, top=max(int(1.1 * max_indels), 10))
        ax.set_xlim(left=0, right=amplicon_length)
        ax.set_xlabel("Position")
        ax.set_ylabel("Indel count")

        with warnings.catch_warnings():
            warnings.simplefilter("ignore", UserWarning)
            fig.savefig(os.path.join(self._output, 'edited_reads_modification_distribution.png'),
                        bbox_inches='tight', dpi=200)
            plt.close(fig)

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
            bar_color = "#39ad48"  # green
        else:
            bar_color = "#db5856"  # red

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
        ax.set_xlim(-0.7, 0.7)

        # Text below each bar plot + y ticks
        ax.set_xticks([0])
        ax.set_xticklabels([self._name])

        ax.set_title(r"Editing Activity with {} % CI".format(self._confidence), weight='bold')

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


