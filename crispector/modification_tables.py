from constants_and_types import ReadsDf, DNASeq, IndelType, ExpType, ModTables, ModTablesP, Path, IsEdit, CIGAR, FREQ, \
    C_TX, C_MOCK, REFERENCE, CS_SHIFT_R, CS_SHIFT_L, ModDist, MOD_TABLE_IDX, INDEL_TYPE, INDEL_LEN, IS_EDIT, \
    POS_IDX_S, POS_IDX_E
from input_processing import InputProcessing
from modification_types import ModificationTypes
import numpy as np
from collections import defaultdict
from matplotlib import pyplot as plt
import matplotlib as mpl
import os
import warnings
import pandas as pd
from utils import Configurator
from copy import deepcopy
from typing import List, Tuple


class ModificationTables:
    """
    A container class for all modification (Indel) tables.
    """
    def __init__(self, tx_reads: ReadsDf, mock_reads: ReadsDf, modifications: ModificationTypes, ref_df_row: pd.Series):
        self._tx_reads = tx_reads
        self._mock_reads = mock_reads
        self._modifications = modifications
        self._amplicon = ref_df_row[REFERENCE]
        self._tables: ModTables = dict()
        self._pointers: ModTablesP = dict()
        self._tx_dist: ModDist = pd.DataFrame() # Tx modification distribution - for visualization only
        self._mock_dist: ModDist = pd.DataFrame()  # Mock modification distribution - for visualization only
        self._create_modification_tables_and_distribution()
        self._n_reads_tx = self._tx_reads[FREQ].sum()
        self._n_reads_mock = self._mock_reads[FREQ].sum()
        self._priors = deepcopy(self._modifications.priors)

        # if cut-site location is ambiguous, use cute-site priors in ambiguous locations as well
        for idx, prior in enumerate(self._priors):
            cut_site_idx = int(len(prior)//2)
            is_ins = self._modifications.types[idx] == IndelType.INS
            if ref_df_row[CS_SHIFT_R]:
                self._priors[idx][cut_site_idx+1] = self._priors[idx][cut_site_idx]
            if ref_df_row[CS_SHIFT_L]:
                # cut-site is left to cut-site index for non insertions, and on cut-site for insertions
                if is_ins:
                    self._priors[idx][cut_site_idx-1] = self._priors[idx][cut_site_idx]
                else:
                    self._priors[idx][cut_site_idx-2] = self._priors[idx][cut_site_idx-1]

    # Getters
    @property
    def priors(self):
        return self._priors

    @property
    def tx_reads(self) -> ReadsDf:
        return self._tx_reads

    @property
    def mock_reads(self) -> ReadsDf:
        return self._mock_reads

    @property
    def tx_dist(self) -> ReadsDf:
        return self._tx_dist

    @property
    def mock_dist(self) -> ReadsDf:
        return self._mock_dist

    @property
    def n_reads_tx(self) -> int:
        return self._n_reads_tx

    @property
    def n_reads_mock(self) -> int:
        return self._n_reads_mock

    @property
    def tables(self) -> ModTables:
        return self._tables

    @property
    def pointers(self) -> ModTablesP:
        return self._pointers

    @property
    def amplicon(self) -> DNASeq:
        return self._amplicon

    def _create_modification_tables_and_distribution(self):
        """
        Function create modification tables and pointers to ReadsDf for all modification types.
        Also - function creates modification distribution data.
        :return: mock table list and mock pointers
        """
        # Create table and pointers
        for idx, indel_type in enumerate(self._modifications.types):
            # Insertions are between positions (so they have an extra item)
            table_size = len(self._amplicon) + 1 if indel_type == IndelType.INS else len(self._amplicon)
            self._tables[idx] = np.zeros((2, table_size), dtype=np.int)
            self._pointers[idx] = defaultdict(list)

        # Fill the tables
        self._tx_dist = self._convert_read_df_to_modifications(self._tx_reads, ExpType.TX)
        self._mock_dist = self._convert_read_df_to_modifications(self._mock_reads, ExpType.MOCK)

    def _convert_read_df_to_modifications(self, read: ReadsDf, exp_type: ExpType) -> ModDist:
        """
        Fill (inplace!) the modification tables from all reads.
        For treatment - fill inplace pointers as well.
        For modification distribution - return the modification distribution
        This function could be split to two different function. The only reason to combine the logic is to reduce
        running time. runtime is O(num_of_reads*num_of_modification)
        :param read: ReadsDf
        :param exp_type: ExpType
        :return: ModDist
        """
        table_row = C_TX if exp_type == ExpType.TX else C_MOCK
        dist_d = defaultdict(list)

        for row_idx, row in read.iterrows():
            pos_idx = 0  # position index
            for length, length_wo_ins, indel_type, mixed_l in InputProcessing.parse_cigar_with_mixed_indels(row[CIGAR]):
                table_idx = self._modifications.find_index(indel_type, length)
                # For a match - continue
                if indel_type == IndelType.MATCH:
                    pos_idx += length_wo_ins
                    continue
                # For a mismatch or deletions - update multi indexes according to length
                elif indel_type in [IndelType.DEL, IndelType.SUB, IndelType.MIXED]:
                    self._tables[table_idx][table_row, pos_idx:pos_idx+length_wo_ins] += row[FREQ]
                    # Update multiple keys
                    if exp_type == ExpType.TX:
                        for pointer_idx in range(pos_idx, pos_idx+length_wo_ins):
                            self._pointers[table_idx][pointer_idx].append(row_idx)
                # For an insertion update single index and don't increase pos_idx
                elif indel_type == IndelType.INS:
                    self._tables[table_idx][table_row, pos_idx] += row[FREQ]
                    if exp_type == ExpType.TX:
                        self._pointers[table_idx][pos_idx].append(row_idx)

                # Add modification to distribution
                indel_l = [(length, indel_type)] if indel_type != IndelType.MIXED else mixed_l
                for dist_length, dist_indel_type in indel_l:
                    dist_d[MOD_TABLE_IDX].append(table_idx)
                    dist_d[INDEL_TYPE].append(dist_indel_type)
                    dist_d[INDEL_LEN].append(dist_length)
                    dist_d[FREQ].append(row[FREQ])
                    dist_d[IS_EDIT].append(False)
                    dist_d[POS_IDX_S].append(pos_idx)
                    end_idx = pos_idx + length_wo_ins if indel_type != IndelType.INS else pos_idx + 1
                    dist_d[POS_IDX_E].append(end_idx)

                # Increase index
                pos_idx += length_wo_ins

        # Convert dict to DataFrame
        dist_df = pd.DataFrame.from_dict(dist_d, orient='columns')
        return dist_df

    def set_tx_dist_is_edit_col(self, edit: IsEdit, table_offset: int, window_size: int):
        """
        Mark all edited modifications from table and algorithm output (edit)
        :param edit: IsEdit
        :param table_offset:
        :param window_size:
        :return:
        """
        edit_size = {IndelType.INS: 2*window_size, IndelType.SUB: 2*window_size-1, IndelType.DEL: 2*window_size-1}
        for row_idx, row in self._tx_dist.iterrows():
            # modification start and end truncated by qualification window
            mod_s = max(row[POS_IDX_S], table_offset)
            mod_e = min(row[POS_IDX_E], table_offset + edit_size[row[INDEL_TYPE]])
            mod_range = range(mod_s - table_offset, mod_e - table_offset)
            if edit[row[MOD_TABLE_IDX]][mod_range].any():
                self._tx_dist.at[row_idx, IS_EDIT] = True

    def plot_tables(self, edit_table: IsEdit, table_offset: int, output: Path, experiment_name: str,
                    table_indexes: List[int], name_suffix = "", figsize: Tuple = (16, 11)):
        """
        Plot all modification tables around cut-site.
        Also display edit events.
        :param edit_table: edit as marked by crispector algorithm
        :param table_offset: reference offset (in bp) to the start of the qualification window
        :param output: output path
        :param experiment_name
        :param table_indexes - list indexes for table plot
        :param name_suffix - prefix for file name
        :return:
        """
        # Set font
        mpl.rcParams.update(mpl.rcParamsDefault)
        mpl.rcParams['font.size'] = 22
        mpl.rcParams['ytick.labelsize'] = 12
        indel_size = 12
        bar_width = 0.4
        mock_color = '#3690c0'  # blue
        tx_color = '#f16a13'  # orange
        edit_color = '#dcdcdc'  # convert_color_to_hex(211, 211, 211) #'drakgrey'# '#e8ebeb' # light grey
        cut_site_color = "red"
        grid_color = 'grey'
        # Create axes
        fig, axes = plt.subplots(nrows=len(table_indexes), ncols=1, figsize=figsize, sharex=True)

        for axes_idx, table_idx in enumerate(table_indexes):

            # Insertions have a slightly different positions marking.
            is_ins = self._modifications.types[table_idx] == IndelType.INS

            # Positions and bars positions
            edit = edit_table[table_idx]
            positions = np.arange(len(edit)) - 0.5 * is_ins
            bar_ind = positions + 0.5
            cut_site = len(positions) / 2

            tx = self._tables[table_idx][C_TX, table_offset:table_offset + len(positions)]
            mock = self._tables[table_idx][C_MOCK, table_offset:table_offset + len(positions)]
            y_max = max(max(tx.max(), mock.max()), 1)

            # Set green color for Edit events
            for pos_idx, pos in enumerate(positions):
                if edit[pos_idx]:
                    axes[axes_idx].axvspan(pos, pos + 1, facecolor=edit_color)

            # Set bp "grid"
            grid_positions = positions if is_ins else np.arange(len(edit) + 1)
            axes[axes_idx].vlines(grid_positions, ymin=0, ymax=y_max, color=grid_color, lw=1)
            axes[axes_idx].spines['bottom'].set_color(grid_color)
            axes[axes_idx].spines['top'].set_color(grid_color)
            axes[axes_idx].spines['right'].set_color(grid_color)
            axes[axes_idx].spines['left'].set_color(grid_color)

            # Set cut-site in red (Insertions cut site is a full bp cell)
            if is_ins:
                axes[axes_idx].vlines([cut_site - 1, cut_site], ymin=0, ymax=y_max, color=cut_site_color, lw=2)
                axes[axes_idx].hlines([0, y_max], xmin=cut_site - 1, xmax=cut_site, color=cut_site_color, lw=4)
            else:
                axes[axes_idx].axvline(cut_site, ymin=0, ymax=y_max, color=cut_site_color, lw=2)

            # Create bar plot
            axes[axes_idx].bar(bar_ind - bar_width / 2, tx, width=bar_width, color=tx_color,
                               label="Tx ({:,} Reads)".format(self._n_reads_tx))
            axes[axes_idx].bar(bar_ind + bar_width / 2, mock, width=bar_width, color=mock_color,
                               label="Mock ({:,} Reads)".format(self._n_reads_mock))

            # Set x, y lim & ticks
            axes[axes_idx].set_xlim(min(positions) - 0.5, max(positions) + 1.5)
            axes[axes_idx].set_xticks([])
            axes[axes_idx].set_yticks([])
            axes[axes_idx].set_ylim(0, y_max)

            # Set y_label - Use text due to alignment issues
            axes[axes_idx].text(x=positions[0] - 5 - 0.5 * (not is_ins), y=y_max / 2,
                                 s=self._modifications.plot_name_at_idx(table_idx), ha="left", va="center")

            # Create legend in the middle of the plot
            if axes_idx == ((len(table_indexes) // 2) - 1):
                axes[axes_idx].bar([0], [0], color=edit_color, label="Edit event", edgecolor='grey')
                axes[axes_idx].plot([], [], color=cut_site_color, label="Cut-Site")
                handles, labels = axes[axes_idx].get_legend_handles_labels()
                order = [1, 2, 3, 0]
                axes[axes_idx].legend([handles[idx] for idx in order], [labels[idx] for idx in order],
                                       bbox_to_anchor=(1.48, 0))

            # Add the reference sequence in the position of the title
            if axes_idx == 0:
                axes[axes_idx].text(0.5, 0.95,
                                    "CRISPECTOR classifier results, by position\n{}".format(experiment_name),
                                    ha="center", va="bottom", weight='bold', size=35,
                                    transform=fig.transFigure, family='serif')
                ref = self._amplicon[table_offset:table_offset + len(positions)]
                for pos_idx, ii in enumerate(bar_ind):
                    axes[axes_idx].text(ii, y_max, ref[pos_idx], ha="center", va="bottom", weight='bold', size=30)
                    axes[axes_idx].text(ii, y_max, ref[pos_idx], ha="center", va="bottom", weight='bold', size=30)
                # red cute-site
                axes[axes_idx].text(cut_site, 1.1 * y_max, "|", ha="center", va="bottom", weight='bold', size=20,
                                    color=cut_site_color)

            # Add bars values (numbers) as text
            for pos_idx, bar_i in enumerate(bar_ind):
                if tx[pos_idx]:
                    axes[axes_idx].text(bar_i - bar_width / 2, 0.05 * y_max, "{:,}".format(tx[pos_idx]), rotation=90,
                                        size=indel_size, ha="center", va="bottom")
                if mock[pos_idx]:
                    axes[axes_idx].text(bar_i + bar_width / 2, 0.05 * y_max, "{:,}".format(mock[pos_idx]), rotation=90,
                                        size=indel_size, ha="center", va="bottom")

        with warnings.catch_warnings():
            warnings.simplefilter("ignore", UserWarning)
            fig.savefig(os.path.join(output, 'classifier_results_by_position_{}.png'.format(name_suffix)),
                        bbox_inches='tight', dpi=200)
            plt.close(fig)

    def dump_tables(self, edit_table: IsEdit, table_offset: int, output: Path):
        """
        Dump all modification tables around cut-site to a csv.
        Also display edit events.
        :param edit_table: edit as marked by crispector algorithm
        :param table_offset: reference offset (in bp) to the start of the qualification window
        :param output: output path
        :return:
        """
        cfg = Configurator.get_cfg()
        win_size = cfg["window_size"]
        ref = self._amplicon[table_offset:table_offset + 2*win_size]
        tables_dict = dict()

        for table_idx in range(self._modifications.size):
            # Gather data
            is_ins = self._modifications.types[table_idx] == IndelType.INS
            edit = edit_table[table_idx]
            pos_len = len(edit)  # length of positions
            tx = self._tables[table_idx][C_TX, table_offset:table_offset + pos_len]
            mock = self._tables[table_idx][C_MOCK, table_offset:table_offset + pos_len]
            cut_site = np.zeros(pos_len, dtype=np.int)
            cut_site[(pos_len // 2) - 1 * (not is_ins):(pos_len // 2) + 1] = 1

            # Insert to dictionary
            table_name = self._modifications.name_at_idx(table_idx)
            tables_dict[table_name + " Tx"] = np.insert(tx, 0, self._n_reads_tx)
            tables_dict[table_name + " Mock"] = np.insert(mock, 0, self._n_reads_mock)
            tables_dict[table_name + " Is edited"] = [None] + edit
            tables_dict[table_name + " Is cut site"] = [None] + list(cut_site)

            tables_df = pd.DataFrame.from_dict(tables_dict, orient='index')
            new_cols = ['number of reads'] + [c for c in ref] + ['insertion at the end']

        tables_df.set_axis(new_cols, axis=1, inplace=True)
        tables_df.to_csv(os.path.join(output, "modification_tables.csv"))
