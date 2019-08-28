from constants import CIGAR, FREQ, C_TX, C_MOCK
from enum_types import ReadsDf, DNASeq, IndelType, ExpType, ModTables, ModTablesP, Path, IsEdit
from input_processing import InputProcessing
from modification_types import ModificationTypes
import numpy as np
from collections import defaultdict
from matplotlib import pyplot as plt  # TODO - add to package requirements
import matplotlib as mpl
import os
import warnings


class ModificationTables:
    """
    A container class for all modification (Indel) tables.
    """
    def __init__(self, tx_reads: ReadsDf, mock_reads: ReadsDf, modifications: ModificationTypes, amplicon: DNASeq):
        self._tx_reads = tx_reads
        self._mock_reads = mock_reads
        self._modifications = modifications
        self._amplicon = amplicon
        self._tables: ModTables = dict()
        self._pointers: ModTablesP = dict()
        self._create_modification_tables()
        self._n_reads_tx = self._tx_reads[FREQ].sum()
        self._n_reads_mock = self._mock_reads[FREQ].sum()

    # Getters
    @property
    def tx_reads(self) -> ReadsDf:
        return self._tx_reads

    @property
    def mock_reads(self) -> ReadsDf:
        return self._mock_reads

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

    def _create_modification_tables(self):
        """
        Function create modification tables and pointers to ReadsDf for all modification types.
        :return: mock table list and mock pointers
        """
        # Create table and pointers
        for idx, indel_type in enumerate(self._modifications.types):
            # Insertions are between positions (so they have an extra item)
            table_size = len(self._amplicon) + 1 if indel_type == IndelType.INS else len(self._amplicon)
            self._tables[idx] = np.zeros((2, table_size), dtype=np.int)
            self._pointers[idx] = defaultdict(list)

        # Fill the tables
        self._fill_modification_table_from_read_df(self._tx_reads, ExpType.TX)
        self._fill_modification_table_from_read_df(self._mock_reads, ExpType.MOCK)

    def _fill_modification_table_from_read_df(self, read: ReadsDf, exp_type: ExpType):
        """
        Fill the modification tables from all reads.
        For treatment - keep pointers as well.
        :param read: ReadsDf
        :param exp_type: ExpType
        :return:
        """
        table_row = C_TX if exp_type == ExpType.TX else C_MOCK
        for row_idx, row in read.iterrows():
            pos_idx = 0  # position index
            for length, indel_type in InputProcessing.parse_cigar(row[CIGAR]):
                table_idx = self._modifications.find_index(indel_type, length)
                # For a match - continue
                if indel_type == IndelType.MATCH:
                    pos_idx += length
                # For a mismatch or deletions - update multi indexes according to length
                elif indel_type in [IndelType.DEL, IndelType.SUB]:
                    self._tables[table_idx][table_row, pos_idx:pos_idx+length] += row[FREQ]
                    # update multiple keys
                    if exp_type == ExpType.TX:
                        for pointer_idx in range(pos_idx, pos_idx+length):
                            self._pointers[table_idx][pointer_idx].append(row_idx)
                    pos_idx += length
                # For an insertion update single index and don't increase pos_idx
                elif indel_type == IndelType.INS:
                    self._tables[table_idx][table_row, pos_idx] += row[FREQ]
                    if exp_type == ExpType.TX:
                        self._pointers[table_idx][pos_idx].append(row_idx)

    def plot_tables(self, edit_table: IsEdit, table_offset: int, output: Path):
        """
        Plot all modification tables around cut-site.
        Also display edit events.
        :param edit_table: edit as marked by crispector algorithm
        :param table_offset: reference offset (in bp) to the start of the qualification window
        :param output: output path
        :return:
        """
        # Set font
        mpl.rcParams['font.size'] = 22
        mpl.rcParams['ytick.labelsize'] = 12
        indel_size = 12
        bar_width = 0.4

        # Create axes
        fig, axes = plt.subplots(nrows=self._modifications.size, ncols=1, figsize=(12, 10), sharex=True)

        for table_idx in range(self._modifications.size):

            # Insertions have a slightly different positions marking.
            is_ins = self._modifications.types[table_idx] == IndelType.INS

            # Positions and bars positions
            edit = edit_table[table_idx]
            positions = np.arange(len(edit)) - 0.5*is_ins
            bar_ind = positions + 0.5
            cut_site = len(positions) / 2

            tx = self._tables[table_idx][C_TX, table_offset:table_offset+len(positions)]
            mock = self._tables[table_idx][C_MOCK, table_offset:table_offset+len(positions)]
            y_max = max(max(tx.max(), mock.max()), 1)

            # Set green color for Edit events
            for pos_idx, pos in enumerate(positions):
                if edit[pos_idx]:
                    axes[table_idx].axvspan(pos, pos + 1, facecolor='green', alpha=0.3)

            # Set bp "grid"
            grid_positions = positions if is_ins else np.arange(len(edit) + 1)
            axes[table_idx].vlines(grid_positions, ymin=0, ymax=y_max, color='grey', lw=1, alpha=0.3)

            # Set cut-site in red (Insertions cut site is a full bp cell)
            if is_ins:
                axes[table_idx].vlines([cut_site-1, cut_site], ymin=0, ymax=y_max, color='red', lw=2)
                axes[table_idx].hlines([0, y_max], xmin=cut_site-1, xmax=cut_site, color='red', lw=4)
            else:
                axes[table_idx].axvline(cut_site, ymin=0, ymax=y_max, color='red', lw=2)

            # Create bar plot
            axes[table_idx].bar(bar_ind - bar_width/2, tx, width=bar_width, color='darkorange',
                                label="Tx ({} Reads)".format(self._n_reads_tx))
            axes[table_idx].bar(bar_ind + bar_width/2, mock, width=bar_width, color='darkgrey',
                                label="Mock ({} Reads)".format(self._n_reads_mock))

            # Set x, y lim & ticks
            axes[table_idx].set_xlim(min(positions) - 0.5, max(positions) + 1.5)
            axes[table_idx].set_xticks([])
            axes[table_idx].set_yticks([0, y_max])
            axes[table_idx].set_yticklabels([0, y_max])
            axes[table_idx].set_ylim(0, y_max)

            # Set y_label - Use text due to alignment issues
            axes[table_idx].text(x=positions[0] - 3.5 - 0.5*(not is_ins), y=y_max/2,
                                 s=self._modifications.plot_name_at_idx(table_idx), ha="left", va="center")

            # Create legend in the middle of the plot
            if table_idx == 2:
                axes[table_idx].bar([0], [0], color='green', alpha=0.3, label="Edit event")
                axes[table_idx].plot([], [], color='r', label="Cut-Site")
                handles, labels = axes[table_idx].get_legend_handles_labels()
                order = [1, 2, 3, 0]
                axes[table_idx].legend([handles[idx] for idx in order], [labels[idx] for idx in order],
                                       bbox_to_anchor=(1.5, 0))

            # Add the reference sequence in the position of the title
            if table_idx == 0:
                ref = self._amplicon[table_offset:table_offset+len(positions)]
                for pos_idx, ii in enumerate(bar_ind):
                    axes[table_idx].text(ii, y_max, ref[pos_idx], ha="center", va="bottom")

            # Add bars values (numbers) as text
            for pos_idx, bar_i in enumerate(bar_ind):
                if tx[pos_idx]:
                    axes[table_idx].text(bar_i - bar_width / 2, 0.05 * y_max, tx[pos_idx], rotation=90,
                                         size=indel_size, ha="center", va="bottom")
                if mock[pos_idx]:
                    axes[table_idx].text(bar_i + bar_width / 2,  0.05 * y_max, mock[pos_idx], rotation=90,
                                         size=indel_size, ha="center", va="bottom")

        with warnings.catch_warnings():
            warnings.simplefilter("ignore", UserWarning)
            fig.savefig(os.path.join(output, 'modification_tables.png'), bbox_inches='tight', dpi=300)
            plt.close(fig)
