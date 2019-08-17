from crispector_constants import CIGAR, FREQ
from crispector_types import ReadsDf, DNASeq, IndelType, ExpType, ModTables, ModTablesP
from input_process_utils import InputProcess
from modification_types import ModificationTypes
import numpy as np
from collections import defaultdict


class ModificationTables:
    """
    A container class for all modification (Indel) tables.
    """
    def __init__(self, tx_reads: ReadsDf, mock_reads: ReadsDf, modifications: ModificationTypes, amplicon: DNASeq):
        self._tx_reads = tx_reads
        self._mock_reads = mock_reads  # TODO add setters
        self._modifications = modifications
        self._amplicon = amplicon
        self._tables: ModTables = dict()
        self._pointers: ModTablesP = dict()
        self._create_modification_tables()

    def _create_modification_tables(self):
        """
        Function create modification tables and pointers to ReadsDf for all modification types.
        :return: mock table list and mock pointers
        """
        # Create table and pointers
        for idx, indel_type in enumerate(self._modifications.types):
            # Insertions are between position (so they have an extra item)
            table_size = len(self._amplicon) + 1 if indel_type == IndelType.INS else len(self._amplicon)
            self._tables[idx] = np.zeros((2, table_size))
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
        table_row = 0 if exp_type == ExpType.TX else 1
        for row_idx, row in read.iterrows():
            pos_idx = 0  # position index
            for length, indel_type in InputProcess.parse_cigar(row[CIGAR]):
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
