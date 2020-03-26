from crispector.input_processing.utils import parse_cigar_with_mixed_indels
from crispector.utils.constants_and_types import ReadsDf, DNASeq, IndelType, ExpType, ModTables, ModTablesP, IsEdit, CIGAR, FREQ, \
    C_TX, C_MOCK, REFERENCE, CS_SHIFT_R, CS_SHIFT_L, ModDist, MOD_TABLE_IDX, INDEL_TYPE, INDEL_LEN, IS_EDIT, \
    POS_IDX_S, POS_IDX_E
from crispector.modifications.modification_types import ModificationTypes
import numpy as np
from collections import defaultdict
import pandas as pd
from copy import deepcopy


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
            for length, length_wo_ins, indel_type, mixed_l in parse_cigar_with_mixed_indels(row[CIGAR]):
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
