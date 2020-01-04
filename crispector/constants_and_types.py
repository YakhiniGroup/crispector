from typing import Dict, List, DefaultDict, Tuple
from enum import Enum
import pandas
import numpy as np


version = "0.1.0"
welcome_msg = "\n\
 CCCCC  RRRRRR  IIIII  SSSSS  PPPPPP  EEEEEEE  CCCCC  TTTTTTT  OOOOO  RRRRRR\n\
CC    C RR   RR  III  SS      PP   PP EE      CC    C   TTT   OO   OO RR   RR\n\
CC      RRRRRR   III   SSSSS  PPPPPP  EEEEE   CC        TTT   OO   OO RRRRRR\n\
CC    C RR  RR   III       SS PP      EE      CC    C   TTT   OO   OO RR  RR\n\
 CCCCC  RR   RR IIIII  SSSSS  PP      EEEEEEE  CCCCC    TTT    OOOOO  RR   RR\n\
\n\
              CRISPR/CAS9 Off-Target Analysis From NGS Data\n\
                            version={}\n".format(version)
# TODO -  Need to read version from setup tools just like CRISPRESSO2


# Enums


class ExpType(Enum):
    TX = 0
    MOCK = 1

    def name(self):
        if self._name_ == "TX":
            return "treatment"
        else:
            return "mock"


class IndelType(Enum):
    DEL = 0
    INS = 1
    MIXED = 2
    SUB = 3
    MATCH = -1

    @property
    def name(self):
        if self._name_ == "DEL":
            return "Deletions"
        elif self._name_ == "INS":
            return "Insertions"
        elif self._name_ == "SUB":
            return "Substitutions"
        elif self._name_ == "MIXED":
            return "Mixed"
        else:
            return "Match"
    @property
    def plot_name(self):
        if self._name_ == "DEL":
            return "Del"
        elif self._name_ == "INS":
            return "Ins"
        elif self._name_ == "SUB":
            return "Sub"
        elif self._name_ == "MIXED":
            return "Mix"
        else:
            return "Match"

    @property
    def color(self):
        if self._name_ == "DEL":
            return 'darkviolet'
        elif self._name_ == "INS":
            return '#32b165'
        elif self._name_ == "SUB":
            return 'b'
        else:
            return "black"

    @classmethod
    def from_cigar(cls, indel_type: str):
        if indel_type == CIGAR_I:
            return cls.INS
        elif indel_type == CIGAR_D:
            return cls.DEL
        elif indel_type == CIGAR_S:
            return cls.SUB
        else:
            return cls.MATCH

# Pandas DataFrames
# pandas data frame with 5 columns: SITE_NAME, REFERENCE, SGRNA, ON_TARGET, CUT_SITE
AmpliconDf = pandas.DataFrame

# pandas data frame with all reads columns :READ, ALIGNMENT_W_INS, ALIGNMENT_W_DEL, CIGAR, CUT_SITE...
ReadsDf = pandas.DataFrame

# pandas data frame with all reads columns :READ, TRANS_NAME, REFERENCE, ALIGNMENT_W_INS, ALIGNMENT_W_DEL, CUT_SITE,
# R_SITE, L_SITE
TransDf = pandas.DataFrame

# pandas data frame with columns: SITE_A, SITE_B, TX_TRANS_READ, MOCK_TRANS_READ, TRANS_PVAL, TRANS_FDR
TransResultDf = pandas.DataFrame

# Others
# A dictionary with key=SITE_NAME and value=ReadsDf - Dict for all reads in the experiments
ReadsDict = Dict[str, ReadsDf]

# numpy array size (2 * reference_sequence length)
ModTable = np.ndarray
# A dict of keys modification table index and value IndelTable
ModTables = Dict[int, ModTable]

# pandas data frame with all modification. columns - INDEL_TYPE, INDEL_LEN, MOD_TABLE_IDX, FREQ, IS_EDIT, POS_IDX_R
# R_SITE, L_SITE
ModDist = pandas.DataFrame

# Dictionary of pointers. The key is position index and the value is row_idx for  ReadsDf
ModTableP = DefaultDict[int, List]
# A dict of keys modification table index and value ModTableP
ModTablesP = Dict[int, ModTableP]

# A dict of keys modification table index and value List of booleans
IsEdit = Dict[int, np.ndarray]

DNASeq = str
Path = str
Pr = float
CigarPath = str

# AlgResult - dictionary with key name_of_property (e.g. 'CI_high') and their value
AlgResult = Dict[str, Dict]
AlgResultDf = pandas.DataFrame

# fastp constants
FASTP_DIR = dict()
FASTP_DIR[ExpType.TX] = "treatment_fastp"
FASTP_DIR[ExpType.MOCK] = "mock_fastp"

# AmpliconDf constants
SITE_NAME, REFERENCE, SGRNA, ON_TARGET, CUT_SITE = 'site_name', 'reference', 'sgRNA', 'Site type', 'cut-site'
F_PRIMER = 'forward_primer'
R_PRIMER = 'reverse_primer'
MAX_SCORE = 'max_score'
CS_SHIFT_R = 'possible_cut_site_shift_to_right'
CS_SHIFT_L = 'possible_cut_site_shift_to_left'
SGRNA_REVERSED = 'sgRNA_reversed'

# ReadDf constants
READ = "read"
ALIGNMENT_W_INS = "alignment_w_ins"
ALIGNMENT_W_DEL = "alignment_w_del"
CIGAR = 'cigar_path'
CIGAR_LEN = 'cigar_length'
ALIGN_SCORE = 'alignment_score'
NORM_SCORE = "normalized_alignment_score"
FREQ = "frequency"
IS_EDIT = "is_edited"
INS_LEN = "inserted_bases_length"
INS_POS = "insertion_position"
INS_BASE = "inserted_bases"
DEL_LEN = "deleted_bases_length"
DEL_START = "deletion_start_position"
DEL_END = "deletion_end_position"
DEL_BASE = "deleted_bases"
SUB_CNT = "substitution_count"
SUB_POS = "substitution_position"
SUB_BASE = "substitution_bases"
REVERSED = "reversed_reads"
L_SITE = "left_site_name"
L_REV = "left_site_reversed"
R_SITE = "right_site_name"
R_REV = "right_site_reversed"
R_READ = "right_primer_read"
L_READ = "left_primer_read"
ALIGN_CUT_SITE = "alignment_cut_site"
INDEL_COLS = [ALIGN_CUT_SITE, DEL_LEN, DEL_START, DEL_END, DEL_BASE, INS_LEN, INS_POS, INS_BASE, SUB_CNT, SUB_POS, SUB_BASE]

# TransDf constants
TRANS_NAME = "translocation_name"
IS_TRANS = "is_translocation"

# ModDist constants
INDEL_TYPE = 'IndelType'
INDEL_LEN = 'indel_length'
MOD_TABLE_IDX = 'modification_table_index'
POS_IDX_S = 'position_index_start'
POS_IDX_E = 'position_index_end'

# Binomial probability estimation
BINOM_AVERAGE_P = 0.9
BINOM_PERCENTILE = 0.95

# General constants
C_TX = 0
C_MOCK = 1
PRIMER_LEN = 20
COMPLEMENT = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
BAD_AMPLICON_THRESHOLD = 500
CIGAR_LEN_THRESHOLD = 8 # This threshold is applied only if alignment score is low
MIN_PRIMER_DIMER_THRESH = 10
READ_LEN_SIDE = 20
OUTPUT_DIR = "crispector_output"

# AlgResult columns
TX_READ_NUM = "Treatment number of reads"
MOCK_READ_NUM = "Mock number of reads"
TX_EDIT = "Edited reads"
EDIT_PERCENT = "Editing Activity"
CI_LOW = "CI_low"
CI_HIGH = "CI_high"
SUMMARY_RESULTS_TITLES = [SITE_NAME, ON_TARGET, MOCK_READ_NUM, TX_READ_NUM, TX_EDIT, EDIT_PERCENT, CI_LOW, CI_HIGH]

#TransResultDf columns
SITE_A = "site_A"
SITE_B = "site_B"
TX_TRANS_READ = "treatment_reads"
MOCK_TRANS_READ = "mock_reads"
TRANS_PVAL = "p_value"
TRANS_FDR = "FDR"

# Cigar path constants
CIGAR_D, CIGAR_I, CIGAR_S, CIGAR_M = "D", "I", "X", "="

# AlignedIndel - A Tuple of indel_type, length, and position in alignment coordinates.
# Only used by input processing module.
AlignedIndel = Tuple[IndelType, int, int]

# plots colors
OFF_TARGET_COLOR = "#db5856"
ON_TARGET_COLOR = "#39ad48"

# html constants
DISCARDED_SITES = "discarded_sites"
EDITING_ACTIVITY = "editing_activity"
W = "width"
H = "height"
PLOT_PATH = "plot_path"
PDF_PATH =  "pdf_path"
TITLE = "title"
PAGE_TITLE = "page_title"
READING_STATS = "reading_statistics"
MAPPING_STATS = "mapping_stats"
MAPPING_PER_SITE = "mapping_per_site"
FASTP_TX_PATH = "fastp_tx_path"
FASTP_MOCK_PATH = "fastp_mock_path"
RESULT_TABLE = "result_table"
TAB_DATA = "tabular_data"
HTML_SITE_NAMES = "site_names"
LOG_PATH = "log_path"