from crispector_types import ExpType

welcome_msg = "\n\
 CCCCC  RRRRRR  IIIII  SSSSS  PPPPPP  EEEEEEE  CCCCC  TTTTTTT  OOOOO  RRRRRR\n\
CC    C RR   RR  III  SS      PP   PP EE      CC    C   TTT   OO   OO RR   RR\n\
CC      RRRRRR   III   SSSSS  PPPPPP  EEEEE   CC        TTT   OO   OO RRRRRR\n\
CC    C RR  RR   III       SS PP      EE      CC    C   TTT   OO   OO RR  RR\n\
 CCCCC  RR   RR IIIII  SSSSS  PP      EEEEEEE  CCCCC    TTT    OOOO0  RR   RR\n\
\n\
              CRISPR/CAS9 Off-Target Analysis From NGS Data\n\
                            version={}\n".format("0.1.0")

""
# TODO -  Need to auto import from somewhere -search in online repositories

# fastp constants
FASTP_DIR = dict()
FASTP_DIR[ExpType.TX] = "treatment_fastp"
FASTP_DIR[ExpType.MOCK] = "mock_fastp"

# AmpliconDf constants
SITE_NAME, REFERENCE, SGRNA, ON_TARGET, CUT_SITE = 'site_name', 'reference', 'sgRNA', 'on_target', 'cut-site'

# ReadDf constants
READ = "read"
ALIGNMENT_W_INS = "alignment_w_ins"
ALIGNMENT_W_DEL = "alignment_w_del"
CIGAR = 'cigar_path'
FREQ = "frequency"
IS_EDIT = "is_edited"

# Algorithm constants
C_TX = 0
C_MOCK = 1

# General constants
COMPLEMENT = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
NONE = "NONE"

# AlgResult columns
TX_READ_NUM = "treatment_number_of_reads"
MOCK_READ_NUM = "mock_number_of_reads"
TX_EDIT = "edited_reads"
EDIT_PERCENT = "editing_activity"
CI_LOW = "CI_low"
CI_HIGH = "CI_high"
