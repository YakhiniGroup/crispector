# from docs.conf import version TODO - can't import from docs. crispector is not a module
from markdown import __version__

welcome_msg = "CRISPECTOR {}\n ADD SOMETHING COOL".format(__version__)


# fastp
TX_FASTP_DIR = "treatment_fastp"
MOCK_FASTP_DIR = "mock_fastp"


SITE_NAME, REFERENCE, SGRNA, ON_TARGET, CUT_SITE = 'site_name', 'reference', 'sgRNA', 'on_target', 'cut-site'

COMPLEMENT = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
