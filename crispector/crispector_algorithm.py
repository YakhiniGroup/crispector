from crispector_constants import FREQ, TX_POS, MOCK_POS, IS_EDIT, TX_READ_NUM, MOCK_READ_NUM, TX_EDIT, EDIT_PERCENT, \
    CI_LOW, CI_HIGH
from crispector_types import IsEdit, IndelType, AlgResult, Pr, Path
from crispector_utils import Configurator
from modification_tables import ModificationTables
from modification_types import ModificationTypes
from typing import List, Tuple
from scipy.stats import norm, binom, hypergeom  # TODO - add to package requirements
import numpy as np


class CrispectorAlgorithm:
    """
    Core algorithm for CRISPECTOR
    """

    def __init__(self, site_name: str, cut_site: int, modification: ModificationTypes, binom_p_l: List[Pr],
                 confidence: Pr, output: Path):
        """
        :param site_name:
        :param cut_site:
        :param modification:
        :param binom_p_l: binom_p list
        :param confidence: confidence interval
        """
        self._name = site_name
        self._cut_site = cut_site
        self._modification = modification
        self._binom_p_l = binom_p_l
        self._confidence = confidence
        self._cfg = Configurator.get_cfg()
        self._win_size = self._cfg["window_size"]
        self._tx_reads_num = 0
        self._mock_reads_num = 0
        self._edit: IsEdit = dict()
        self._output = output
        self._tx_df = None
        self._mock_df = None

    def evaluate(self, tables: ModificationTables) -> AlgResult:
        """
        Evaluate editing activity from ModificationTables, and:
        TODO - Add description
        :param tables:
        :return: Dict with results
        """
        # set algorithm properties
        self._tx_df = tables.tx_reads
        self._mock_df = tables.mock_reads
        self._tx_reads_num = self._tx_df[FREQ].sum()
        self._mock_reads_num = self._mock_df[FREQ].sum()
        edited_indexes = []

        # Run evaluation on each modification type
        for table_idx, (indel_type, prior) in enumerate(zip(self._modification.types, self._modification.priors)):
            table = tables.tables[table_idx]
            pointers = tables.pointers[table_idx]
            binom_p = self._binom_p_l[table_idx]
            offset = self._cut_site - self._win_size  # position offset
            eval_size = (2*self._win_size+1) if indel_type == IndelType.INS else 2*self._win_size
            self._edit[table_idx] = eval_size*[False]

            # Run evaluation on every position
            for pos_idx in range(eval_size):
                tx_indels, mock_indels = table[[TX_POS, MOCK_POS], offset+pos_idx]
                is_edit = self._classify_position(tx_indels, mock_indels, prior[pos_idx], binom_p)
                self._edit[table_idx][pos_idx] = is_edit
                if is_edit:
                    edited_indexes += pointers[offset+pos_idx]

        # Store is_edit in modification tables
        # TODO - Add it + call for excel + plot dump plot_modification_table + modification_table_to_excel
        # TODO - Add bar plot for editing activity
        # TODO - Mutation histogram plot using cigar path.
        # TODO - Dump both tx and mock read  tables, with n_deleted, mutated and so on, only in the 10 windows base.
        # maybe do it with 40+45 in range(130,150).
        # TODO - Add most frequent allels plot . Tx up (10) and Mock down (10).

        # Compute editing activity.
        result_dict = self._compute_editing_activity(edited_indexes)

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
        total_reads = self._tx_reads_num + self._mock_reads_num
        no_edit_prior = 1 - edit_prior

        # if no indels in treatment, return True (same results as if run the classifier)
        if tx_indels == 0:
            return False

        # Likelihood function computation
        no_edit_likelihood = hypergeom.pmf(tx_indels, total_reads, total_edit, self._tx_reads_num)
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
        editing_activity = edited_reads / self._tx_reads_num

        # Compute CI
        CI_low, CI_high = self._compute_confidence_interval(editing_activity)

        # Prepare return dict
        result_d: AlgResult = dict()
        result_d[TX_READ_NUM] = self._tx_reads_num
        result_d[MOCK_READ_NUM] = self._mock_reads_num
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
        half_len_CI = confidence_inv * np.sqrt((editing_activity*(1-editing_activity))/self._tx_reads_num)
        return max(0, editing_activity - half_len_CI), editing_activity + half_len_CI

