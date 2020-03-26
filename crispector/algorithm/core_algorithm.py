from crispector.utils.constants_and_types import IsEdit, IndelType, AlgResult, Pr, FREQ, IS_EDIT, C_TX, C_MOCK, TX_READ_NUM, \
    MOCK_READ_NUM, TX_EDIT, EDIT_PERCENT, CI_LOW, CI_HIGH, ReadsDf
from crispector.utils.exceptions import ClassificationFailed
from crispector.utils.logger import LoggerWrapper
from crispector.utils.configurator import Configurator
from crispector.modifications.modification_tables import ModificationTables
from crispector.modifications.modification_types import ModificationTypes
from typing import List, Tuple, Dict
from scipy.stats import norm, binom, hypergeom
import numpy as np


class CoreAlgorithm:
    """
    Core algorithm for CRISPECTOR
    """

    def __init__(self, cut_site: int, modification: ModificationTypes, binom_p_l: List[Pr], confidence: Pr, on_target: bool):
        """
        :param cut_site:
        :param modification:
        :param binom_p_l: binom_p list
        :param confidence: confidence interval
        :param on_target: is on target flag
        """
        self._cut_site = cut_site
        self._modifications = modification
        self._binom_p_l = binom_p_l
        self._confidence = confidence
        self._cfg = Configurator.get_cfg()
        self._win_size = self._cfg["NHEJ_inference"]["window_size"]
        self._tables_offset = self._cut_site - self._win_size  # position offset
        self._n_reads_tx = 0
        self._n_reads_mock = 0
        self._edit: IsEdit = dict()
        self._tx_df: ReadsDf = None
        self._mock_df: ReadsDf = None
        self._logger = LoggerWrapper.get_logger()
        self._on_target = on_target

    def evaluate(self, tables: ModificationTables) -> Dict:
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

        # Mark all edited modifications as an edit in modification distribution
        tables.set_tx_dist_is_edit_col(self._edit, self._tables_offset, self._win_size)

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
            total_indels = tx_indels + mock_indels

            # Likelihood function computation
            no_edit_likelihood = hypergeom.pmf(tx_indels, total_reads, total_indels, self._n_reads_tx)
            edit_likelihood = binom.pmf(k=tx_indels, n=total_indels, p=binom_p)

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

    # Getters
    @property
    def edited(self) -> IsEdit:
        return self._edit

    @property
    def tables_offset(self) -> int:
        return self._tables_offset

    @property
    def cut_site(self) -> int:
        return self._cut_site

    @property
    def win_size(self) -> int:
        return self._win_size

    @property
    def confidence(self) -> int:
        return self._confidence

    @property
    def is_on_target(self) -> bool:
        return self._on_target
