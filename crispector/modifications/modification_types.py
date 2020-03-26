from crispector.utils.exceptions import PriorPositionHasWrongLength
from crispector.utils.constants_and_types import IndelType
from typing import List
from collections import defaultdict
from crispector.utils.configurator import Configurator



class ModificationTypes:
    """
    A container class for all modification (Indel) types
    """
    def __init__(self, indel_type: List[IndelType], min_l: List[int], max_l: List[int], priors: List[List[float]],
                 max_indel_size):
        """
        :param indel_type: IndelType
        :param min_l: indel_type length min
        :param max_l: indel_type length min
        :param priors: prior for each modification type
        :param max_indel_size: max number of indel (from config file)
        """
        self._type = indel_type
        self._range = []
        for indel_min, indel_max in zip(min_l, max_l):
            self._range.append(range(indel_min, indel_max+1))
        self._priors = priors
        self._max_indel_size = max_indel_size
        self._type_d = defaultdict(list)
        for idx, (indel_type, indel_range) in enumerate(zip(self._type, self._range)):
            self._type_d[indel_type].append((indel_range, idx))

    @property
    def types(self):
        return self._type

    @property
    def ranges(self):
        return self._range

    @property
    def priors(self):
        return self._priors

    @property
    def size(self):
        return len(self._type)

    def find_index(self, indel_type: IndelType, length: int) -> int:
        """
        Function returns the index of the modification type
        :param indel_type:
        :param length:
        :return: modification type index
        """
        if indel_type == IndelType.MATCH:
            return -1

        index = 0
        for indel_range, indel_idx in self._type_d[indel_type]:
            if length in indel_range:
                index = indel_idx
                break
        return index

    @classmethod
    def init_from_cfg(cls, enable_substitutions: bool):
        """
        Parse the config file to ModificationTypes class
        :param enable_substitutions: Flag
        :return: ModificationTypes.
        """
        cfg = Configurator.get_cfg()["NHEJ_inference"]
        indel_types_cfg = cfg["IndelTypes"]
        window_size = cfg["window_size"]
        m_list = []
        m_min = []
        m_max = []
        m_prior = []
        for indel_type in IndelType:
            if indel_type == IndelType.MATCH:
                continue
            for key, indel_dict in indel_types_cfg[indel_type.name].items():
                m_list.append(indel_type)
                m_min.append(indel_dict["min"])
                m_max.append(indel_dict["max"])

                # Extract priors
                pos_prior = indel_dict["pos_prior"]
                # If enable_substitutions is false, change priors to zero
                if not enable_substitutions and (indel_type == IndelType.SUB):
                    pos_prior = len(pos_prior)*[0.0]

                prior_expected_len = (2*window_size + 1) if indel_type == IndelType.INS else 2*window_size
                if len(pos_prior) != prior_expected_len:
                    raise PriorPositionHasWrongLength(prior_expected_len, len(pos_prior), indel_type, m_min[-1],
                                                      m_max[-1])
                m_prior.append(pos_prior)

        return cls(m_list, m_min, m_max, m_prior, cfg["max_indel_size"])

    def name_at_idx(self, idx: int) -> str:
        """
        Return modification name by index
        :param idx:
        :return:
        """
        max_range = "INF" if self._max_indel_size == max(self._range[idx]) else max(self._range[idx])
        return "{} ({}-{})".format(self._type[idx].name, min(self._range[idx]), max_range)

    def plot_name_at_idx(self, idx: int) -> str:
        """
        Return modification plot name (identical size) by index
        :param idx:
        :return:
        """
        if self._max_indel_size == max(self._range[idx]):
            return r"{} of len $\geq {}$".format(self._type[idx].plot_name, min(self._range[idx]))
        elif min(self._range[idx]) == max(self._range[idx]):
            return r"{} of len ${}$".format(self._type[idx].plot_name, min(self._range[idx]))
        else:
            return r"{} of len ${}:{}$".format(self._type[idx].plot_name, min(self._range[idx]), max(self._range[idx]))



