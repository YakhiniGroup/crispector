from crispector_types import IndelType
from typing import List, Dict
from collections import defaultdict


# TODO - maybe need to add iterator for this class

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
        :param priors: prior for eac modification type
        :param max_indel_size: max number of indel (from config file)
        """
        self._type = indel_type
        self._range = []
        for indel_min, indel_max in zip(min_l, max_l):
            self._range.append(range(indel_min, indel_max+1))
        self._priors = priors
        self.max_indel_size = max_indel_size
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

    @priors.setter
    def priors(self, value: List[List]):
        self._priors = value

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
    def init_from_cfg(cls, cfg: Dict):
        """
        Parse the config file to ModificationTypes class
        :param cfg: config (Dict)
        :return: ModificationTypes.
        """
        indel_types_cfg = cfg["IndelTypes"]
        m_list = []
        m_min = []
        m_max = []
        m_prior = []
        for indel_type in IndelType:
            if indel_type == IndelType.MATCH:
                continue
            for key, indel_dict in indel_types_cfg[indel_type.name()].items():
                m_list.append(indel_type)
                m_min.append(indel_dict["min"])
                m_max.append(indel_dict["max"])
                m_prior.append(indel_dict["prior"])

        return cls(m_list, m_min, m_max, m_prior, cfg["max_indel_size"])


