from typing import Dict, List, DefaultDict
from enum import Enum
import pandas
import numpy as np


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
    SUB = 2
    MATCH = -1

    @property
    def name(self):
        if self._name_ == "DEL":
            return "Deletions"
        elif self._name_ == "INS":
            return "Insertions"
        elif self._name_ == "SUB":
            return "Substitutions"
        else:
            return "Match"

    @classmethod
    def from_cigar(cls, indel_type: str):
        if indel_type == "I":
            return cls.INS
        elif indel_type == "D":
            return cls.DEL
        elif indel_type == "X":
            return cls.SUB
        else:
            return cls.MATCH


# Pandas DataFrames

# pandas data frame with 5 columns: SITE_NAME, REFERENCE, SGRNA, ON_TARGET, CUT_SITE
AmpliconDf = pandas.DataFrame

# TODO - Add all columns
# pandas data frame with XXX columns READ, ALIGNMENT_W_INS, ALIGNMENT_W_DEL, CIGAR, CUT_SITE,
ReadsDf = pandas.DataFrame

# Others

# A dictionary with key=SITE_NAME and value=ReadsDf - Dict for all reads in the experiments
ReadsDict = Dict[str, ReadsDf]

# numpy array size (2 * reference_sequence length)
ModTable = np.ndarray
# A dict of keys modification table index and value IndelTable
ModTables = Dict[int, ModTable]

# Dictionary of pointers. The key is position index and the value is row_idx for  ReadsDf
ModTableP = DefaultDict[int, List]
# A dict of keys modification table index and value ModTableP
ModTablesP = Dict[int, ModTableP]

# A dict of keys modification table index and value List of booleans
IsEdit = Dict[int, List[bool]]

DNASeq = str
Path = str
Pr = float

# TODO - change description
# AlgResult - dictionary with key name_of_property (e.g. 'CI_high') and their value
AlgResult = Dict[str, float]




