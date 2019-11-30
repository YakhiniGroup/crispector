
class FastpRunTimeError(Exception):
    pass


class Bowtie2BuildRunTimeError(Exception):
    pass


class Bowtie2RunTimeError(Exception):
    pass


class NoneValuesInAmpliconsCSV(Exception):
    pass


class SgRNANotInReferenceSequence(Exception):
    def __init__(self, site_name):
        self.site_name = site_name


class CantOpenDemultiplexedSamFile(Exception):
    def __init__(self, path):
        self.path = path


class CantOpenMergedFastqFile(Exception):
    def __init__(self, path):
        self.path = path


class ConfiguratorIsCalledBeforeInitConfigPath(Exception):
    pass


class BadReferenceAmpliconChar(Exception):
    pass

class BadSgRNAChar(Exception):
    pass

class AlignerSubstitutionDoesntExist(Exception):
    def __init__(self, name):
        self.name = name


class ClassificationFailed(Exception):
    pass

class PriorPositionHasWrongLength(Exception):
    def __init__(self, expected_len, actual_len, indel_type, indel_min, indel_max):
        self.expected_len = expected_len
        self.actual_len = actual_len
        self.indel_type = indel_type
        self.indel_min = indel_min
        self.indel_max = indel_max


class UnknownAlignmentChar(Exception):
    pass
