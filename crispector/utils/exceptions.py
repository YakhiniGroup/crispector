
class FastpRunTimeError(Exception):
    pass


class BadInputError(Exception):
    def __init__(self, msg):
        self.msg = msg


class SgRNANotInReferenceSequence(Exception):
    def __init__(self, site_name):
        self.site_name = site_name


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
