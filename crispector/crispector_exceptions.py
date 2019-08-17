
class FastpRunTimeError(Exception):
    pass


class NoneValuesInAmpliconsCSV(Exception):
    pass


class SgRNANotInReferenceSequence(Exception):
    def __init__(self, site_name):
        self.site_name = site_name


class CantOpenMergedFastqFile(Exception):
    def __init__(self, path):
        self.path = path

class ConfiguratorIsCalledBeforeInitConfigPath(Exception):
    pass
