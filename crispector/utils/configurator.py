import os

import yaml

from utils.constants_and_types import Path
from utils.exceptions import ConfiguratorIsCalledBeforeInitConfigPath


class Configurator:
    """
    Singleton YAML configurator based on yaml package.
    """

    _config_file = None
    _CONFIG_PATH = None

    @classmethod
    def set_cfg_path(cls, path: Path):
        if path is None:
            path = os.path.join(os.path.dirname(os.path.abspath(__file__)), '../config/default_config.yml')
        cls._CONFIG_PATH = path

    @classmethod
    def get_cfg(cls):
        if cls._config_file is None:
            if cls._CONFIG_PATH is None:
                raise ConfiguratorIsCalledBeforeInitConfigPath()

            # Read YAML file
            with open(cls._CONFIG_PATH, 'r') as stream:
                cls._config_file = yaml.safe_load(stream)

        return cls._config_file