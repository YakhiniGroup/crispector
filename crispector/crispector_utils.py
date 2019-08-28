import logging
import os
import yaml # TODO - add to project requirements
from crispector_exceptions import ConfiguratorIsCalledBeforeInitConfigPath
from crispector_types import Path


class Logger:
    """
    Singleton logger based on logging package.
    Dump all messages both to shell and crispector_main.log.
    """
    _configured = False
    _OUTPUT_DIR = None
    _logger_level = logging.DEBUG

    @classmethod
    def get_logger(cls):
        # Create a custom logger
        logger = logging.getLogger("default")
        logger.level = cls._logger_level

        if not cls._configured:
            cls._configured = True

            # Create handlers
            c_handler = logging.StreamHandler()
            logger_path = os.path.join(cls._OUTPUT_DIR, 'crispector_main.log')
            if os.path.exists(logger_path):
                os.remove(logger_path)
            f_handler = logging.FileHandler(logger_path)

            f_handler.setLevel(cls._logger_level)
            c_handler.setLevel(cls._logger_level)

            # Create formatters and add it to handlers
            f_format = logging.Formatter('%(asctime)s %(levelname)s\t %(message)s')
            c_format = logging.Formatter('%(asctime)s %(levelname)s\t %(message)s')
            f_handler.setFormatter(f_format)
            c_handler.setFormatter(c_format)

            # Add handlers to the logger
            logger.addHandler(f_handler)
            logger.addHandler(c_handler)
        return logger

    @classmethod
    def set_log_path(cls, path: Path):
        cls._OUTPUT_DIR = path

    @classmethod
    def set_logger_level(cls, mode):
        cls._logger_level = mode


class Configurator:
    """
    Singleton YAML configurator based on yaml package.
    """

    _config_file = None
    _CONFIG_PATH = None

    @classmethod
    def set_cfg_path(cls, path: Path):
        if path is None:
            path = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'config/default_config.yml')
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
