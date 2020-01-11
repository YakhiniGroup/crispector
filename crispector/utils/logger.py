import logging
import os

from utils.constants_and_types import Path


class Logger:
    """
    Singleton logger based on logging package.
    Dump all messages both to shell and crispector_main.log.
    """
    _configured = False
    _OUTPUT_DIR = None
    _logger_level = logging.DEBUG
    _logger_path = None
    logger_name = 'crispector_main.log'

    @classmethod
    def get_logger(cls):
        # Create a custom logger
        logger = logging.getLogger("default")
        logger.level = cls._logger_level

        if not cls._configured:
            cls._configured = True

            # Create handlers
            c_handler = logging.StreamHandler()
            cls._logger_path = os.path.join(cls._OUTPUT_DIR, cls.logger_name)
            if os.path.exists(cls._logger_path):
                os.remove(cls._logger_path)
            f_handler = logging.FileHandler(cls._logger_path)

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
    def set_output_dir(cls, path: Path):
        cls._OUTPUT_DIR = path

    @classmethod
    def get_log_path(cls) -> Path:
        return cls._logger_path

    @classmethod
    def set_logger_level(cls, mode):
        cls._logger_level = mode