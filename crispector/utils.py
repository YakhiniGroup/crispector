import logging
import os


class Logger:
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
            if os.path.exists(cls._OUTPUT_DIR + 'crispector_main.log'):
                os.remove(cls._OUTPUT_DIR + 'crispector_main.log')
            f_handler = logging.FileHandler(os.path.join(cls._OUTPUT_DIR, 'crispector_main.log'))

            f_handler.setLevel(cls._logger_level)
            c_handler.setLevel(cls._logger_level)

            # Create formatters and add it to handlers
            # TODO - change formating by remove filename and funcName
            f_format = logging.Formatter('%(asctime)s %(filename)s, %(funcName)s()\t %(levelname)s\t %(message)s')
            c_format = logging.Formatter('%(asctime)s %(filename)s, %(funcName)s()\t %(levelname)s\t %(message)s')
            f_handler.setFormatter(f_format)
            c_handler.setFormatter(c_format)

            # Add handlers to the logger
            logger.addHandler(f_handler)
            logger.addHandler(c_handler)
        return logger

    @classmethod
    def set_log_path(cls, path):
        cls._OUTPUT_DIR = path

    @classmethod
    def set_logger_level(cls, mode):
        cls._logger_level = mode


if __name__ == "__main__":
    print(os.path)


