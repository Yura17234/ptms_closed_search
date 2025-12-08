import logging
import os
from datetime import datetime

def setup_logger(log_dir, name: str, level=logging.INFO) -> logging.Logger:
    os.makedirs(log_dir, exist_ok=True)

    logger = logging.getLogger(name)
    logger.setLevel(logging.DEBUG)
    logger.propagate = False

    if logger.hasHandlers():
        logger.handlers.clear()

    timestamp = datetime.now().strftime('%H%M%S')
    pid = os.getpid()
    log_path = log_dir / f'{name}_{timestamp}_{pid}.log'

    file_handler = logging.FileHandler(log_path, encoding='utf-8')
    file_handler.setLevel(level)
    file_formatter = logging.Formatter('[%(asctime)s] | %(levelname)s | %(name)s | %(message)s',
                                       datefmt='%H:%M:%S')
    file_handler.setFormatter(file_formatter)
    logger.addHandler(file_handler)

    console_handler = logging.StreamHandler()
    console_handler.setLevel(level)
    console_formatter = logging.Formatter('[%(asctime)s] | %(levelname)s | %(name)s | %(message)s',
                                          datefmt='%H:%M:%S')
    console_handler.setFormatter(console_formatter)
    logger.addHandler(console_handler)

    return logger
