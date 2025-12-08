import argparse
from ptm_search.config import Config
from ptm_search.setup_logging import setup_logger
from ptm_search.search.multiple_search import multiple_search

def main():
    parser = argparse.ArgumentParser(description="Run IdentiPy PTM search")
    parser.add_argument('--config', type=str, default='parameters.cfg', help='Path to config file')
    args = parser.parse_args()

    config = Config(args.config)
    logger = setup_logger(config.ptm_search_dir / 'logs', name="multiple_search")

    logger.info('Run multiple_search !')
    multiple_search(config)
    logger.info('multiple_search -- Done !')
