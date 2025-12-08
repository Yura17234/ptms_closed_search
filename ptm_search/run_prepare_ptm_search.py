import argparse
from ptm_search.config import Config
from ptm_search.setup_logging import setup_logger
from ptm_search.preprocessing.prepare_ptm_search import prepare_ptm_search

def main():
    parser = argparse.ArgumentParser(description="Prepare PTM database and configs")
    parser.add_argument('--config', type=str, default='parameters.cfg', help='Path to config file')
    args = parser.parse_args()

    config = Config(args.config)
    logger = setup_logger(config.ptm_search_dir / 'logs', name="prepare_ptm_search")

    logger.info('Run prepare_ptm_search !')
    prepare_ptm_search(config)
    logger.info('prepare_ptm_search -- Done !')
