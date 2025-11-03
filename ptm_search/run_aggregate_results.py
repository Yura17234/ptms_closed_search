import argparse
from ptm_search.config import Config
from ptm_search.setup_logging import setup_logger
from ptm_search.postprocessing.aggregate_results import aggregate_results

def main():
    parser = argparse.ArgumentParser(description="Aggregate PTM search results")
    parser.add_argument('--config', type=str, default='pipeline.cfg', help='Path to config file')
    args = parser.parse_args()

    config = Config(args.config)
    logger = setup_logger(config.ptm_search_dir / 'logs', name="aggregate_results")

    logger.info('Run aggregate_results !')
    aggregate_results(config)
    logger.info('aggregate_results -- Done !')
