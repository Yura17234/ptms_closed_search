import argparse
from ptm_search.config import Config
from ptm_search.postprocessing.aggregate_results import aggregate_results

def main():
    parser = argparse.ArgumentParser(description="Aggregate PTM search results")
    parser.add_argument('--config', type=str, default='pipeline.cfg', help='Path to config file')
    args = parser.parse_args()

    config = Config(args.config)
    print('Run aggregate result !')
    aggregate_results(config)
