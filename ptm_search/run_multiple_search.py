import argparse
from ptm_search.config import Config
from ptm_search.search.multiple_search import multiple_search

def main():
    parser = argparse.ArgumentParser(description="Run IdentiPy PTM search")
    parser.add_argument('--config', type=str, default='pipeline.cfg', help='Path to config file')
    args = parser.parse_args()

    config = Config(args.config)
    print('Run multiple search !')
    multiple_search(config)
