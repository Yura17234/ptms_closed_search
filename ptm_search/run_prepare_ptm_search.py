import argparse
from ptm_search.config import Config
from ptm_search.preprocessing.prepare_ptm_search import prepare_ptm_search

def main():
    parser = argparse.ArgumentParser(description="Prepare PTM database and configs")
    parser.add_argument('--config', type=str, default='pipeline.cfg', help='Path to config file')
    args = parser.parse_args()

    config = Config(args.config)
    print('Run prepare ptm search !')
    prepare_ptm_search(config)
