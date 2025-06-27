from pathlib import Path
import pandas as pd
import os
import glob

from ptm_search.preprocessing.parsing_human_proteom import parsing_human_proteom
from ptm_search.preprocessing.adding_ptm_info_from_db_ptm import adding_ptm_info_from_db_ptm
from ptm_search.preprocessing.make_ptms_df import make_ptms_df
from ptm_search.preprocessing.make_fasta_file_for_searches import make_fasta_file_for_searches
from ptm_search.preprocessing.parse_config_file import parse_config_file
from ptm_search.preprocessing.make_mgfs_for_ptm import make_mgfs_for_ptm

def prepare_ptm_search(config):
    '''
        Prepare PTM search
    '''

    mgf_dir = config.st_search_dir / "full_mgf_files"

    os.makedirs(config.ptm_search_dir, exist_ok=True)
    os.makedirs(config.st_search_dir, exist_ok=True)
    os.makedirs(mgf_dir, exist_ok=True)

    for ext in ['*.mgf', '*.tsv', '*.pep.xml', '*.png']:
        for file in glob.glob(str(config.work_dir / ext)):
            dest = mgf_dir if ext == '*.mgf' else st_search_path
            shutil.move(file, dest)

    st_search_df = pd.read_csv(config.st_search_dir / f"union_proteins.tsv", sep='\t')
    print(config.st_search_dir / f"union_proteins.tsv")

    ''' 1 '''
    grouped_prots_by_ptms_dict, acc_to_names_dict = parsing_human_proteom(config, st_search_df)

    # print(f'Количество белков из стандартного начального поиска: {len(list(dict_acc_to_names.keys()))}')
    # print(st_search_df.head())
    grouped_prots_by_ptms_dict = adding_ptm_info_from_db_ptm(grouped_prots_by_ptms_dict, list(acc_to_names_dict.keys()))

    ''' 2 '''
    ptm_info_df = make_ptms_df(grouped_prots_by_ptms_dict, acc_to_names_dict, config)
    # print(ptm_info_df.head())
    # print()

    ptm_counts = ptm_info_df.drop_duplicates(subset=['accession', 'PTM'], keep='first')['PTM'].value_counts(
        sort=True).head(15)
    ptm_counts.plot(kind='barh', title="Top PTMs (UniProt)", figsize=(5, 8), color=(0.1, 0.7, 0.1, 0.7),
                    edgecolor='black')
    plt.xlabel("Protein count")
    plt.savefig(config.ptm_search_dir / "PTMs_count_of_potential_prots_ENG.png", dpi=200,
                bbox_inches='tight')
    plt.close()

    ''' 3 '''
    with open(config.fasta_path, 'r') as fasta_file:
        make_fasta_file_for_searches(list(ptm_info_df['accession'].unique()), grouped_prots_by_ptms_dict, fasta_file, config)

    ''' 4 '''
    ptm_list = list(ptm_info_df['PTM'].unique())
    for mode in [1, 2]:
        parse_config_file(ptm_list, config, mode)

    ''' 5 '''
    # Создание mgf-файлов из не идентифицированных спектров для PTM-поиска
    if all(['_for_PTM.mgf' not in file for file in os.listdir(f'{config.ptm_search_dir}')]):
        make_mgfs_for_ptm(mgf_dir, config)

    print(config.ptm_search_dir)
    print(config.st_search_dir)
    print(mgf_dir)
    print('prepare_ptm_search -- connected')
