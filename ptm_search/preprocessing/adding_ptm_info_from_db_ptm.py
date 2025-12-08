import json
import gzip
from tqdm import tqdm
from pathlib import Path
import logging
logger = logging.getLogger("prepare_ptm_search")

def merge_dicts(dict_of_grouped_by_ptms_proteins, dict_of_proteins_by_ptms_additional, proteins_from_msms) -> dict[str, dict[str, list[str]]]:

    for ptm_name in tqdm(list(dict_of_grouped_by_ptms_proteins.keys()) + ['Ubiquitinlysine', 'Sumoyllysine']):
        if ptm_name not in dict_of_proteins_by_ptms_additional.keys():
            continue
        for protein_with_ptm in dict_of_proteins_by_ptms_additional[ptm_name].keys():
            if protein_with_ptm not in proteins_from_msms:
                continue
            if ptm_name not in dict_of_grouped_by_ptms_proteins.keys():
                dict_of_grouped_by_ptms_proteins[ptm_name] = {}
            if protein_with_ptm not in dict_of_grouped_by_ptms_proteins[ptm_name].keys():
                dict_of_grouped_by_ptms_proteins[ptm_name][protein_with_ptm] = dict_of_proteins_by_ptms_additional[ptm_name][protein_with_ptm]
                continue
            if protein_with_ptm in dict_of_grouped_by_ptms_proteins[ptm_name].keys():
                dict_of_grouped_by_ptms_proteins[ptm_name][protein_with_ptm] = list(set(dict_of_grouped_by_ptms_proteins[ptm_name][protein_with_ptm] + dict_of_proteins_by_ptms_additional[ptm_name][protein_with_ptm]))
                continue

    return dict_of_grouped_by_ptms_proteins

def adding_ptm_info_from_db_ptm(dict_of_grouped_by_ptms_proteins_dbptm: dict[str, dict[str, list[str]]], proteins_from_msms_dbptm: list[str]) -> dict[str, dict[str, list[str]]]:
    module_dir = Path(__file__).parent.parent.resolve()
    with gzip.open(module_dir / 'data' / 'dict_of_proteins_by_ptms_from_dbPTM.json.gz', 'rt', encoding='utf-8') as dbptm_info_file:
        dict_of_proteins_by_ptms_from_dbptm_json = dbptm_info_file.read()
    dict_of_proteins_by_ptms_from_dbptm = json.loads(dict_of_proteins_by_ptms_from_dbptm_json)

    text1_2 = ' Adding information from the dbPTM database '
    number1_2 = int(round((200 - len(text1_2)) / 2, 0))
    logger.info(f'\n{text1_2:.^{number1_2}}\n')
    logger.info('\n'.join(map(str, [f'{ptm_key} -- {len(dict_of_grouped_by_ptms_proteins_dbptm[ptm_key].keys())}' for ptm_key in dict_of_grouped_by_ptms_proteins_dbptm.keys()])))

    dict_of_grouped_by_ptms_proteins_dbptm = merge_dicts(dict_of_grouped_by_ptms_proteins_dbptm, dict_of_proteins_by_ptms_from_dbptm, proteins_from_msms_dbptm)

    logger.info('\n'.join(map(str, [f'{ptm_key} -- {len(dict_of_grouped_by_ptms_proteins_dbptm[ptm_key].keys())}' for ptm_key in dict_of_grouped_by_ptms_proteins_dbptm.keys()])))

    return dict_of_grouped_by_ptms_proteins_dbptm

def adding_ptm_info_from_additional_lists(dict_of_grouped_by_ptms_proteins_add_lists: dict[str, dict[str, list[str]]], proteins_from_msms_add_lists: list[str], config) -> dict[str, dict[str, list[str]]]:
    with open(config.additional_lists_path, 'r') as additional_lists_file:
        additional_lists_text = additional_lists_file.read()
    logger.info(config.additional_lists_path)

    dict_of_proteins_by_ptms_from_add_lists = {}
    for ptm_add in additional_lists_text.split('>'):
        ptm_list = ptm_add.split('\n')
        dict_of_proteins_by_ptms_from_add_lists[ptm_list[0]] = dict.fromkeys(ptm_list[1:], [])

    text1_3 = ' Adding information from additional lists '
    number1_3 = int(round((200 - len(text1_3)) / 2, 0))
    logger.info(f'\n{text1_3:.^{number1_3}}\n')

    dict_of_grouped_by_ptms_proteins_add_lists = merge_dicts(dict_of_grouped_by_ptms_proteins_add_lists, dict_of_proteins_by_ptms_from_add_lists, proteins_from_msms_add_lists)

    return dict_of_grouped_by_ptms_proteins_add_lists
