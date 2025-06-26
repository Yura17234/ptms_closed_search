import json
import os
from tqdm import tqdm

def adding_ptm_info_from_db_ptm(dict_of_grouped_by_ptms_proteins, proteins_from_msms):
    with open(f"{os.getcwd()}/PTM_search_6_0/dict_of_proteins_by_ptms_from_dbPTM.json", "r") as dbptm_info_file:
        dict_of_proteins_by_ptms_from_dbptm_json = dbptm_info_file.read()# Достаем информацию о PTM белков из dbPTM
    dict_of_proteins_by_ptms_from_dbptm = json.loads(dict_of_proteins_by_ptms_from_dbptm_json)

    print('\n.......... Добавление информации из базы данных dbPTM ..........\n')
    print('\n'.join(map(str, [f'{ptm_key} -- {len(dict_of_grouped_by_ptms_proteins[ptm_key].keys())}' for ptm_key in dict_of_grouped_by_ptms_proteins.keys()])))
    for ptm_name in tqdm(list(dict_of_grouped_by_ptms_proteins.keys()) + ['Ubiquitinlysine', 'Sumoyllysine']):
        if ptm_name not in dict_of_proteins_by_ptms_from_dbptm.keys():
            continue
        for protein_with_ptm in dict_of_proteins_by_ptms_from_dbptm[ptm_name].keys():
            if protein_with_ptm not in proteins_from_msms:
                continue
            if ptm_name not in dict_of_grouped_by_ptms_proteins.keys():
                dict_of_grouped_by_ptms_proteins[ptm_name] = {}
            if protein_with_ptm not in dict_of_grouped_by_ptms_proteins[ptm_name].keys():
                dict_of_grouped_by_ptms_proteins[ptm_name][protein_with_ptm] = dict_of_proteins_by_ptms_from_dbptm[ptm_name][protein_with_ptm]
                continue
            if protein_with_ptm in dict_of_grouped_by_ptms_proteins[ptm_name].keys():
                dict_of_grouped_by_ptms_proteins[ptm_name][protein_with_ptm] = list(set(dict_of_grouped_by_ptms_proteins[ptm_name][protein_with_ptm] + dict_of_proteins_by_ptms_from_dbptm[ptm_name][protein_with_ptm]))
                continue

    print('\n'.join(map(str, [f'{ptm_key} -- {len(dict_of_grouped_by_ptms_proteins[ptm_key].keys())}' for ptm_key in dict_of_grouped_by_ptms_proteins.keys()])))

    return dict_of_grouped_by_ptms_proteins
