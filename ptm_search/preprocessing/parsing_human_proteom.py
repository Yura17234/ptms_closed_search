from typing import TextIO
from pandas import DataFrame
from tqdm import tqdm
import random
random.seed(42)
from ptm_search.preprocessing.modification_name_correction import smaller_groups
import logging
logger = logging.getLogger("prepare_ptm_search")

'''
                    /--- smaller_groups
    parsing_human_proteom <-- get_PTMs_lists <-- get_all_accs_and_names
'''

from ptm_search.find_prot_name_sequence import get_protein_name

'''
    The code searches for modifications in the protein sequence for each protein found
'''

# ---------------------------/ The function creates a dictionary accession - protein name /-----------------------------
def get_all_accs_and_names(df: DataFrame) -> dict[str, str]:
    accs_and_names_dict = {}
    for list_accs in tqdm(df["dbname"]):
        for acc in list(list_accs[2:-2].split("', '")):
            if acc.split('|')[1] in accs_and_names_dict.keys():
                continue
            accs_and_names_dict[acc.split('|')[1]] = get_protein_name(acc.split('|')[1])

    logger.info(f'Вид словаря accession - название белка:')
    logger.info( '\n'.join(map(str, [ f'{dict_elm} : {accs_and_names_dict[dict_elm]}' for dict_elm in random.sample(list(accs_and_names_dict.keys()), 5)] )) )
    return accs_and_names_dict

# -----------------------------/ The function generates lists of proteins for each PTM /--------------------------------
def get_PTMs_lists(query_text: str, modres_file: TextIO, list_of_msms_proteins: dict[str, str]) -> dict[str, dict[str, list[str]]]:
    PTMs_groups = {}
    query_text_list = query_text.split('\n//\n')
    count = 0
    list_of_missed_modifs = []
    for description in query_text_list:
        if 'MOD_RES' not in description:
            continue
        count += 1
        description_list = description.split('\n')
        acc_of_protein = ''
        for i in range(0, len(description_list)):
            if 'ID   ' in description_list[i]:
                continue
            if 'AC   ' in description_list[i]:
                for elm in description_list[i].split()[1:]:
                    if elm[:-1] in list_of_msms_proteins.keys():
                        acc_of_protein = elm[:-1]
                        modres_file.write('// ' + f'{acc_of_protein}|{list_of_msms_proteins[acc_of_protein]}' + '\n')
                        modres_file.write('\n')
                        break
            if acc_of_protein == '':
                break  # if no accession was found in the experiment, proceed to the next query
            if 'FT   MOD_RES' in description_list[i]:
                name_of_modification = ''
                # cycle is needed to collect the entire name of the protein modification from the string
                for u in range(28, len(description_list[i + 1])):
                    if description_list[i + 1][u] != '\"':
                        name_of_modification += description_list[i + 1][u]
                    else:
                        break
                correct_prot_modif_name = smaller_groups(name_of_modification)
                if correct_prot_modif_name == None:
                    list_of_missed_modifs.append(name_of_modification)
                    continue
                if ':' in description_list[i].split()[2]:
                    position = description_list[i].split()[2].split(':')[1]
                else:
                    position = description_list[i].split()[2]

                # adding information to the dictionary
                if correct_prot_modif_name not in PTMs_groups.keys():
                    PTMs_groups[correct_prot_modif_name] = {}

                if acc_of_protein not in PTMs_groups[correct_prot_modif_name].keys():
                    # dictionary with modification coordinates
                    PTMs_groups[correct_prot_modif_name][acc_of_protein] = [position]
                    modres_file.write(correct_prot_modif_name + '\n')
                    modres_file.write('\n')
                    continue

                # if the position is not in the list
                if position not in PTMs_groups[correct_prot_modif_name][acc_of_protein]:
                    # dictionary with modification coordinates (adding coordinates)
                    PTMs_groups[correct_prot_modif_name][acc_of_protein] += [position]
                modres_file.write(correct_prot_modif_name + '\n')
                modres_file.write('\n')

    logger.info(f'Number of protein annotations: {count}')
    modres_file.close()

    logger.info(f'Dictionary of modifications and their proteins:')
    logger.info('\n'.join(map(str, [f'{dict_elm} : {[*PTMs_groups[dict_elm].keys()][0:5]}' for dict_elm in random.sample(list(PTMs_groups.keys()), 5)])))
    logger.info(f'The number of modifications that the built-in dictionary did not recognize: {len(set(list_of_missed_modifs))}\n')
    logger.info(f'The list of modifications that the built-in dictionary did not recognize: {set(list_of_missed_modifs)}\n')
    return PTMs_groups

# -----------------/ Opening the necessary files for reading and writing. Launching internal functions /----------------
def parsing_human_proteom(config, dataframe: DataFrame) -> tuple[dict[str, dict[str, list[str]]], dict[str, str]]:
    text1 = ' Proteome parsing in UniProt to search for PTMs '
    number1 = int(round((200 - len(text1)) / 2, 0))
    logger.info(f'{text1:.^{number1}}')
    # proteome from UniProt
    with open(config.uniprot_query_path, 'r') as query_file:
        query_text = query_file.read()
    logger.info(config.uniprot_query_path)

    # writing protein modifications in a text file
    modres_file = open(config.ptm_search_dir / f'{config.experiment_name}_MOD_RES_{config.analysis_index}.txt', 'w')
    logger.info(config.ptm_search_dir / f'{config.experiment_name}_MOD_RES_{config.analysis_index}.txt')

    prot_acc_and_names = get_all_accs_and_names(dataframe)
    return get_PTMs_lists(query_text, modres_file, prot_acc_and_names), prot_acc_and_names
