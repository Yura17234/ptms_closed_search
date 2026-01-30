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
    df['accession'] = df["dbname"].str.split('|').str[1]
    df['protein_name'] = df['accession'].apply(lambda x: get_protein_name(x))
    accs_and_names_dict = df.set_index('accession')['protein_name'].to_dict()

    logger.info(f'Dict: accession - protein name')
    sampled_keys = random.sample(accs_and_names_dict.keys(), 5)
    logger.info('\n'.join((f'{dict_elm} : {accs_and_names_dict[dict_elm]}' for dict_elm in sampled_keys)))
    return accs_and_names_dict

def uniprot_query_generator(query_f: TextIO):
    buf = []
    for line in query_f:
        if line.strip() == '//':
            if buf:
                yield ''.join(buf)
                buf = []
        else:
            buf.append(line)
    if buf:
        yield ''.join(buf)

# -----------------------------/ The function generates lists of proteins for each PTM /--------------------------------
def get_PTMs_lists(query_file: TextIO, modres_file: TextIO, list_of_msms_proteins: dict[str, str]) -> dict[str, dict[str, set[str]]]:
    PTMs_groups, count, list_of_missed_modifs = {}, 0, []
    for description in uniprot_query_generator(query_file):
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
                    if elm[:-1] in list_of_msms_proteins:
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
                if correct_prot_modif_name not in PTMs_groups:
                    PTMs_groups[correct_prot_modif_name] = {}

                if acc_of_protein not in PTMs_groups[correct_prot_modif_name]:
                    # dictionary with modification coordinates
                    PTMs_groups[correct_prot_modif_name][acc_of_protein] = {position}
                    modres_file.write(correct_prot_modif_name + '\n')
                    modres_file.write('\n')
                    continue

                # if the position is not in the list
                # dictionary with modification coordinates (adding coordinates)
                PTMs_groups[correct_prot_modif_name][acc_of_protein].add(position)
                modres_file.write(correct_prot_modif_name + '\n')
                modres_file.write('\n')

    logger.info(f'Number of protein annotations: {count}')
    modres_file.close()

    logger.info(f'Dictionary of modifications and their proteins:')
    random_ptms_groups = random.sample(list(PTMs_groups.keys()), 5)
    logger.info('\n'.join((f'{dict_elm} : {[*PTMs_groups[dict_elm].keys()][0:5]}' for dict_elm in random_ptms_groups)))
    logger.info(f'The number of modifications that the built-in dictionary did not recognize: {len(set(list_of_missed_modifs))}\n')
    logger.info(f'The list of modifications that the built-in dictionary did not recognize: {set(list_of_missed_modifs)}\n')
    return PTMs_groups

# -----------------/ Opening the necessary files for reading and writing. Launching internal functions /----------------
def parsing_proteom_query(config, dataframe: DataFrame) -> tuple[dict[str, dict[str, set[str]]], dict[str, str]]:
    text1 = ' Proteome parsing in UniProt to search for PTMs '
    number1 = int(round((200 - len(text1)) / 2, 0))
    logger.info(f'{text1:.^{number1}}')

    # writing protein modifications in a text file
    logger.info(config.ptm_search_dir / f'{config.experiment_name}_MOD_RES_{config.analysis_index}.txt')
    modres_file = open(config.ptm_search_dir / f'{config.experiment_name}_MOD_RES_{config.analysis_index}.txt', 'w')

    prot_acc_and_names = get_all_accs_and_names(dataframe)

    logger.info(config.uniprot_query_path)
    # proteome from UniProt
    with open(config.uniprot_query_path, 'r') as query_file:
        ptms_groups = get_PTMs_lists(query_file, modres_file, prot_acc_and_names)
    return ptms_groups, prot_acc_and_names
