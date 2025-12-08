import pandas as pd
from tqdm import tqdm
import random
import os
import logging
logger = logging.getLogger("prepare_ptm_search")

'''
    make_ptms_df <-- make_df_of_PTMs
'''

from typing import TextIO

'''
   The script groups proteins by modification
'''

# ---------------------------------------------------/ --- /------------------------------------------------------------
def make_df_of_PTMs(group_by_mod_res: dict[str, dict[str, list[str]]], dict_acc_to_names: dict[str, str], grouped_modres_file: TextIO, config) -> pd.DataFrame:
    logger.info(f'\nView of PTM lists:')

    logger.info( '\n'.join(map(str, [f'{dict_elm} : {[*group_by_mod_res[dict_elm]][0:5]}' for dict_elm in random.sample(list(group_by_mod_res.keys()), 5)] )) )
    logger.info(f'Число ПТМ: {len(group_by_mod_res.keys())}\n')

    # creating a text file and dataframe for protein groups by their modifications
    df4 = pd.DataFrame(columns=['accession', 'Protein name', 'PTM', 'coordinate'])
    for modification in tqdm(group_by_mod_res.keys()):
        grouped_modres_file.write('\n')
        grouped_modres_file.write('// ' + modification + ' :' + '\n')
        grouped_modres_file.write('\n')

        for protein in group_by_mod_res[modification].keys():

            df4.loc[len(df4)] = [protein, dict_acc_to_names[protein], modification, '|'.join(map(str, group_by_mod_res[modification][protein] ))]
            grouped_modres_file.write(f'{protein}|{dict_acc_to_names[protein]}' + '\n')

    grouped_modres_file.close()
    df4.to_csv(config.ptm_search_dir / f'{config.experiment_name}_PTM_info_from_UniProt_{config.analysis_index}.csv', encoding='utf-8')
    return df4

# ----------------/ Opening the necessary files for reading and writing. Launching an internal function /---------------
def make_ptms_df(group_by_mod_res_: dict[str, dict[str, list[str]]], dict_acc_to_names_: dict[str, str], config) -> pd.DataFrame:
    text2 = ' Creating a dataframe of protein lists for each post-translational modification '
    number2 = int(round((200 - len(text2)) / 2, 0))
    logger.info(f'{text2:.^{number2}}')
    grouped_modres_file = open(config.ptm_search_dir / f'{config.experiment_name}_Group_by_MOD_RES_{config.analysis_index}.txt', 'w')
    return make_df_of_PTMs(group_by_mod_res_, dict_acc_to_names_, grouped_modres_file, config)
