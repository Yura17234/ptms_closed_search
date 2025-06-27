'''
   Скрипт выполняет группировку белков по модификации
'''
import pandas as pd
from tqdm import tqdm
import random
import os

'''
    make_ptms_df <-- make_df_of_PTMs
'''

# ---------------------------------------------------/ --- /------------------------------------------------------------
def make_df_of_PTMs(group_by_mod_res, dict_acc_to_names, grouped_modres_file, config):
    print(f'\nВид списков ПТМ:')

    print( '\n'.join(map(str, [f'{dict_elm} : {[*group_by_mod_res[dict_elm]][0:5]}' for dict_elm in random.sample(group_by_mod_res.keys(), 5)] )) )
    print(f'Число ПТМ: {len(group_by_mod_res.keys())}\n')

    # Создание текстового файла и датафрейма по группам белков по их модификациям
    dict = {'accession': [],
            'Protein name': [],
            'PTM': [],
            'coordinate': []
            }
    df4 = pd.DataFrame(dict)
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

# ---------------------/ Открытие необходимых файлов на чтение и запись. Запуск внутренных функций /--------------------
def make_ptms_df(group_by_mod_res_, dict_acc_to_names_, config):
    text2 = ' Создание датафрейма списков белков по каждой пост-трансляционной модификации '
    number2 = int(round((200 - len(text2)) / 2, 0))
    print(f'{text2:.^{number2}}')
    grouped_modres_file = open(config.ptm_search_dir / f'{config.experiment_name}_Group_by_MOD_RES_{config.analysis_index}.txt', 'w')
    return make_df_of_PTMs(group_by_mod_res_, dict_acc_to_names_, grouped_modres_file, config)
