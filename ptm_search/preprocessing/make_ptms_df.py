'''
   Скрипт выполняет группировку белков по модификации
'''
import pandas as pd
from tqdm import tqdm
import random

'''
    main2 <-- make_df_of_PTMs
'''

# ===================================================/ === /============================================================
def make_df_of_PTMs(Group_by_MOD_RES, dict_acc_to_names, file_3, path_, x_, x2_):
    print(f'\nВид списков ПТМ:')

    print( '\n'.join(map(str, [f'{dict_elm} : {[*Group_by_MOD_RES[dict_elm]][0:5]}' for dict_elm in random.sample(Group_by_MOD_RES.keys(), 5)] )) )
    print(f'Число ПТМ: {len(Group_by_MOD_RES.keys())}\n')

    # Создание текстового файла и датафрейма по группам белков по их модификациям
    dict = {'accession': [],
            'Protein name': [],
            'PTM': [],
            'coordinate': []
            }
    df4 = pd.DataFrame(dict)
    for modification in tqdm(Group_by_MOD_RES.keys()):
        file_3.write('\n')
        file_3.write('// ' + modification + ' :' + '\n')
        file_3.write('\n')

        for protein in Group_by_MOD_RES[modification].keys():

            df4.loc[len(df4)] = [protein, dict_acc_to_names[protein], modification, '|'.join(map(str, Group_by_MOD_RES[modification][protein] ))]
            file_3.write(f'{protein}|{dict_acc_to_names[protein]}' + '\n')

    file_3.close()
    df4.to_csv(f'{path_}/{x2_}_{x_}_PTM_info_from_UniProt.csv', encoding='utf-8')
    return df4

# =====================/ Открытие необходимых файлов на чтение и запись. Запуск внутренных функций /====================
def main2(Group_by_MOD_RES_, dict_acc_to_names_, path, x, x2):
    text2 = ' Создание датафрейма списков белков по каждой пост-трансляционной модификации '
    number2 = int(round((200 - len(text2)) / 2, 0))
    print(f'{text2:.^{number2}}')
    file4 = open(f'{path}/{x2}_Group_by_MOD_RES_{x}.txt', 'w')
    return make_df_of_PTMs(Group_by_MOD_RES_, dict_acc_to_names_, file4, path, x, x2)