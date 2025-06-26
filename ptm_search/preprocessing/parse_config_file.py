'''
    Парсер config-фалов конфигурации запуска Identipy
'''
import os

import configparser
from configupdater import ConfigUpdater

'''
    main3 <-- parse_config_file
'''

import json
with open(f"{os.getcwd()}/PTM_search_6_0/ptm_name_to_config_ptm_name_dict.json", "r") as ptm_name_to_config_ptm_name_file:
    ptm_name_to_config_ptm_name_dict_json = ptm_name_to_config_ptm_name_file.read()
dict_of_modifications = json.loads(ptm_name_to_config_ptm_name_dict_json)

# ===================================================/ === /============================================================
def make_config_files(list_of_MOD_RES, path, x, x2, path2, full_path_to_first_search_Config, variant_fo_search_):
    print(f'Список всех модификаций:\n{list_of_MOD_RES}\n')

    # Config-файл с параметрами стандартного запуска identipy
    config = configparser.ConfigParser()
    config.sections()
    config.read(full_path_to_first_search_Config)

    try:
        if variant_fo_search_ == 1:
            os.mkdir(f'{path}/{x}_Configs_full_search')
        elif variant_fo_search_ == 2:
            os.mkdir(f'{path}/{x}_Configs_fast_search')
    except OSError as error:
        print(error)

    f = []
    count = 0
    for modif in list_of_MOD_RES:
        f.append(modif)
        try:
            mification_name_for_Identipy = dict_of_modifications[modif]
        except:
            count += 1
            continue

        # Config-файл для PTM поиска
        updater = ConfigUpdater()
        updater.read(f"{os.getcwd()}/PTM_search_6_0/PTM_default.cfg")

        # Добавляем название модификации, которое поймет Identipy
        if config["modifications"]["variable"] == '':
            updater["modifications"]["variable"].value = mification_name_for_Identipy
        elif config["modifications"]["variable"] != '':
            updater["modifications"]["variable"].value = f'{config["modifications"]["variable"]}, {mification_name_for_Identipy}'

        if variant_fo_search_ == 1: # Добавление полного пути к полной базе данных поиска
            updater["input"]["database"].value = f'{path2}/{x2}_PTMs_search_{x}_reverse.fasta'
        elif variant_fo_search_ == 2: # Добавление полного пути к малым базам данных поиска
            modif_for_fasta = modif.replace(' ', '_').replace(';', '').replace('/', '_')
            updater["input"]["database"].value = f'{path2}/{x}_fasta_for_fast_search/{modif_for_fasta}_{x}_reverse.fasta'

        # Перенос остальных параметров стандартного поиска
        updater["search"]["number of missed cleavages"].value = config["search"]["number of missed cleavages"]

        updater["performance"]["processes"].value = config["performance"]["processes"]

        # Сохранение config-файла
        modif2 = modif.replace(' ', '_').replace(';', '').replace('/', '_')
        if variant_fo_search_ == 1:
            file = open(f'{path}/{x}_Configs_full_search/{modif2}_{x2}.cfg', 'w')
            file.write(str(updater))
            file.close()
        if variant_fo_search_ == 2:
            file = open(f'{path}/{x}_Configs_fast_search/{modif2}_{x2}.cfg', 'w')
            file.write(str(updater))
            file.close()

    f2 = set(f)
    f2 = [i for i in f2 if i not in dict_of_modifications.keys()]
    print(f'Список модификаций не из словаря:\n{f2}\n')
    print(len(f2), len(f))
    print(f'Число неправильно прописанных для Identipy модификаций: {count}')
    print(f'Число модификаций {len(list_of_MOD_RES)}')

# =====================/ Открытие необходимых файлов на чтение и запись. Запуск внутренных функций /====================
def main4(list_of_MOD_RES_, path_, x_, x2_, path2_, full_path_to_first_search_Config_, variant_fo_search):
    text4 = ' Создание конфигурационных файлов для поиска PTM '
    number4 = int(round((200 - len(text4)) / 2, 0))
    print(f'{text4:.^{number4}}')
    make_config_files(list_of_MOD_RES_, path_, x_, x2_, path2_, full_path_to_first_search_Config_, variant_fo_search)