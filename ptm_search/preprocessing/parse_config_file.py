'''
    Парсер config-фалов конфигурации запуска Identipy
'''
import os
from pathlib import Path
import configparser
from configupdater import ConfigUpdater
import json

'''
    parse_config_file <-- make_config_files
'''

# ---------------------------------------------------/ --- /------------------------------------------------------------
def make_config_files(list_of_MOD_RES, config, variant_of_search_):
    print(f'Список всех модификаций:\n{list_of_MOD_RES}\n')

    module_dir = Path(__file__).parent.parent.resolve()
    with open(module_dir / 'data' / 'ptm_name_to_config_ptm_name_dict.json', 'r') as ptm_name_to_config_ptm_name_file:
        ptm_name_to_config_ptm_name_dict_json = ptm_name_to_config_ptm_name_file.read()
    dict_of_modifications = json.loads(ptm_name_to_config_ptm_name_dict_json)

    # Config-файл с параметрами стандартного запуска identipy
    st_search_config = configparser.ConfigParser()
    st_search_config.sections()
    st_search_config.read(str(config.base_config_path))
    print(config.base_config_path)
    print(Path(config.base_config_path).exists())

    configs_full_search_dir = '_'
    configs_fast_search_dir = '_'
    if variant_of_search_ == 1:
        configs_full_search_dir = config.ptm_search_dir / f'{config.experiment_name}_Configs_full_search'
        os.makedirs(configs_full_search_dir, exist_ok=True)
    elif variant_of_search_ == 2:
        configs_fast_search_dir = config.ptm_search_dir / f'{config.experiment_name}_Configs_fast_search'
        os.makedirs(configs_fast_search_dir, exist_ok=True)

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
        updater.read(module_dir / 'configs' / f'PTM_default.cfg')

        # Добавляем название модификации, которое поймет Identipy
        if st_search_config["modifications"]["variable"] == '':
            updater["modifications"]["variable"].value = mification_name_for_Identipy
        elif st_search_config["modifications"]["variable"] != '':
            updater["modifications"]["variable"].value = f'{st_search_config["modifications"]["variable"]}, {mification_name_for_Identipy}'

        if variant_of_search_ == 1: # Добавление полного пути к полной базе данных поиска
            updater["input"]["database"].value = f'{str(config.ptm_search_dir)}/{config.experiment_name}_PTMs_search_{config.analysis_index}_reverse.fasta'
        elif variant_of_search_ == 2: # Добавление полного пути к малым базам данных поиска
            modif_for_fasta = modif.replace(' ', '_').replace(';', '').replace('/', '_')
            updater["input"]["database"].value = f'{str(config.ptm_search_dir)}/{config.experiment_name}_fasta_for_fast_search/{modif_for_fasta}_{config.analysis_index}_reverse.fasta'

        # Перенос остальных параметров стандартного поиска
        updater["search"]["number of missed cleavages"].value = st_search_config["search"]["number of missed cleavages"]

        updater["performance"]["processes"].value = st_search_config["performance"]["processes"]

        # Сохранение config-файла
        modif2 = modif.replace(' ', '_').replace(';', '').replace('/', '_')
        if variant_of_search_ == 1:
            with open(configs_full_search_dir / f'{modif2}_{config.analysis_index}.cfg', 'w') as file:
                file.write(str(updater))
        if variant_of_search_ == 2:
            with open(configs_fast_search_dir / f'{modif2}_{config.analysis_index}.cfg', 'w') as file:
                file.write(str(updater))

    f2 = set(f)
    f2 = [i for i in f2 if i not in dict_of_modifications.keys()]
    print(f'Список модификаций не из словаря:\n{f2}\n')
    print(len(f2), len(f))
    print(f'Число неправильно прописанных для Identipy модификаций: {count}')
    print(f'Число модификаций {len(list_of_MOD_RES)}')

# ---------------------/ Открытие необходимых файлов на чтение и запись. Запуск внутренных функций /--------------------
def parse_config_file(list_of_MOD_RES_, config, variant_of_search):
    text4 = ' Создание конфигурационных файлов для поиска PTM '
    number4 = int(round((200 - len(text4)) / 2, 0))
    print(f'{text4:.^{number4}}')
    make_config_files(list_of_MOD_RES_, config, variant_of_search)
