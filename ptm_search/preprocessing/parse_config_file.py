from typing import NoReturn
import os
from pathlib import Path
import configparser
from configupdater import ConfigUpdater
import json
import logging
logger = logging.getLogger("prepare_ptm_search")

'''
    parse_config_file <-- make_config_files
'''

'''
    Parser of config files for Identipy PTM search
'''

# ---------------------------------------------------/ --- /------------------------------------------------------------
def make_config_files(list_of_MOD_RES: list[str], config, variant_of_search_: int) -> NoReturn:
    logger.info(f'List of all modifications:\n{list_of_MOD_RES}\n')

    module_dir = Path(__file__).parent.parent.resolve()
    with open(module_dir / 'data' / 'ptm_name_to_config_ptm_name_dict.json', 'r') as ptm_name_to_config_ptm_name_file:
        dict_of_modifications = json.load(ptm_name_to_config_ptm_name_file)

    # config file with parameters of the identipy standard search
    st_search_config = configparser.ConfigParser()
    st_search_config.sections()
    st_search_config.read(config.base_config_path)

    configs_full_search_dir = None
    configs_fast_search_dir = None
    if variant_of_search_ == 1:
        configs_full_search_dir = config.ptm_search_dir / f'{config.analysis_index}_Configs_full_search'
        os.makedirs(configs_full_search_dir, exist_ok=True)
    elif variant_of_search_ == 2:
        configs_fast_search_dir = config.ptm_search_dir / f'{config.analysis_index}_Configs_fast_search'
        os.makedirs(configs_fast_search_dir, exist_ok=True)

    f, missing_modifs, count = set(), set(), 0
    for modif in list_of_MOD_RES:
        try:
            mification_name_for_Identipy = dict_of_modifications[modif]
            f.add(modif)
        except KeyError:
            count += 1
            missing_modifs.add(modif)
            continue

        # config file for PTM search
        updater = ConfigUpdater()
        updater.read(module_dir / 'configs' / f'PTM_default.cfg')

        # adding modification name for identity
        if st_search_config["modifications"]["variable"] == '':
            updater["modifications"]["variable"].value = mification_name_for_Identipy
        elif st_search_config["modifications"]["variable"] != '':
            updater["modifications"]["variable"].value = f'{st_search_config["modifications"]["variable"]}, {mification_name_for_Identipy}'

        if variant_of_search_ == 1:  # adding full path to the full search database
            updater["input"]["database"].value = f'{str(config.ptm_search_dir)}/{config.experiment_name}_PTMs_search_{config.analysis_index}_reverse.fasta'
        elif variant_of_search_ == 2:  # adding full path to the small search database
            modif_for_fasta = modif.replace(' ', '_').replace(';', '').replace('/', '_')
            updater["input"]["database"].value = f'{str(config.ptm_search_dir)}/{config.analysis_index}_fasta_for_fast_search/{modif_for_fasta}_{config.analysis_index}_reverse.fasta'

        # transferring remaining parameters of the standard search
        updater["search"]["number of missed cleavages"].value = st_search_config["search"]["number of missed cleavages"]

        updater["performance"]["processes"].value = st_search_config["performance"]["processes"]

        updater["search"]["precursor accuracy unit"].value = st_search_config["search"]["precursor accuracy unit"]
        updater["search"]["precursor accuracy left"].value = st_search_config["search"]["precursor accuracy left"]
        updater["search"]["precursor accuracy right"].value = st_search_config["search"]["precursor accuracy right"]
        updater["search"]["product accuracy unit"].value = st_search_config["search"]["product accuracy unit"]
        updater["search"]["product accuracy"].value = st_search_config["search"]["product accuracy"]

        # saving the config file
        modif2 = modif.replace(' ', '_').replace(';', '').replace('/', '_')
        if variant_of_search_ == 1:
            with open(configs_full_search_dir / f'{modif2}_{config.analysis_index}.cfg', 'w') as file:
                file.write(str(updater))
        if variant_of_search_ == 2:
            with open(configs_fast_search_dir / f'{modif2}_{config.analysis_index}.cfg', 'w') as file:
                file.write(str(updater))

    logger.info(f'The list of modifications is not from the dictionary:\n{missing_modifs}\n')
    logger.info(f'{len(missing_modifs)}')
    logger.info(f'The number of modifications incorrectly registered for Identipy: {count}')
    logger.info(f'Number of modifications {len(f)}')

# -----------------/ Opening the necessary files for reading and writing. Launching internal functions /----------------
def parse_config_file(list_of_MOD_RES_: list[str], config, variant_of_search: int) -> NoReturn:
    text4 = ' Creating configuration files for PTM search '
    number4 = int(round((200 - len(text4)) / 2, 0))
    logger.info(f'{text4:.^{number4}}')
    make_config_files(list_of_MOD_RES_, config, variant_of_search)
