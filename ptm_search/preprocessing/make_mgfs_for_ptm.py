from typing import NoReturn
import pyteomics.mgf
import pandas as pd
import os
from multiprocessing import Pool
import logging
logger = logging.getLogger("prepare_ptm_search")

'''
    Creating mgf files only with unidentified spectra for PTM search
'''

# info about all and identified spectra and getting as global variables
def init_pool(union_PSMs_df0, mgf_dir0, config0):
    global union_PSMs_df
    union_PSMs_df = union_PSMs_df0
    global mgf_dir
    mgf_dir = mgf_dir0
    global config
    config = config0

# writing a list of spectra in a new mgf file
def make_mgf_files_for_ptm(file_name0) -> NoReturn:
    fname = file_name0.split(".")[0]
    mgf_file = pyteomics.mgf.MGF(str(mgf_dir / f'{fname}.mgf'))
    temporary_spectra_set = set(union_PSMs_df[union_PSMs_df['file'] == f'{fname}.pep.xml']['spectrum'])

    outmgf = open(config.work_dir / f'{file_name0.split(".")[0]}_for_PTM.mgf', 'w')
    for spectrum in mgf_file:

        if spectrum['params']['title'] not in list(temporary_spectra_set):

            outmgf.write('BEGIN IONS\n')
            outmgf.write('TITLE=%s\n' % (spectrum['params']['title']))
            outmgf.write(f'PEPMASS={spectrum["params"]["pepmass"][0]:.{len(str(spectrum["params"]["pepmass"][0]).split(".")[1])}f}\n')
            outmgf.write(f'RTINSECONDS={float(spectrum["params"]["rtinseconds"]):.{len(str(spectrum["params"]["rtinseconds"]).split(".")[1])}f}\n')
            try:
                outmgf.write('CHARGE=%d+\n' % (spectrum['params']['charge'][0],))
            except:
                pass
            outmgf.write('SCANS=%s\n' % (spectrum['params']['scans']))

            for m_z, intensity in zip(spectrum['m/z array'], spectrum['intensity array']):
                outmgf.write(f'{m_z:.{len(str(m_z).split(".")[1])}f} {intensity:.{len(str(intensity).split(".")[1])}f} \n')
            outmgf.write('END IONS\n\n')
    outmgf.close()

    return file_name0

# using parallelization to write a list of spectra to a new mgf file
def make_mgf_files_for_ptm_multiprocessing(mgf_files_list0: list[str], union_PSMs_df0: pd.DataFrame, mgf_dir0: str, config0) -> NoReturn:

    with Pool(initializer=init_pool, initargs=(union_PSMs_df0, mgf_dir0, config0,), processes=8) as p:
        results = p.imap_unordered(make_mgf_files_for_ptm, mgf_files_list0)

        for filename in results:
            logger.info(f'{filename.split(".")[0]}_for_PTM.mgf --> Done!\n')

# -----------------/ Opening the necessary files for reading and writing. Launching internal functions /----------------
def make_mgfs_for_ptm(mgf_dir: str, config) -> NoReturn:
    text5 = ' Creating mgf files only with unidentified spectra for PTM search '
    number5 = int(round((200 - len(text5)) / 2, 0))
    logger.info(f'\n{text5:.^{number5}}')
    union_PSMs_df = pd.read_csv(config.st_search_dir / 'union_PSMs.tsv', sep='\t',
                                   usecols=['file', 'spectrum'], dtype={'file': str, 'spectrum': str})
    mgf_files_list = [file for file in os.listdir(mgf_dir) if '.mgf' in file]
    make_mgf_files_for_ptm_multiprocessing(mgf_files_list, union_PSMs_df, mgf_dir, config)
