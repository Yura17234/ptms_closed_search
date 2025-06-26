'''
    Создание mgf-файлов только с не идентифицированными спектрами для PTM-поиска
'''
import pyteomics.mgf
from tqdm import tqdm
import pandas as pd
import os
from multiprocessing import Pool
import params

# Информация обо всех и идентифицированных спектрах представляется как глобальные переменные
def init_pool(union_PSMs_df0):
    global union_PSMs_df
    union_PSMs_df = union_PSMs_df0

# Осуществление записи списка спектров в новый mgf-файл
def make_mgf_files_for_ptm(file_name0):

    dict_mgf0 = pyteomics.mgf.IndexedMGF(f'{params.path}/{params.x2}_Standard_search/full_mgf_files/{file_name0.split(".")[0]}.mgf')
    temporary_spectra_list0 = union_PSMs_df.query(f'file == "{file_name0.split(".")[0]}.pep.xml"')['spectrum']

    outmgf = open(f'{params.path}/{file_name0.split(".")[0]}_for_PTM.mgf', 'w')
    for spectrum in dict_mgf0:

        if spectrum['params']['title'] not in list(temporary_spectra_list0):

            outmgf.write('BEGIN IONS\n')
            outmgf.write('TITLE=%s\n' % (spectrum['params']['title']))
            outmgf.write(f'PEPMASS={spectrum["params"]["pepmass"][0]:.{len(str(spectrum["params"]["pepmass"][0]).split(".")[1])}f}\n')
            outmgf.write(f'RTINSECONDS={float(spectrum["params"]["rtinseconds"]):.{len(str(spectrum["params"]["rtinseconds"]).split(".")[1])}f}\n')
            outmgf.write('CHARGE=%d+\n' % (spectrum['params']['charge'][0],))
            outmgf.write('SCANS=%s\n' % (spectrum['params']['scans']))

            for m_z, intensity in zip(spectrum['m/z array'], spectrum['intensity array']):
                outmgf.write(f'{m_z:.{len(str(m_z).split(".")[1])}f} {intensity:.{len(str(intensity).split(".")[1])}f} \n')
            outmgf.write('END IONS\n\n')
    outmgf.close()

    return file_name0

# Использование параллелизации записи списка спектров в новый mgf-файл
def make_mgf_files_for_ptm_multiprocessing(mgf_files_list0, union_PSMs_df0):

    with Pool(initializer=init_pool, initargs=(union_PSMs_df0,), processes=8) as p:
        results = p.imap_unordered(make_mgf_files_for_ptm, mgf_files_list0)

        for filename in results:
            print(f'{filename.split(".")[0]}_for_PTM.mgf --> Done!\n')

# =====================/ Открытие необходимых файлов на чтение и запись. Запуск внутренных функций /====================
def main5(path, x2):
    text5 = ' Создание mgf-файлов только с не идентифицированными спектрами для PTM-поиска '
    number5 = int(round((200 - len(text5)) / 2, 0))
    print(f'{text5:.^{number5}}')
    union_PSMs_df = pd.read_csv(f'{path}/{x2}_Standard_search/union_PSMs.tsv', sep='\t')
    mgf_files_list = [file for file in os.listdir(f'{path}/{x2}_Standard_search/full_mgf_files') if '.mgf' in file]
    make_mgf_files_for_ptm_multiprocessing(mgf_files_list, union_PSMs_df)