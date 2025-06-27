'''
    Получение fasta-файла
'''
from tqdm import tqdm
import os
from random import sample
import json
from pathlib import Path
import gzip

'''
    make_fasta_file_for_searches <-- make_fasta_file
'''
# -------------------------------/ Получаем fasta-файл по accession скачивая из UniProt /-------------------------------
import requests as r
import io
import pyteomics.fasta

def make_fasta_file(proteins_from_msms_with_ptms_, ptm_groups_, fasta_file_, config):
    baseUrl="http://www.uniprot.org/uniprot/"

    module_dir = Path(__file__).parent.parent.resolve()
    with gzip.open(module_dir / 'data' / 'protein_fastas_dict.json.gz', 'r') as my_file:
        protein_fastas_dict_json = my_file.read()# Достаем информацию о последовательностях fasta-формата
    protein_fastas_dict = json.loads(protein_fastas_dict_json)

    handle = io.open(config.ptm_search_dir / f'{config.experiment_name}_PTMs_search_{config.analysis_index}_reverse.fasta', 'w', newline='\n')
    fasta_file_text = fasta_file_.read()
    handle.write(fasta_file_text)# Создаем полную базу поиска
    handle.close()

    # Получаем список идентификаторов белков всего протеома
    proteom_accessions_ = []
    for fasta_line in fasta_file_text.splitlines():
        if '>' in fasta_line:
            proteom_accessions_.append(fasta_line.split('|')[1])
    proteom_accessions_ = list(set(proteom_accessions_))
    print(f'Количество белков в изначальной базе поиска fasta-файла: {len(proteom_accessions_)}\nКоличество белков, которые будут проверяться на модификации: {len(proteins_from_msms_with_ptms_)}')

    # Загрузка 3000 последовательностей (target'ов)
    Additional_proteins_for_search = []
    # Дополнительные 3000 белков не должны быть белками для PTM проверки
    proteom_accessions_ = [p for p in proteom_accessions_ if p not in proteins_from_msms_with_ptms_]

    for accession in sample(list(set(proteom_accessions_)), 3000):
        try:
            Additional_proteins_for_search.append(protein_fastas_dict[accession])
        except KeyError:
            currentUrl1 = baseUrl + accession + ".fasta"
            response1 = r.post(currentUrl1)
            Additional_proteins_for_search.append(response1.text)
        except:
            continue
    Additional_proteins_for_search_text = ''.join(map(str, Additional_proteins_for_search))

    fasta_for_fast_search_dir = config.ptm_search_dir / f'{config.analysis_index}_fasta_for_fast_search'
    os.makedirs(fasta_for_fast_search_dir, exist_ok=True)

    for modif in tqdm(ptm_groups_.keys()):
        modif2 = modif.replace(' ', '_').replace(';', '').replace('/', '_')
        with open(fasta_for_fast_search_dir / f'{modif2}_{config.analysis_index}.fasta', 'w') as handle2:
            for protein2 in ptm_groups_[modif].keys():
                try:
                    handle2.write(protein_fastas_dict[protein2])
                except KeyError:
                    currentUrl2 = baseUrl + protein2 + ".fasta"
                    response2 = r.post(currentUrl2)
                    handle2.write(response2.text)
                except:
                    continue

        with open(fasta_for_fast_search_dir / f'{modif2}_{config.analysis_index}.fasta', 'a') as handle3: # _reverse
            # Добавление 3000 последовательностей (target'ов)
            handle3.write(Additional_proteins_for_search_text)

        # Добавление в малые базы поиска decoy-последовательности
        pyteomics.fasta.write_decoy_db(source=fasta_for_fast_search_dir / f'{modif2}_{config.analysis_index}.fasta',
                                       output=fasta_for_fast_search_dir / f'{modif2}_{config.analysis_index}_reverse.fasta',
                                       mode='reverse', prefix='DECOY_', decoy_only=False)

# ---------------------/ Открытие необходимых файлов на чтение и запись. Запуск внутренных функций /--------------------
def make_fasta_file_for_searches(proteins_from_msms_with_ptms, ptm_groups, fasta_file, config):
    text3 = ' Создание базы данных для поиска PTM '
    number3 = int(round((200 - len(text3)) / 2, 0))
    print(f'{text3:.^{number3}}')
    print(f'\nКоличество белков для которых предполагаются PTM: {len(proteins_from_msms_with_ptms)}')

    make_fasta_file(proteins_from_msms_with_ptms, ptm_groups, fasta_file, config)
