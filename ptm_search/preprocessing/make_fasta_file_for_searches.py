from typing import NoReturn, TextIO
from tqdm import tqdm
import os
from random import sample
import json
from pathlib import Path
import gzip
import logging
logger = logging.getLogger("prepare_ptm_search")

'''
    make_fasta_file_for_searches <-- make_fasta_file
'''
# --------------------------/ Downloading fasta format sequences by accessions from UniProt /---------------------------
import requests as r
import io
import pyteomics.fasta

'''
    Getting a fasta files
'''

def make_fasta_file(proteins_from_msms_with_ptms_: list[str], ptm_groups_: dict[str, dict[str, list[str]]], fasta_file_: TextIO, config) -> NoReturn:
    baseUrl="http://www.uniprot.org/uniprot/"

    module_dir = Path(__file__).parent.parent.resolve()
    with gzip.open(module_dir / 'data' / 'protein_fastas_dict.json.gz', 'r') as my_file:
        protein_fastas_dict_json = my_file.read()  # get info about fasta format sequences
    protein_fastas_dict = json.loads(protein_fastas_dict_json)

    handle = io.open(config.ptm_search_dir / f'{config.experiment_name}_PTMs_search_{config.analysis_index}_reverse.fasta', 'w', newline='\n')
    fasta_file_text = fasta_file_.read()
    handle.write(fasta_file_text)  # creating a complete search database
    handle.close()

    # get a list of protein identifiers for the entire proteome
    proteom_accessions_ = []
    for fasta_line in fasta_file_text.splitlines():
        if '>' in fasta_line:
            proteom_accessions_.append(fasta_line.split('|')[1])
    proteom_accessions_ = list(set(proteom_accessions_))
    logger.info(f'The number of proteins in the original fasta file search database: {len(proteom_accessions_)}\nThe number of proteins that will be tested for modification: {len(proteins_from_msms_with_ptms_)}')

    # loading 5000 target sequences
    Additional_proteins_for_search = []
    # the additional 5,000 proteins do not need to be proteins for PTM identification
    proteom_accessions_ = [p for p in proteom_accessions_ if p not in proteins_from_msms_with_ptms_]

    for accession in sample(list(set(proteom_accessions_)), 5000):
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

        if len(ptm_groups_[modif].keys()) < 5000:
            with open(fasta_for_fast_search_dir / f'{modif2}_{config.analysis_index}.fasta', 'a') as handle3: # _reverse
                # adding 5000 target sequences
                handle3.write(Additional_proteins_for_search_text)

        # adding decoy sequences to small search databases
        pyteomics.fasta.write_decoy_db(source=str(fasta_for_fast_search_dir / f'{modif2}_{config.analysis_index}.fasta'),
                                       output=str(fasta_for_fast_search_dir / f'{modif2}_{config.analysis_index}_reverse.fasta'),
                                       mode='reverse', prefix='DECOY_', decoy_only=False)

# -----------------/ Opening the necessary files for reading and writing. Launching internal functions /----------------
def make_fasta_file_for_searches(proteins_from_msms_with_ptms: list[str], ptm_groups: dict[str, dict[str, list[str]]], fasta_file: TextIO, config) -> NoReturn:
    text3 = ' Creating search databases for PTM searches '
    number3 = int(round((200 - len(text3)) / 2, 0))
    logger.info(f'{text3:.^{number3}}')
    logger.info(f'\nThe number of proteins for which PTMs are assumed: {len(proteins_from_msms_with_ptms)}')

    make_fasta_file(proteins_from_msms_with_ptms, ptm_groups, fasta_file, config)
