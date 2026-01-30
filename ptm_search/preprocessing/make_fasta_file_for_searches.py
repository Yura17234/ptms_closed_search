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
import pyteomics.fasta

'''
    Getting a fasta files
'''

def make_fasta_file(proteins_from_msms_with_ptms_: set[str], ptm_groups_: dict[str, dict[str, list[str]]], fasta_file_: TextIO, config) -> NoReturn:
    baseUrl="http://www.uniprot.org/uniprot/"

    module_dir = Path(__file__).parent.parent.resolve()
    with gzip.open(module_dir / 'data' / 'protein_fastas_dict.json.gz', 'rt', encoding='utf-8') as my_file:
        protein_fastas_dict = json.load(my_file)

    # get a list of protein identifiers for the entire proteome
    proteom_accessions_ = set()
    out_path = config.ptm_search_dir / f'{config.experiment_name}_PTMs_search_{config.analysis_index}_reverse.fasta'
    with open(out_path, 'w') as handle:
        for fasta_line in fasta_file_:
            handle.write(fasta_line)
            if fasta_line.startswith('>'):
                proteom_accessions_.add(fasta_line.split('|')[1])
    logger.info(f'The number of proteins in the original fasta file search database: {len(proteom_accessions_)}\nThe number of proteins that will be tested for modification: {len(proteins_from_msms_with_ptms_)}')

    # loading 5000 target sequences
    Additional_proteins_for_search = []
    # the additional 5,000 proteins do not need to be proteins for PTM identification
    proteom_accessions_ = [p for p in proteom_accessions_ if p not in proteins_from_msms_with_ptms_]

    session = r.Session()
    timeout = (5, 15)
    cache_hits = 0
    cache_misses = 0
    http_ok = 0
    http_fail = 0

    for accession in sample(proteom_accessions_, 5000):
        try:
            Additional_proteins_for_search.append(protein_fastas_dict[accession])
            cache_hits += 1
        except KeyError:
            cache_misses += 1
            currentUrl1 = baseUrl + accession + '.fasta'
            try:
                response1 = session.get(currentUrl1, timeout=timeout)
                response1.raise_for_status()
                Additional_proteins_for_search.append(response1.text)
                http_ok += 1
            except r.exceptions.RequestException as e:
                http_fail += 1
                logger.debug(f'Failed to fetch {accession}: {e}')
                continue
    Additional_proteins_for_search_text = ''.join(Additional_proteins_for_search)

    fasta_for_fast_search_dir = config.ptm_search_dir / f'{config.analysis_index}_fasta_for_fast_search'
    os.makedirs(fasta_for_fast_search_dir, exist_ok=True)

    for modif in tqdm(ptm_groups_):
        modif2 = modif.replace(' ', '_').replace(';', '').replace('/', '_')
        with open(fasta_for_fast_search_dir / f'{modif2}_{config.analysis_index}.fasta', 'w') as handle2:
            for protein2 in ptm_groups_[modif]:
                try:
                    handle2.write(protein_fastas_dict[protein2])
                    cache_hits += 1
                except KeyError:
                    cache_misses += 1
                    currentUrl2 = baseUrl + protein2 + '.fasta'
                    try:
                        response2 = session.get(currentUrl2, timeout=timeout)
                        response2.raise_for_status()
                        handle2.write(response2.text)
                        http_ok += 1
                    except r.exceptions.RequestException as e:
                        http_fail += 1
                        logger.debug(f'Failed to fetch {protein2}: {e}')
                        continue

        if len(ptm_groups_[modif]) < 5000:
            with open(fasta_for_fast_search_dir / f'{modif2}_{config.analysis_index}.fasta', 'a') as handle3: # _reverse
                # adding 5000 target sequences
                handle3.write(Additional_proteins_for_search_text)

        # adding decoy sequences to small search databases
        pyteomics.fasta.write_decoy_db(source=str(fasta_for_fast_search_dir / f'{modif2}_{config.analysis_index}.fasta'),
                                       output=str(fasta_for_fast_search_dir / f'{modif2}_{config.analysis_index}_reverse.fasta'),
                                       mode='reverse', prefix='DECOY_', decoy_only=False)

    logger.info(
        f'FASTA fetch stats: '
        f'cache_hits={cache_hits}, '
        f'cache_misses={cache_misses}, '
        f'http_ok={http_ok}, '
        f'http_fail={http_fail}')

# -----------------/ Opening the necessary files for reading and writing. Launching internal functions /----------------
def make_fasta_file_for_searches(proteins_from_msms_with_ptms: set[str], ptm_groups: dict[str, dict[str, list[str]]], fasta_file: TextIO, config) -> NoReturn:
    text3 = ' Creating search databases for PTM searches '
    number3 = int(round((200 - len(text3)) / 2, 0))
    logger.info(f'{text3:.^{number3}}')
    logger.info(f'\nThe number of proteins for which PTMs are assumed: {len(proteins_from_msms_with_ptms)}')

    make_fasta_file(proteins_from_msms_with_ptms, ptm_groups, fasta_file, config)
