'''
    Get protein name and sequence
'''
import gzip
import json
from unipressed import UniprotkbClient
from pathlib import Path

module_dir = Path(__file__).parent.resolve()
with open(module_dir / 'data' / 'protein_names_dict.json', 'r', encoding='utf-8') as my_file0:
    protein_names_dict = json.load(my_file0)

def get_protein_name(accession: str) -> str:
    try:
        name = protein_names_dict[accession]
        return name
    except KeyError:
        try:
            dict_query = UniprotkbClient.fetch_one(accession)
            return f"{dict_query['proteinDescription']['recommendedName']['fullName']['value']}"
        except:
            return 'Name not found !'

# ----------------------------------------------------------------------------------------------------------------------

with gzip.open(module_dir / 'data' / 'protein_sequences_dict.json.gz', 'rt', encoding='utf-8') as my_file1:
    protein_sequences_dict = json.load(my_file1)

def get_protein_sequence(accession1: str) -> str:
    try:
        name1 = protein_sequences_dict[accession1]
        return name1
    except KeyError:
        try:
            dict_query1 = UniprotkbClient.fetch_one(accession1)
            return f"{dict_query1['sequence']['value']}"
        except:
            return 'Sequence not found !'
