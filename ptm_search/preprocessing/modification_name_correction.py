import json
from pathlib import Path
import os
# ------------------------------/ Функция парсит в строке записи UniProt'а название белка /-----------------------------

def smaller_groups(modif: str) -> str | None:
    modif_lower = modif.lower()

    module_dir = Path(__file__).parent.parent.resolve()
    with open(module_dir / 'data' / 'ptm_keyword_map_dict.json', 'r') as my_file:
        ptm_keyword_map_dict_json = my_file.read()
    ptm_keyword_map_dict = json.loads(ptm_keyword_map_dict_json)

    for keyword, result in ptm_keyword_map_dict.items():
        if keyword in modif_lower:
            return result
    return None
