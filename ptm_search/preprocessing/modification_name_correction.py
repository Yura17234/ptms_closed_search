import json
from pathlib import Path
# ------------------------------/ Функция парсит в строке записи UniProt'а название белка /-----------------------------

module_dir = Path(__file__).parent.parent.resolve()
with open(module_dir / 'data' / 'ptm_keyword_map_dict.json', 'r') as my_file:
    ptm_keyword_map_dict = json.load(my_file)

def smaller_groups(modif: str) -> str | None:
    modif_lower = modif.lower()

    for keyword, result in ptm_keyword_map_dict.items():
        if keyword in modif_lower:
            return result
    return None
