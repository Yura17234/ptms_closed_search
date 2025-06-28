'''
    Код выполняет поиск модификаций и
    замен в белковой последовательности по каждому найденному белку
'''

from tqdm import tqdm
import random
random.seed(42)
from ptm_search.preprocessing.modification_name_correction import smaller_groups

'''
                    /--- smaller_groups
    parsing_human_proteom <-- get_PTMs_lists <-- get_all_accs_and_names
'''

from ptm_search.find_prot_name_sequence import get_protein_name

# -------------------------------/ Функция создает словарь accession - название белка /---------------------------------
def get_all_accs_and_names(df):
    accs_and_names_dict = {}
    for list_accs in tqdm(df["dbname"]):
        for acc in list(list_accs[2:-2].split("', '")):
            if acc.split('|')[1] in accs_and_names_dict.keys():
                continue
            accs_and_names_dict[acc.split('|')[1]] = get_protein_name(acc.split('|')[1])

    print(f'Вид словаря accession - название белка:')
    print( '\n'.join(map(str, [ f'{dict_elm} : {accs_and_names_dict[dict_elm]}' for dict_elm in random.sample(list(accs_and_names_dict.keys()), 5)] )) )
    return accs_and_names_dict

# ---------------------------------/ Функция формирует списки белков по каждой PTM /------------------------------------
def get_PTMs_lists(query_text, modres_file, list_of_msms_proteins):
    PTMs_groups = {}
    query_text_list = query_text.split('\n//\n')
    count = 0
    list_of_missed_modifs = []
    for description in query_text_list:
        if 'MOD_RES' not in description:
            continue
        count += 1
        description_list = description.split('\n')
        acc_of_protein = ''
        for i in range(0, len(description_list)):
            if 'ID   ' in description_list[i]:
                continue
            if 'AC   ' in description_list[i]:
                for elm in description_list[i].split()[1:]:
                    if elm[:-1] in list_of_msms_proteins.keys():
                        acc_of_protein = elm[:-1]
                        modres_file.write('// ' + f'{acc_of_protein}|{list_of_msms_proteins[acc_of_protein]}' + '\n')
                        modres_file.write('\n')
                        break
            if acc_of_protein == '':
                break  # Если ни один acc не был найден в эксперименте, переходим к следующему Query
            if 'FT   MOD_RES' in description_list[i]:
                name_of_modification = ''
                # Этот цикл нужен, чтоб собрать всё название модификации белка из строки
                for u in range(28, len(description_list[i + 1])):
                    if description_list[i + 1][u] != '\"':
                        name_of_modification += description_list[i + 1][u]
                    else:
                        break
                correct_prot_modif_name = smaller_groups(name_of_modification)
                if correct_prot_modif_name == None:
                    list_of_missed_modifs.append(name_of_modification)
                    continue
                if ':' in description_list[i].split()[2]:
                    position = description_list[i].split()[2].split(':')[1]
                else:
                    position = description_list[i].split()[2]

                # Добавление информации в словарь
                if correct_prot_modif_name not in PTMs_groups.keys():
                    PTMs_groups[correct_prot_modif_name] = {}

                if acc_of_protein not in PTMs_groups[correct_prot_modif_name].keys():
                    # основной словарь с координатами модификации
                    PTMs_groups[correct_prot_modif_name][acc_of_protein] = [position]
                    modres_file.write(correct_prot_modif_name + '\n')
                    modres_file.write('\n')
                    continue

                # если позиции нет в списке
                if position not in PTMs_groups[correct_prot_modif_name][acc_of_protein]:
                    # основной словарь с координатами модификации (добавление координаты)
                    PTMs_groups[correct_prot_modif_name][acc_of_protein] += [position]
                modres_file.write(correct_prot_modif_name + '\n')
                modres_file.write('\n')

    print(f'Число антотаций белков: {count}')
    modres_file.close()

    print(f'Словарь модификаций и их белков:')
    print('\n'.join(map(str, [f'{dict_elm} : {[*PTMs_groups[dict_elm].keys()][0:5]}' for dict_elm in random.sample(PTMs_groups.keys(), 5)])))
    print(f'Число модификаций, которые не узнал встроенный словарь: {len(set(list_of_missed_modifs))}\n')
    print(f'Список модификаций, которые не узнал встроенный словарь: {set(list_of_missed_modifs)}\n')
    return PTMs_groups

# ---------------------/ Открытие необходимых файлов на чтение и запись. Запуск внутренных функций /--------------------
def parsing_human_proteom(config, dataframe):
    text1 = ' Парсинг протеома в UniProt для поиска ПТМ '
    number1 = int(round((200 - len(text1)) / 2, 0))
    print(f'{text1:.^{number1}}')
    # Протеом человека из UniProt
    with open(config.uniprot_query_path, 'r') as query_file:
        query_text = query_file.read()
    print(config.uniprot_query_path)

    # Запись модификаций белков в текстовый файл
    modres_file = open(config.ptm_search_dir / f'{config.experiment_name}_MOD_RES_{config.analysis_index}.txt', 'w')
    print(config.ptm_search_dir / f'{config.experiment_name}_MOD_RES_{config.analysis_index}.txt')

    prot_acc_and_names = get_all_accs_and_names(dataframe)
    return get_PTMs_lists(query_text, modres_file, prot_acc_and_names), prot_acc_and_names
