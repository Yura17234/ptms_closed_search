'''
    Скрипт выполняющий стандартную обработку результатов поиска PTM
'''

import json
import pandas as pd
import numpy as np
import math
import matplotlib.pyplot as plt
import seaborn as sns

from pathlib import Path
from tqdm import tqdm
from pyteomics.parser import coverage

def clear_acc_and_add_alt_acc(df_psms):
    new_df_psms = pd.DataFrame(columns=list(df_psms.columns) + ['accession_of_protein'])
    accs_list = []
    for index, row in df_psms.iterrows():
        accs_list.append(row['protein'].split('|')[1])

        for acc in list(row['protein'][2:-2].split("', '"))[1:]:
            new_df_psms.loc[len(new_df_psms)] = list(row) + [acc.split('|')[1]]

    df_psms['accession_of_protein'] = accs_list
    new_df_psms = pd.concat([df_psms, new_df_psms], ignore_index=True)
    return new_df_psms

from ptm_search.find_prot_name_sequence import (get_protein_sequence, get_protein_name)
# ----------------------------------------------------------------------------------------------------------------------

def get_count_peptides_for_proteins(dict_prot):
    df_count_detect_peptides = pd.DataFrame(columns=['accession', 'peptides', 'count_of_peptides', 'coverage_%'])
    for prot3 in dict_prot.keys():
        # Подсчет покрытия
        try:
            seq_coverage = coverage(get_protein_sequence(prot3), list(set(dict_prot[prot3]))) * 100
        except:  # Если последовательность не нашлась
            seq_coverage = 0.0

        df_count_detect_peptides.loc[len(df_count_detect_peptides)] = [prot3, ";".join(map(str, dict_prot[prot3])), len(dict_prot[prot3]), seq_coverage]

    return df_count_detect_peptides
# ----------------------------------------------------------------------------------------------------------------------

def give_names(accs):
    names_list = []
    dict_of_protein_names = {}

    for acc3 in accs:

        if acc3 in dict_of_protein_names.keys():
            names_list.append(dict_of_protein_names[acc3])
            continue

        name = get_protein_name(acc3)
        dict_of_protein_names[acc3] = name
        names_list.append(name)

    return names_list

# ----------------------------------------------------------------------------------------------------------------------
'''
                     /--- clear_acc0
                    /--- give_names <-- get_protein_name
    get_plots_from_result_of_analysis <-- get_count_peptides_for_proteins <-- make_dictionary_of_proteins
'''
def get_plots_from_result_of_analysis(SS_and_PTM_PSMs, SS_peptides, config, fdr_analysis_dir):
    PTM_PSMs = SS_and_PTM_PSMs.query('Search != "Standard search" & modifications != "[]"')

    # ------------------------------------------------------------------------------------------------------------------
    SS_peptides = clear_acc_and_add_alt_acc(SS_peptides)
    PTM_PSMs = clear_acc_and_add_alt_acc(PTM_PSMs)
    print(PTM_PSMs.shape)

    info_from_uniprot = pd.read_csv(config.ptm_search_dir / f'{config.experiment_name}_PTM_info_from_UniProt_{config.analysis_index}.csv')
    if config.search_mode == 'fast_search':
        new_PTM_PSMs = pd.DataFrame(columns=list(PTM_PSMs.columns))
        for ptm in list(set(PTM_PSMs['Search'])):
            temporary_PTM_PSMs = PTM_PSMs.query(f'Search == "{ptm}"')
            temporary_PTM_PSMs = temporary_PTM_PSMs[temporary_PTM_PSMs["accession_of_protein"].isin(
                list(info_from_uniprot.query(f'PTM == "{ptm}"')["accession"]))]
            new_PTM_PSMs = pd.concat([new_PTM_PSMs, temporary_PTM_PSMs], ignore_index=True)
        PTM_PSMs = new_PTM_PSMs
        print(PTM_PSMs.shape)

    dict_PTM_prot_peptides = PTM_PSMs.groupby('accession_of_protein')['peptide'].apply(list).to_dict()
    dict_prot_peptides = SS_peptides.groupby('accession_of_protein')['peptide'].apply(list).to_dict()
    dict_SS_prot_peptides = dict_prot_peptides

    dict_PTM_and_SS_prot_peptides = {}
    for prot in dict_PTM_prot_peptides.keys():
        if prot in dict_prot_peptides.keys():
            dict_PTM_and_SS_prot_peptides[prot] = list(set(dict_prot_peptides[prot] + dict_PTM_prot_peptides[prot]))
        elif prot not in dict_prot_peptides.keys():
            dict_PTM_and_SS_prot_peptides[prot] = dict_PTM_prot_peptides[prot]

    df_SS_prot_peptides = get_count_peptides_for_proteins(dict_SS_prot_peptides)
    df_PTM_and_SS_prot_peptides = get_count_peptides_for_proteins(dict_PTM_and_SS_prot_peptides)

    df5_prot_PTM = pd.merge(df_PTM_and_SS_prot_peptides, df_SS_prot_peptides, how="outer", on=["accession"], suffixes=('_PTM', '_FS'))

    '''Увеличение покрытия белков в msms идентификаии (нахождение новых пептидов для ранее определенных белков)'''
    df5_prot_PTM["count_of_peptides_delta"] = df5_prot_PTM["count_of_peptides_PTM"] - df5_prot_PTM["count_of_peptides_FS"]
    df5_prot_PTM["Protein_Name"] = give_names(df5_prot_PTM["accession"])
    # Убираем все белки без нормального ID (но по сути это наверно все из списка контоминант)
    df5_prot_PTM = df5_prot_PTM[df5_prot_PTM["accession"].notna()]
    df5_prot_PTM = df5_prot_PTM.sort_values(by=["count_of_peptides_delta"], ascending=False)
    df5_prot_PTM = df5_prot_PTM[df5_prot_PTM["count_of_peptides_delta"].notna()]
    df5_prot_PTM.to_csv(fdr_analysis_dir / f'{config.experiment_name}_increasing_coverage_of_peptides_{config.analysis_index}.csv',
                        encoding='utf-8', index=False)
    print()
    print(f'{df5_prot_PTM.head(10)}\n')


    plt.figure(figsize=(15, 15))
    sns.set(font_scale=2.5)
    sns.set_style({"grid.color": ".1", "grid.linestyle": ":"})
    plt.subplot(1, 1, 1)
    # Coverage of standard\n and PTM searches | Покрытие стандартного\nи PTM поисков
    sns.barplot(data=df5_prot_PTM[0:10], y='Protein_Name', x="coverage_%_PTM", color="#f28e29",
                label="Coverage of standard\nand PTM searches", edgecolor="black", linewidth=0.9)
    plt.legend(bbox_to_anchor=(0.60, 0.25), loc='upper left')

    plt.subplot(1, 1, 1)
    # Standard search\ncoverage | Покрытие стандартного\nпоиска
    sns.barplot(data=df5_prot_PTM[0:10], y="Protein_Name", x="coverage_%_FS", color="#b51816",
                label="Standard search\ncoverage", edgecolor="black", linewidth=0.9)
    plt.legend(bbox_to_anchor=(0.60, 0.25), loc='upper left')
    plt.ylabel("")  # Proteins
    plt.xlabel("Coverage of proteins (%)")  # Coverage of proteins | Покрытие белков
    plt.savefig(fdr_analysis_dir / f"{config.experiment_name}_increasing_coverage_of_peptides.png",
                dpi=100, bbox_inches='tight')

    # ------------------------------------------------------------------------------------------------------------------
    if config.search_mode == 'fast_search':
        filtered_msms_PTM = PTM_PSMs[PTM_PSMs["accession_of_protein"].isin(list(df5_prot_PTM["accession"]))]
    else:
        filtered_msms_PTM = PTM_PSMs

    filtered_msms_PTM["Protein_Name"] = give_names(filtered_msms_PTM["accession_of_protein"])

    filtered_msms_PTM_psms = filtered_msms_PTM.drop_duplicates(subset=[elm for elm in list(filtered_msms_PTM.columns) if elm != 'accession_of_protein'], keep='first')

    '''Число PSM найденных по каждой модификации'''
    plt.figure(figsize=(10, 10))
    sns.set(font_scale=1)
    sns.set_style({"grid.color": ".6", "grid.linestyle": ":"})
    sns.countplot(data=filtered_msms_PTM_psms, y='Search', hue='Search', edgecolor="black", linewidth=0.9,
                  order=filtered_msms_PTM_psms['Search'].value_counts().index)
    plt.xlabel("Number of modified PSMs")  # Number of modified PSMs | Число PSM с PTM
    plt.ylabel("")  # Modifications | Модификации
    plt.savefig(fdr_analysis_dir / f"{config.experiment_name}_number_of_modified_PSMs.png", dpi=100,
                bbox_inches='tight')
    filtered_msms_PTM_psms['Search'].value_counts(sort=True).to_excel(fdr_analysis_dir / f'{config.experiment_name}_number_of_modified_PSMs_all.xlsx', index=True,
                                                                 header=True)

    Number_of_modified_PSMs = filtered_msms_PTM_psms.groupby(
        ['accession_of_protein', 'Protein_Name', 'modified_peptide', 'Search']).size().to_frame('Count')
    Number_of_modified_PSMs.to_excel(fdr_analysis_dir / f'{config.experiment_name}_list_of_modified_PSMs_{config.analysis_index}.xlsx',
        index=True, header=True)

    # ------------------------------------------------------------------------------------------------------------------
    ''' Число пептидов найденных по каждой модификации '''
    filtered_msms_PTM_peptides = filtered_msms_PTM.drop_duplicates(
        subset=['accession_of_protein', 'modified_peptide', 'Search'], keep='first')
    plt.figure(figsize=(10, 10))
    sns.set(font_scale=1)
    sns.set_style({"grid.color": ".6", "grid.linestyle": ":"})
    sns.countplot(data=filtered_msms_PTM_peptides, y='Search', hue='Search', edgecolor="black", linewidth=0.9,
                  order=filtered_msms_PTM_peptides['Search'].value_counts().index)
    plt.xlabel("Number of modified peptides")  # Number of modified peptides | Число модифицированных пептидов
    plt.ylabel("")  # Modifications | Модификации
    plt.savefig(fdr_analysis_dir / f"{config.experiment_name}_number_of_modified_peptides.png", dpi=100,
                bbox_inches='tight')
    filtered_msms_PTM_peptides['Search'].value_counts(sort=True).to_excel(fdr_analysis_dir / f'{config.experiment_name}_number_of_modified_peptides_all.xlsx',
                                                                          index=True, header=True)

    Number_of_modified_peptides = filtered_msms_PTM_peptides.groupby(['accession_of_protein', 'Protein_Name', 'Search']).size().to_frame('Count')
    Number_of_modified_peptides.to_excel(fdr_analysis_dir / f'{config.experiment_name}_list_of_modified_peptides_{config.analysis_index}.xlsx', index=True, header=True)

    # ------------------------------------------------------------------------------------------------------------------
    '''Число белков найденных по каждой модификации'''
    filtered_msms_PTM_proteins = filtered_msms_PTM.drop_duplicates(subset=['Search', 'accession_of_protein'],
                                                                   keep='first')
    plt.figure(figsize=(10, 10))
    sns.set(font_scale=1)
    sns.set_style({"grid.color": ".6", "grid.linestyle": ":"})
    sns.countplot(data=filtered_msms_PTM_proteins, y='Search', hue='Search', edgecolor="black", linewidth=0.9,
                  order=filtered_msms_PTM_proteins['Search'].value_counts().index)
    plt.xlabel("Count of modified proteins")  # Count of modified proteins | Число модифицированных белков
    plt.ylabel("")  # Modifications | Модификации
    plt.savefig(fdr_analysis_dir / f"{config.experiment_name}_number_of_modified_proteins.png", dpi=100,
                bbox_inches='tight')
    filtered_msms_PTM_proteins['Search'].value_counts(sort=True).to_excel(fdr_analysis_dir / f'{config.experiment_name}_number_of_modified_proteins_{config.analysis_index}.xlsx',
                                                                          index=True, header=True)

    Number_of_modified_proteins = filtered_msms_PTM_proteins.groupby(['accession_of_protein', 'Protein_Name', 'Search']).size().to_frame('Count')
    Number_of_modified_proteins.to_excel(fdr_analysis_dir / f'{config.experiment_name}_list_of_modified_proteins_{config.analysis_index}.xlsx', index=True, header=True)

    # ------------------------------------------------------------------------------------------------------------------
    ''' Сохранение информации по результату поиска PTM '''
    with open(fdr_analysis_dir / f'{config.experiment_name}_result_info.txt', 'w', encoding="utf-8") as file:
        FS_mean_coverage = df5_prot_PTM['coverage_%_FS'].mean()
        file.write(f'Среднее покрытие белков в стандартном поиске: {FS_mean_coverage} %\n')
        PTM_mean_coverage = df5_prot_PTM['coverage_%_PTM'].mean()
        file.write(f'Среднее покрытие белков в стандартном поиске и поиске PTM: {PTM_mean_coverage} %\n\n')

        file.write(f'Количество PSM, для которых были найдены PTM: {filtered_msms_PTM.shape[0]}\n')
        file.write(f'Количество белков, для которых были найдены PTM: {df5_prot_PTM.shape[0]}\n')
        file.write(f'Количество пептидов, для которых были найдены PTM: {filtered_msms_PTM_peptides.shape[0]}\n\n')

        number_of_identified_PTMs = len(list(set(filtered_msms_PTM['Search'])))
        file.write(f'Количество PTM, с которыми были идентифицированы пептиды белков: {number_of_identified_PTMs}\n\n')
        list_of_identified_PTMs = '\n'.join(map(str, list(set(filtered_msms_PTM['Search']))))
        file.write(list_of_identified_PTMs)
    # ------------------------------------------------------------------------------------------------------------------

    df_locations_of_PTMs = PTM_PSMs[['file_name', 'spectrum', 'accession_of_protein', 'modified_peptide', 'peptide', 'modifications', 'Search']]

    df_locations_of_PTMs = df_locations_of_PTMs.drop_duplicates(subset=['file_name', 'spectrum', 'accession_of_protein', 'modified_peptide', 'peptide', 'modifications', 'Search'], keep='first')

    ''' Создание колонки с позициями модифицированных аминокислот в идентифицированных пептидах '''
    positions_list = []
    for ptm_location in df_locations_of_PTMs['modifications']:
        ptm_location = ptm_location[2:-2].split("', '")
        ptm_location = [x.split('@')[1] for x in ptm_location if "160.031@" not in x and "147.035@" not in x]
        positions_list.append('_'.join(map(str, ptm_location)))
    df_locations_of_PTMs['PTM_positions'] = positions_list

    df_locations_of_PTMs = df_locations_of_PTMs.query('PTM_positions != ""')

    ''' Создание словаря с информацией для каждого белка по каждому модифицированному пептиду '''
    dict_of_locations_of_PTMs = {}
    for index, row in df_locations_of_PTMs.iterrows():
        if row['accession_of_protein'] not in dict_of_locations_of_PTMs.keys():
            dict_of_locations_of_PTMs[row['accession_of_protein']] = [[row['peptide'], row['PTM_positions'],
                                                                       row['Search'], row['modified_peptide'],
                                                                       row['file_name'], row['spectrum']]]
            continue
        dict_of_locations_of_PTMs[row['accession_of_protein']] += [[row['peptide'], row['PTM_positions'],
                                                                    row['Search'], row['modified_peptide'],
                                                                    row['file_name'], row['spectrum']]]

    ''' Создание словарей последовательностей и названий белков '''
    dict_of_full_sequences = {}
    dict_of_prot_names = {}
    for protein0 in tqdm(dict_of_locations_of_PTMs.keys()):
        dict_of_full_sequences[protein0] = get_protein_sequence(protein0)
        dict_of_prot_names[protein0] = get_protein_name(protein0)

    ''' Словарь с укороченными названиями модификаций '''
    module_dir = Path(__file__).parent.parent.resolve()
    with open(module_dir / 'data' / f"ptm_name_to_config_ptm_name_dict.json", "r") as ptm_name_to_config_ptm_name_file:
        ptm_name_to_config_ptm_name_dict_json = ptm_name_to_config_ptm_name_file.read()
    dict_PTM_name = json.loads(ptm_name_to_config_ptm_name_dict_json)

    palette = [
        "#FF5733", "#33FF57", "#3357FF", "#FF33A5", "#A533FF", "#33FFF2", "#FFC300", "#581845",
        "#C70039", "#900C3F", "#DAF7A6", "#FFC300", "#FF5733", "#C70039", "#900C3F", "#581845",
        "#1ABC9C", "#2ECC71", "#3498DB", "#9B59B6", "#34495E", "#16A085", "#27AE60", "#2980B9",
        "#8E44AD", "#2C3E50", "#F1C40F", "#E67E22", "#E74C3C", "#ECF0F1", "#95A5A6", "#F39C12",
        "#D35400", "#C0392B", "#BDC3C7", "#7F8C8D", "#FF8C00", "#8A2BE2", "#A52A2A", "#DEB887",
        "#5F9EA0", "#7FFF00", "#D2691E", "#FF7F50", "#6495ED", "#FFF8DC", "#DC143C", "#00FFFF",
        "#00008B", "#008B8B", "#B8860B", "#A9A9A9", "#006400", "#BDB76B", "#8B008B", "#556B2F",
        "#FF8C00", "#9932CC", "#8B0000", "#E9967A", "#8FBC8F", "#483D8B", "#2F4F4F", "#00CED1",
        "#9400D3", "#FF1493", "#00BFFF", "#696969", "#1E90FF", "#B22222", "#FFFAF0", "#228B22",
        "#FF00FF", "#DCDCDC", "#F8F8FF", "#FFD700", "#DAA520", "#808080", "#008000", "#ADFF2F",
        "#F0FFF0", "#FF69B4", "#CD5C5C", "#4B0082", "#FFFFF0", "#F0E68C", "#E6E6FA", "#FFF0F5",
        "#7CFC00", "#FFFACD", "#ADD8E6", "#F08080", "#E0FFFF", "#FAFAD2", "#D3D3D3", "#90EE90"]

    l = 0
    dict_PTM_color = {}
    for PTM in dict_PTM_name.keys():
        dict_PTM_color[PTM] = palette[l]
        l += 1

    full_result = ''
    df_of_modified_positions = pd.DataFrame(columns=['file_name', 'spectrum', 'accession', 'peptide', 'PTM_positions',
                                                     'PTM', 'modified_peptide'])
    for protein in tqdm(dict_of_locations_of_PTMs.keys()):

        if dict_of_full_sequences[protein] == 'No sequence':
            continue
        full_sequence = dict_of_full_sequences[protein]
        full_sequence_array = list(full_sequence)

        temporary_dict = {}
        for peptide in dict_of_locations_of_PTMs[protein]:

            if full_sequence.rfind(peptide[0]) < 0:
                continue
            for pos in peptide[1].split('_'):
                pos = full_sequence.rfind(peptide[0]) + int(pos)

                if peptide[2] not in temporary_dict.keys():
                    temporary_dict[peptide[2]] = []
                temporary_dict[peptide[2]] = list(set(temporary_dict[peptide[2]] + [pos]))

                df_of_modified_positions.loc[len(df_of_modified_positions)] = [peptide[4], peptide[5], protein, peptide[0], pos, peptide[2], peptide[3]]
                full_sequence_array = (full_sequence_array[:pos-1] + [f'<span style="background-color:{dict_PTM_color[peptide[2]]};">' + dict_PTM_name[peptide[2]][-1] + '</span>'] + full_sequence_array[pos:])

        if len(temporary_dict.keys()) == 0:
            continue

        width = 60
        number = math.ceil(len(full_sequence_array) / width)
        seq_array = np.array([full_sequence_array+[' ']*(width*number - len(full_sequence_array))]).reshape(number, width)

        full_result += f'{protein} | {dict_of_prot_names[protein]}<br/>'
        full_result += '<br/>'.join([''.join([str(c) for c in lst]) for lst in seq_array])
        for temporary in set(temporary_dict.keys()):
            full_result += f'<br/><span style="background-color:{dict_PTM_color[temporary]};">&nbsp;</span> - {temporary} - {len(list(temporary_dict[temporary]))}: {temporary_dict[temporary]}<br/>'

        full_result += '<br/><br/>'

    html_template = f""" 
    <html> 
    <h2>{config.experiment_name}_PTM_Search</h2> 
    <body style="font-family: arial;font-size: 15px;font-family: monospace;"> 
    <p>{full_result}</p> 
      
    </body> 
    </html> 
    """
    f = open(fdr_analysis_dir / f'PTM_localization.html', 'w')
    f.write(html_template)
    f.close()

    positions_info_from_db = {}
    for index2, row2 in info_from_uniprot.iterrows():
        if row2['accession'] not in positions_info_from_db.keys():
            positions_info_from_db[row2['accession']] = []
        for coordinate in row2['coordinate'].split('|'):
            positions_info_from_db[row2['accession']].append(f'{row2["PTM"]}_{coordinate}')

    match_with_bd_info = []
    for index3, row3 in df_of_modified_positions.iterrows():
        if row3['accession'] not in positions_info_from_db.keys():
            match_with_bd_info.append(False)
            continue
        if f'{row3["PTM"]}_{row3["PTM_positions"]}' in positions_info_from_db[row3['accession']]:
            match_with_bd_info.append(True)
            continue
        match_with_bd_info.append(False)
    df_of_modified_positions['match_with_bd_info'] = match_with_bd_info
    df_of_modified_positions.to_excel(fdr_analysis_dir / f'{config.experiment_name}_modified_positions_{config.analysis_index}.xlsx', index=True, header=True)