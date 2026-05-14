from typing import NoReturn
import json
import pandas as pd
import numpy as np
import math
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
from tqdm import tqdm
from pyteomics.parser import coverage
import ast
import logging
logger = logging.getLogger("aggregate_results")

def clear_acc_and_add_alt_acc(df_psms: pd.DataFrame) -> pd.DataFrame:
    accs_list = []
    temp_rows_list = []
    for _, row in df_psms.iterrows():
        row_prots_list = ast.literal_eval(row['protein'])
        accs_list.append(row_prots_list[0].split('|')[1])

        for acc in row_prots_list[1:]:
            temp_rows_list.append(list(row) + [acc.split('|')[1]])

    df_psms['accession_of_protein'] = accs_list
    new_rows_df = pd.DataFrame(temp_rows_list, columns=list(df_psms.columns))
    return pd.concat([df_psms, new_rows_df], ignore_index=True)

from ptm_search.find_prot_name_sequence import (get_protein_sequence, get_protein_name)
# ----------------------------------------------------------------------------------------------------------------------

def get_count_peptides_for_proteins(dict_prot: dict[str, list[str]]) -> pd.DataFrame:
    count_detect_peptides_list = []
    for prot3 in dict_prot:
        try:
            seq_coverage = coverage(get_protein_sequence(prot3), list(set(dict_prot[prot3]))) * 100
        except:
            seq_coverage = 0.0

        count_detect_peptides_list.append([prot3, ";".join(dict_prot[prot3]), len(dict_prot[prot3]), seq_coverage])

    return pd.DataFrame(count_detect_peptides_list, columns=['accession', 'peptides', 'count_of_peptides', 'coverage_%'])

# ----------------------------------------------------------------------------------------------------------------------

dict_of_protein_names = {}
def give_names(accs: list[str]) -> list[str]:
    names_list = []

    for acc3 in accs:
        if acc3 in dict_of_protein_names:
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

'''
    Plotting of PTM search results
'''

def get_plots_from_result_of_analysis(SS_and_PTM_PSMs: pd.DataFrame, SS_peptides: pd.DataFrame, config, fdr_analysis_dir: Path) -> NoReturn:
    PTM_PSMs = SS_and_PTM_PSMs.query('Search != "Standard search" & modifications != "[]"')

    # ------------------------------------------------------------------------------------------------------------------
    SS_peptides = clear_acc_and_add_alt_acc(SS_peptides)
    PTM_PSMs = clear_acc_and_add_alt_acc(PTM_PSMs)
    logger.info(PTM_PSMs.shape)

    info_from_uniprot = pd.read_csv(config.ptm_search_dir / f'{config.experiment_name}_PTM_info_from_UniProt_{config.analysis_index}.csv')
    if config.search_mode == 'fast_search':
        valid_pairs = info_from_uniprot[['PTM', 'accession']].drop_duplicates().rename(columns={'PTM': 'Search', 'accession': 'accession_of_protein'})
        PTM_PSMs = PTM_PSMs.merge(valid_pairs, on=['Search', 'accession_of_protein'], how='inner')
        logger.info(PTM_PSMs.shape)

    dict_PTM_prot_peptides = PTM_PSMs.groupby('accession_of_protein')['peptide'].apply(list).to_dict()
    dict_SS_prot_peptides = SS_peptides.groupby('accession_of_protein')['peptide'].apply(list).to_dict()

    dict_PTM_and_SS_prot_peptides = {
        prot: list(set(dict_PTM_prot_peptides[prot] + dict_SS_prot_peptides.get(prot, [])))
        for prot in dict_PTM_prot_peptides}

    df_SS_prot_peptides = get_count_peptides_for_proteins(dict_SS_prot_peptides)
    df_PTM_and_SS_prot_peptides = get_count_peptides_for_proteins(dict_PTM_and_SS_prot_peptides)

    df5_prot_PTM = pd.merge(df_PTM_and_SS_prot_peptides, df_SS_prot_peptides, how="outer", on=["accession"], suffixes=('_PTM', '_FS'))

    ''' Increasing protein coverage in MSMS identification (finding new peptides for previously idefined proteins) '''
    df5_prot_PTM["count_of_peptides_delta"] = df5_prot_PTM["count_of_peptides_PTM"] - df5_prot_PTM["count_of_peptides_FS"]
    df5_prot_PTM["Protein_Name"] = give_names(df5_prot_PTM["accession"])
    # remove all proteins without a normal ID (all from the list of contaminants)
    df5_prot_PTM = df5_prot_PTM[df5_prot_PTM["accession"].notna()]
    df5_prot_PTM = df5_prot_PTM.sort_values(by=["count_of_peptides_delta"], ascending=False)
    df5_prot_PTM = df5_prot_PTM[df5_prot_PTM["count_of_peptides_delta"].notna()]
    df5_prot_PTM.to_csv(fdr_analysis_dir / f'{config.experiment_name}_increasing_coverage_of_peptides_{config.analysis_index}.csv',
                        encoding='utf-8', index=False)
    logger.info(f'\n{df5_prot_PTM.head(10)}\n')

    with sns.plotting_context('paper', font_scale=2.5), sns.axes_style('darkgrid',
                                                                       {'grid.color': '.1', 'grid.linestyle': ':'}):

        fig, ax = plt.subplots(figsize=(15, 15))
        df_plot = df5_prot_PTM.iloc[:10]

        sns.barplot(data=df_plot, y='Protein_Name', x='coverage_%_PTM', color='#f28e29',
                    label='Coverage of standard\nand PTM searches',
                    edgecolor='black', linewidth=0.9, ax=ax)

        sns.barplot(data=df_plot, y='Protein_Name', x='coverage_%_FS', color='#b51816',
                    label='Standard search\ncoverage',
                    edgecolor='black', linewidth=0.9, ax=ax)

        ax.set_ylabel("")
        ax.set_xlabel("Coverage of proteins (%)")
        ax.legend(bbox_to_anchor=(0.60, 0.25), loc='upper left')

        fig.savefig(fdr_analysis_dir / f'{config.experiment_name}_increasing_coverage_of_peptides.png',
                    dpi=300, bbox_inches='tight')
        plt.close(fig)

    # ------------------------------------------------------------------------------------------------------------------
    if config.search_mode == 'fast_search':
        filtered_msms_PTM = PTM_PSMs[PTM_PSMs['accession_of_protein'].isin(list(df5_prot_PTM['accession']))]
    else:
        filtered_msms_PTM = PTM_PSMs

    filtered_msms_PTM['Protein_Name'] = give_names(filtered_msms_PTM['accession_of_protein'])
    # temporary filter
    filtered_msms_PTM = filtered_msms_PTM[filtered_msms_PTM['Search'] != 'ADP-ribosylcysteine']

    ''' The number of PSMs, peptides and proteins found for each modification '''
    filtered_msms_PTM_psms = filtered_msms_PTM.drop_duplicates(
        subset=[c for c in filtered_msms_PTM.columns if c != 'accession_of_protein'],
        keep='first')
    filtered_msms_PTM_peptides = filtered_msms_PTM.drop_duplicates(
        subset=['accession_of_protein', 'modified_peptide', 'Search'],
        keep='first')
    filtered_msms_PTM_proteins = filtered_msms_PTM.drop_duplicates(
        subset=['Search', 'accession_of_protein'],
        keep='first')

    with sns.plotting_context("paper", font_scale=1.5), sns.axes_style('darkgrid',
                                                                        {"grid.color": ".6", "grid.linestyle": ":"}):
        fig, axes = plt.subplots(3, 1, figsize=(10, 14))
        panels = [
            (filtered_msms_PTM_psms, "Number of modified PSMs", "A"),
            (filtered_msms_PTM_peptides, "Number of modified peptides", "B"),
            (filtered_msms_PTM_proteins, "Number of modified proteins", "C")]

        for ax, (df_plot, xlabel, letter) in zip(axes, panels):
            order = df_plot['Search'].value_counts().head(15).index

            sns.countplot(data=df_plot[df_plot['Search'].isin(order)], y='Search', hue='Search', edgecolor='black',
                          linewidth=0.9, order=order, ax=ax)
            ax.set_xlabel(xlabel)
            ax.set_ylabel('')
            legend = ax.get_legend()
            if legend is not None:
                legend.remove()

            ax.annotate(letter, xy=(-0.20, 1.05), xycoords='axes fraction', fontsize=22, fontweight='bold',
                        ha='left', va='top')

        fig.tight_layout()
        fig.savefig(fdr_analysis_dir / f'{config.experiment_name}_PTM_number_of_modified_PSM_peptides_proteins.png',
                    dpi=300, bbox_inches='tight')
        plt.close(fig)

    Number_of_modified_PSMs = filtered_msms_PTM_psms.groupby(['accession_of_protein', 'Protein_Name', 'modified_peptide', 'Search']).size().to_frame('Count')
    Number_of_modified_peptides = filtered_msms_PTM_peptides.groupby(['accession_of_protein', 'Protein_Name', 'Search']).size().to_frame('Count')
    Number_of_modified_proteins = filtered_msms_PTM_proteins.groupby(['accession_of_protein', 'Protein_Name', 'Search']).size().to_frame('Count')

    # ------------------------------------------------------------------------------------------------------------------
    summary_xlsx = (fdr_analysis_dir / f'{config.experiment_name}_PTM_counts_summary_{config.analysis_index}.xlsx')

    details_xlsx = (fdr_analysis_dir / f'{config.experiment_name}_PTM_lists_detailed_{config.analysis_index}.xlsx')

    with pd.ExcelWriter(summary_xlsx, engine="xlsxwriter") as writer:

        filtered_msms_PTM_psms['Search'].value_counts(sort=True).to_frame(
            name='Count').to_excel(writer, sheet_name='PSMs_by_PTM')

        filtered_msms_PTM_peptides['Search'].value_counts(sort=True).to_frame(
            name='Count').to_excel(writer, sheet_name='Peptides_by_PTM')

        filtered_msms_PTM_proteins['Search'].value_counts(sort=True).to_frame(
            name='Count').to_excel(writer, sheet_name='Proteins_by_PTM')

    with pd.ExcelWriter(details_xlsx, engine="xlsxwriter") as writer:
        Number_of_modified_PSMs.to_excel(writer, sheet_name='PSMs_detailed')

        Number_of_modified_peptides.to_excel(writer, sheet_name='Peptides_detailed')

        Number_of_modified_proteins.to_excel(writer, sheet_name='Proteins_detailed')

    # ------------------------------------------------------------------------------------------------------------------
    ''' Saving information based on the PTM search result '''
    with open(fdr_analysis_dir / f'{config.experiment_name}_result_info.txt', 'w', encoding="utf-8") as file:
        FS_mean_coverage = df5_prot_PTM['coverage_%_FS'].mean()
        file.write(f'Average protein coverage in the standard search: {FS_mean_coverage} %\n')
        PTM_mean_coverage = df5_prot_PTM['coverage_%_PTM'].mean()
        file.write(f'Average protein coverage in standard search and PTM search: {PTM_mean_coverage} %\n\n')

        file.write(f'The number of PSMs for which PTMs were found: {filtered_msms_PTM.shape[0]}\n')
        file.write(f'The number of peptides for which PTMs were found: {filtered_msms_PTM_peptides.shape[0]}\n')
        file.write(f'The number of proteins for which PTMs were found: {df5_prot_PTM.shape[0]}\n\n')

        number_of_identified_PTMs = len(list(set(filtered_msms_PTM['Search'])))
        file.write(f'The number of PTMs with which peptides have been identified: {number_of_identified_PTMs}\n\n')
        list_of_identified_PTMs = '\n'.join(map(str, list(set(filtered_msms_PTM['Search']))))
        file.write(list_of_identified_PTMs)
    # ------------------------------------------------------------------------------------------------------------------

    df_locations_of_PTMs = PTM_PSMs[['file_name', 'spectrum', 'accession_of_protein', 'modified_peptide', 'peptide', 'modifications', 'Search']]
    df_locations_of_PTMs = df_locations_of_PTMs.drop_duplicates(subset=['file_name', 'spectrum', 'accession_of_protein', 'modified_peptide', 'peptide', 'modifications', 'Search'], keep='first')

    ''' Creating a column with positions of modified amino acids in identified peptides '''
    positions_list = []
    BACKGROUND_MODS = ('160.031@', '147.035@')
    for ptm_location in df_locations_of_PTMs['modifications']:
        ptm_location = ast.literal_eval(ptm_location)
        ptm_location = [x.split('@')[1] for x in ptm_location if not any(x.startswith(p) for p in BACKGROUND_MODS)]
        positions_list.append('_'.join(ptm_location))
    df_locations_of_PTMs['PTM_positions'] = positions_list

    df_locations_of_PTMs = df_locations_of_PTMs.query('PTM_positions != ""')

    ''' Creating a dictionary with information for each protein for each modified peptide '''
    dict_of_locations_of_PTMs = {}
    for _, row in df_locations_of_PTMs.iterrows():
        if row['accession_of_protein'] not in dict_of_locations_of_PTMs:
            dict_of_locations_of_PTMs[row['accession_of_protein']] = [[row['peptide'], row['PTM_positions'],
                                                                       row['Search'], row['modified_peptide'],
                                                                       row['file_name'], row['spectrum']]]
            continue
        dict_of_locations_of_PTMs[row['accession_of_protein']] += [[row['peptide'], row['PTM_positions'],
                                                                    row['Search'], row['modified_peptide'],
                                                                    row['file_name'], row['spectrum']]]

    ''' Creating dictionaries of protein sequences and names '''
    dict_of_full_sequences = {}
    dict_of_prot_names = {}
    for protein0 in tqdm(dict_of_locations_of_PTMs.keys()):
        dict_of_full_sequences[protein0] = get_protein_sequence(protein0)
        dict_of_prot_names[protein0] = get_protein_name(protein0)

    ''' Dictionary with short modifications names '''
    module_dir = Path(__file__).parent.parent.resolve()
    with open(module_dir / 'data' / f'ptm_name_to_config_ptm_name_dict.json', 'r') as ptm_name_to_config_ptm_name_file:
        dict_PTM_name = json.load(ptm_name_to_config_ptm_name_file)

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

    dict_PTM_color = {}
    for index, PTM in enumerate(dict_PTM_name):
        dict_PTM_color[PTM] = palette[index]

    full_result, modified_positions_list = [], []
    for protein in tqdm(dict_of_locations_of_PTMs):

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

                modified_positions_list.append([peptide[4], peptide[5], protein, peptide[0], pos, peptide[2], peptide[3]])
                full_sequence_array = (full_sequence_array[:pos-1] + [f'<span style="background-color:{dict_PTM_color[peptide[2]]};">' + dict_PTM_name[peptide[2]][-1] + '</span>'] + full_sequence_array[pos:])

        if len(temporary_dict.keys()) == 0:
            continue

        width = 60
        number = math.ceil(len(full_sequence_array) / width)
        seq_array = np.array([full_sequence_array+[' ']*(width*number - len(full_sequence_array))]).reshape(number, width)

        full_result.append(f'{protein} | {dict_of_prot_names[protein]}<br/>')
        full_result.append('<br/>'.join([''.join([str(c) for c in lst]) for lst in seq_array]))
        for temporary in set(temporary_dict.keys()):
            full_result.append(f'<br/><span style="background-color:{dict_PTM_color[temporary]};">&nbsp;</span> - {temporary} - {len(list(temporary_dict[temporary]))}: {temporary_dict[temporary]}<br/>')

        full_result.append('<br/><br/>')

    full_result = ''.join(full_result)

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
    df_of_modified_positions = pd.DataFrame(modified_positions_list, columns=['file_name', 'spectrum', 'accession',
                                                                              'peptide', 'PTM_positions', 'PTM',
                                                                              'modified_peptide'])

    positions_info_from_db = {}
    for _, row2 in info_from_uniprot.iterrows():
        if row2['accession'] not in positions_info_from_db:
            positions_info_from_db[row2['accession']] = set()
        for coordinate in row2['coordinate'].split('|'):
            positions_info_from_db[row2['accession']].add(f'{row2["PTM"]}_{coordinate}')

    match_with_bd_info = []
    for _, row3 in df_of_modified_positions.iterrows():
        if row3['accession'] not in positions_info_from_db:
            match_with_bd_info.append(False)
            continue
        if f'{row3["PTM"]}_{row3["PTM_positions"]}' in positions_info_from_db[row3['accession']]:
            match_with_bd_info.append(True)
            continue
        match_with_bd_info.append(False)
    df_of_modified_positions['match_with_bd_info'] = match_with_bd_info
    df_of_modified_positions.to_excel(fdr_analysis_dir / f'{config.experiment_name}_modified_positions_{config.analysis_index}.xlsx', index=True, header=True)
