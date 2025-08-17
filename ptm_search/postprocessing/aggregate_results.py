import pandas as pd
import numpy as np
import ast

from ptm_search.postprocessing.fdr_filtration import (
    threshold_calculation_identipy,
    threshold_calculation_for_PTM_by_ranks,
)
from ptm_search.postprocessing.analysis_of_result import get_plots_from_result_of_analysis

pd.options.mode.chained_assignment = None

def clear_ptm_psms(df_psms, list_accs_ptm_info):
    acc_was_in_ptm_info = []
    for list_accs in df_psms['accessions_list']:
        if any(acc in list_accs_ptm_info for acc in list_accs.split('|')):
            acc_was_in_ptm_info.append(True)
        else:
            acc_was_in_ptm_info.append(False)

    df_psms['acc_was_in_ptm_info'] = acc_was_in_ptm_info
    df_psms = df_psms.query('acc_was_in_ptm_info == True')
    df_psms = df_psms.drop(columns=['acc_was_in_ptm_info'])
    return df_psms

def mark_variable_modifications(modifications_column_PSMs):
    list_of_markers_of_variable = []
    for line in modifications_column_PSMs:
        result = '-'
        for line2 in line[2:-2].split("', '"):
            if "160.031@" not in line2 and "147.035@" not in line2:
                result = '+'
                break
        list_of_markers_of_variable.append(result)
    return list_of_markers_of_variable

def mark_decoys_and_targets(psms_description):
    return ["'sp" not in n for n in psms_description]

def calculate_threshold(decoys, targets, log_file, config, ptm_name, log_dir):
    if config.fdr_strategy == 'transferred_fdr':
        return threshold_calculation_for_PTM_by_ranks(
            decoys[['PTM', 'rank']],
            targets[['PTM', 'rank']],
            log_dir,
            log_file,
            config,
            ptm_name)
    else:
        return threshold_calculation_identipy(
            decoys[['PTM', 'log_hyperscore']].query("PTM == '+'"),
            targets[['PTM', 'log_hyperscore']].query("PTM == '+'"),
            log_file)

def aggregate_results(config):
    '''
        Aggregate results
    '''

    fdr_result_dir = config.ptm_search_dir / f'Result_of_PTM_search_{config.search_mode}_{config.fdr_strategy}'
    fdr_result_dir.mkdir(exist_ok=True)

    fdr_result_file_path = fdr_result_dir / f'result_PTM_{config.analysis_index}.csv'

    if not fdr_result_file_path.exists():
        log_dir = config.ptm_search_dir / f'log_info_of_filtration_{config.search_mode}_{config.fdr_strategy}'
        log_dir.mkdir(exist_ok=True)
        log_file_path = log_dir / 'log_filtration_info.txt'
        log_file = open(log_file_path, 'w')

        ss_psms = pd.read_csv(config.st_search_dir / 'union_PSMs_full.tsv', sep='\t')
        ss_psms['Search'] = 'Standard search'
        ss_psms = ss_psms.query("modifications == '[]'")
        ss_psms['PTM'] = '-'
        ss_psms['variable'] = '-'

        ptm_uniprot_info_df = pd.read_csv(config.ptm_search_dir / f"{config.experiment_name}_PTM_info_from_UniProt_{config.analysis_index}.csv")

        all_fdr_ptm_psms = pd.DataFrame()
        raw_result_dir = config.ptm_search_dir / f'{config.analysis_index}_result_{config.search_mode}'
        for index, ptm_file_path in enumerate(raw_result_dir.glob('*_result.csv')):
            ptm_name = ptm_file_path.stem.split(config.analysis_index)[0][:-1].replace('_', ' ')
            print(f'{index} / {len(list(raw_result_dir.iterdir()))} | {ptm_name}')
            if not ptm_file_path.exists():
                continue

            ptm_df = pd.read_csv(ptm_file_path)
            ptm_df['Search'] = ptm_name
            ptm_df = ptm_df.query("modifications != '[]'")
            ptm_df['variable'] = mark_variable_modifications(ptm_df['modifications'])
            ptm_df = ptm_df.query("variable == '+'")
            ptm_df['decoy'] = mark_decoys_and_targets(ptm_df['protein'])
            ptm_df['PTM'] = '+'

            full_df = pd.concat([ss_psms, ptm_df], ignore_index=True).sort_values("hyperscore")
            full_df['rank'] = range(1, len(full_df) + 1)

            target = full_df.query("decoy == False")# & PTM == '+'")
            decoy = full_df.query("decoy == True")# & PTM == '+'")

            try:
                threshold, q_values = calculate_threshold(decoy, target, log_file, config, ptm_name, log_dir)
            except:
                print(f'Размер результата анализа после фильтрации по {ptm_name} : {0}')
                continue

            ptm_df = full_df.query(f"PTM == '+' & rank >= {threshold}")
            for rank_val, q_val in sorted(q_values.items()):
                ptm_df.loc[ptm_df['rank'] > rank_val, 'q_value'] = q_val
            print(f'Размер результата анализа после фильтрации по {ptm_name} : {ptm_df.shape}')
            log_file.write(f'Размер результата анализа после фильтрации по {ptm_name} : {ptm_df.shape}\n')

            ptm_df['accessions_list'] = ptm_df['protein'].apply(lambda x: '|'.join([p.split('|')[1] for p in ast.literal_eval(x)]))

            if config.search_mode == 'fast_search':
                valid_accs = ptm_uniprot_info_df.query(f'PTM == "{ptm_name}"')['accession'].tolist()
                ptm_df = clear_ptm_psms(ptm_df, valid_accs)
                print(f'Размер результата анализа после фильтрации на наличие информации в UniProt по {ptm_name} : {ptm_df.shape}')
                log_file.write(f'Размер результата анализа после фильтрации на наличие информации в UniProt по {ptm_name} : {ptm_df.shape}\n\n')

            all_fdr_ptm_psms = pd.concat([all_fdr_ptm_psms, ptm_df], ignore_index=True)

        all_fdr_ptm_psms.to_csv(fdr_result_file_path, index=False)

        log_file.write(str(all_fdr_ptm_psms.shape))
        log_file.close()
    else:
        all_fdr_ptm_psms = pd.read_csv(fdr_result_file_path)

    fdr_analysis_dir = config.ptm_search_dir / f"Result_of_PTM_analysis_{config.search_mode}_{config.fdr_strategy}"
    fdr_analysis_dir.mkdir(exist_ok=True)

    ss_peptides = pd.read_csv(config.st_search_dir / 'union_peptides.tsv', sep='\t')
    get_plots_from_result_of_analysis(all_fdr_ptm_psms, ss_peptides, config, fdr_analysis_dir)

    print('Aggregate results -- Done !')
