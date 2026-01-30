from typing import NoReturn
import pandas as pd
from pathlib import Path
import ast
import logging
logger = logging.getLogger("aggregate_results")

from ptm_search.postprocessing.fdr_filtration import (
    threshold_calculation_separate_fdr,
    threshold_calculation_transferred_fdr,
    threshold_calculation_transferred_fdr_linear_reg)

from ptm_search.postprocessing.analysis_of_result import get_plots_from_result_of_analysis

pd.options.mode.chained_assignment = None

def clear_ptm_psms(df_psms: pd.DataFrame, list_accs_ptm_info: set[str]) -> pd.DataFrame:
    mask = df_psms["accessions_list"].str.split("|").apply(lambda xs: any(a in list_accs_ptm_info for a in xs))
    return df_psms.loc[mask].copy()

BACKGROUND_MODS = ('160.031@', '147.035@')
def mark_variable_modifications(modifications_desc: str) -> bool:
    try:
        mods = ast.literal_eval(modifications_desc)
    except Exception:
        return False

    for mod in mods:
        if not any(mod.startswith(p) for p in BACKGROUND_MODS):
            return True
    return False

def mark_decoys_and_targets(psm_description: str) -> bool:
    return ("'sp|" not in psm_description) and ("'tr|" not in psm_description)

def calculate_threshold(decoys: pd.DataFrame, targets: pd.DataFrame, config, ptm_name: str, ratio_info_dir: Path) -> tuple[float, dict[float, float]]:
    if config.fdr_strategy == 'transferred_fdr':
        return threshold_calculation_transferred_fdr(
            decoys[['PTM', 'rank']],
            targets[['PTM', 'rank']],
            ratio_info_dir,
            ptm_name)
    elif config.fdr_strategy == 'transferred_fdr_linear_reg':
        return threshold_calculation_transferred_fdr_linear_reg(
            decoys[['PTM', 'rank']],
            targets[['PTM', 'rank']],
            ratio_info_dir,
            ptm_name)
    else:
        return threshold_calculation_separate_fdr(
            decoys[['PTM', 'rank']],
            targets[['PTM', 'rank']])

def aggregate_results(config) -> NoReturn:
    '''
        Aggregate results
    '''

    fdr_result_dir = config.ptm_search_dir / f'Result_of_PTM_search_{config.search_mode}_{config.fdr_strategy}'
    fdr_result_dir.mkdir(exist_ok=True)

    fdr_result_file_path = fdr_result_dir / f'result_PTM_{config.analysis_index}.csv'

    if not fdr_result_file_path.exists():
        ratio_info_dir = config.ptm_search_dir / f'ratio_info_of_filtration_{config.search_mode}_{config.fdr_strategy}'
        ratio_info_dir.mkdir(exist_ok=True)

        ss_psms = pd.read_csv(config.st_search_dir / 'union_PSMs_full.tsv', sep='\t')
        ss_psms['Search'] = 'Standard search'
        ss_psms = ss_psms.query("modifications == '[]'")
        ss_psms['PTM'] = False
        ss_psms['variable'] = False

        ptm_uniprot_info_df = pd.read_csv(config.ptm_search_dir / f"{config.experiment_name}_PTM_info_from_UniProt_{config.analysis_index}.csv")

        all_fdr_ptm_psms = pd.DataFrame()
        raw_result_dir = config.ptm_search_dir / f'{config.analysis_index}_result_{config.search_mode}'
        ptm_files = list(raw_result_dir.glob("*_result.csv"))
        for index, ptm_file_path in enumerate(ptm_files, start=1):
            ptm_name = ptm_file_path.stem.split(config.analysis_index)[0][:-1].replace('_', ' ')
            logger.info(f'{index} / {len(ptm_files)} | {ptm_name}')
            if not ptm_file_path.exists():
                continue

            ptm_df = pd.read_csv(ptm_file_path)
            ptm_df['Search'] = ptm_name
            ptm_df = ptm_df.query("modifications != '[]'")
            ptm_df['variable'] = ptm_df['modifications'].apply(mark_variable_modifications)
            ptm_df = ptm_df.query("variable == True")
            ptm_df['decoy'] = ptm_df['protein'].apply(mark_decoys_and_targets)
            ptm_df['PTM'] = True

            full_df = pd.DataFrame()
            if config.fdr_strategy == 'transferred_fdr' or config.fdr_strategy == 'transferred_fdr_linear_reg':
                full_df = pd.concat([ss_psms, ptm_df], ignore_index=True).sort_values('hyperscore')
                full_df['rank'] = range(1, len(full_df) + 1)
                target = full_df.query('decoy == False')
                decoy = full_df.query('decoy == True')
            else:
                ptm_df = ptm_df.sort_values('hyperscore')
                ptm_df['rank'] = range(1, len(ptm_df) + 1)
                target = ptm_df.query('decoy == False')
                decoy = ptm_df.query('decoy == True')

            try:
                threshold, q_values = calculate_threshold(decoy, target, config, ptm_name, ratio_info_dir)
            except:
                logger.error(f'The size of the analysis result after filtering by {ptm_name} : {0}')
                continue

            if config.fdr_strategy == 'transferred_fdr' or config.fdr_strategy == 'transferred_fdr_linear_reg':
                ptm_df = full_df.query(f'PTM == True & rank >= {threshold}')
            else:
                ptm_df = ptm_df.query(f'rank >= {threshold}')

            for rank_val, q_val in sorted(q_values.items()):
                ptm_df.loc[ptm_df['rank'] > rank_val, 'q_value'] = q_val
            logger.info(f'The size of the analysis result after filtering by {ptm_name} : {ptm_df.shape}')

            ptm_df['accessions_list'] = ptm_df['protein'].apply(lambda x: '|'.join([p.split('|')[1] for p in ast.literal_eval(x)]))

            if config.search_mode == 'fast_search':
                valid_accs = set(ptm_uniprot_info_df.query(f'PTM == "{ptm_name}"')['accession'].unique())
                ptm_df = clear_ptm_psms(ptm_df, valid_accs)
                logger.info(f'The size of the analysis result after filtering for information in UniProt by {ptm_name} : {ptm_df.shape}')

            all_fdr_ptm_psms = pd.concat([all_fdr_ptm_psms, ptm_df], ignore_index=True)

        all_fdr_ptm_psms['scan'] = all_fdr_ptm_psms['spectrum'].str.split('scan=').str[1]
        all_fdr_ptm_psms['file_scan_peptide_id'] = all_fdr_ptm_psms['file_name']+'.'+all_fdr_ptm_psms['scan']+'.'+all_fdr_ptm_psms['peptide']
        all_fdr_ptm_psms = all_fdr_ptm_psms.loc[all_fdr_ptm_psms.groupby(['file_scan_peptide_id'])['hyperscore'].idxmax()]

        all_fdr_ptm_psms.to_csv(fdr_result_file_path, index=False)

        logger.info(f'Result dataframe :\n{all_fdr_ptm_psms.shape}')
    else:
        all_fdr_ptm_psms = pd.read_csv(fdr_result_file_path)

    fdr_analysis_dir = config.ptm_search_dir / f'Result_of_PTM_analysis_{config.search_mode}_{config.fdr_strategy}'
    fdr_analysis_dir.mkdir(exist_ok=True)

    ss_peptides = pd.read_csv(config.st_search_dir / 'union_peptides.tsv', sep='\t')
    get_plots_from_result_of_analysis(all_fdr_ptm_psms, ss_peptides, config, fdr_analysis_dir)
