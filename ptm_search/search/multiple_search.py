import subprocess
import pandas as pd
from pyteomics import pepxml

def multiple_search(config):
    '''
        Multiple search
    '''

    configs_dir = config.ptm_search_dir / f"{config.analysis_index}_Configs_{config.search_mode}"
    results_dir = config.ptm_search_dir / f"{config.analysis_index}_result_{config.search_mode}"
    results_dir.mkdir(exist_ok=True)

    config_files = sorted(configs_dir.glob("*.cfg"))
    print(f'\nНайдено конфигураций:\n{len(config_files)}\n')

    already_done = {d.name for d in results_dir.iterdir() if d.is_dir()}

    for count, cfg_path in enumerate(config_files, 1):
        mod_name = cfg_path.stem
        if mod_name in already_done:
            print(f"Pass {mod_name} ({count}/{len(config_files)})")
            continue

        print(f"{'':-^{50}}\n{count}/{len(config_files)} | {mod_name}\n{'':-^{50}}")

        # === Запуск IdentiPy ===
        # mgf_pattern = str(config.work_dir / "*.mgf")
        mgf_files = list(config.work_dir.glob('*.mgf'))
        identipy_cmd = [
            "identipy",
            mgf_files,
            "-cfg", str(cfg_path)
        ]

        try:
            subprocess.run(identipy_cmd, check=True)
        except subprocess.CalledProcessError as e:
            print(f"\nОшибка при запуске IdentiPy:\n{e}")
            continue

        # === Сбор результатов .pep.xml ===
        ptm_all_df = pd.DataFrame()
        for pepxml_file in config.work_dir.glob("*.pep.xml"):
            file_name = pepxml_file.stem
            df_pep = pepxml.DataFrame(str(pepxml_file))
            df_pep['file_name'] = file_name
            ptm_all_df = pd.concat([ptm_all_df, df_pep], ignore_index=True)

        # === Сохраняем результаты ===
        mod_result_dir = results_dir / mod_name
        mod_result_dir.mkdir(exist_ok=True)

        result_csv_path = results_dir / mod_name / "result.csv"
        ptm_all_df.to_csv(result_csv_path, index=False)

        # Удаление .pep.xml
        for f in config.work_dir.glob("*.pep.xml"):
            f.unlink()

    print('multiple_search -- Done !')
