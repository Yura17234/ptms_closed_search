from pathlib import Path

def prepare_ptm_search(config):
    '''
        Prepare PTM search
    '''

    current_dir = Path.cwd()
    home_dir = Path.home()
    module_dir = Path(__file__).parent.resolve()

    ptm_search_path = home_dir / config.work_dir / f"{config.experiment_name}_PTM_search"
    st_search_path = home_dir / config.work_dir / f"{config.experiment_name}_Standard_search"
    mgf_path = st_search_path / "full_mgf_files"

    # os.makedirs(ptm_search_path, exist_ok=True)
    # so.makedirs(st_search_path, exist_ok=True)
    # os.makedirs(mgf_path, exist_ok=True)
    #
    # for ext in ['*.mgf', '*.tsv', '*.pep.xml', '*.png']:
    #     for file in glob.glob(os.path.join(home_dir, config.work_dir, ext)):
    #         dest = mgf_path if ext == '*.mgf' else st_search_path
    #         shutil.move(file, dest)
    #
    # st_search_df = pd.read_csv(os.path.join(st_search_path, f"union_protein.tsv"), sep='\t')

    ''' 1 '''
    # list_of_grouped_prots_by_ptms, dict_acc_to_names = parsing_human_proteom(ptm_search_path, config, st_search_df)

    # print(f'Количество белков из стандартного начального поиска: {len(list(dict_acc_to_names.keys()))}')
    # print(st_search_df.head())
    # list_of_grouped_prots_by_ptms = adding_ptm_info_from_db_ptm(list_of_grouped_prots_by_ptms,
    #                                                             list(dict_acc_to_names.keys()))
    #
    # ''' 2 '''
    # ptm_info_df = make_ptms_df(list_of_grouped_prots_by_ptms, dict_acc_to_names, ptm_search_path, config)
    # print(ptm_info_df.head())
    # print()
    #
    # ptm_counts = ptm_info_df.drop_duplicates(subset=['accession', 'PTM'], keep='first')['PTM'].value_counts(
    #     sort=True).head(15)
    # ptm_counts.plot(kind='barh', title="Top PTMs (UniProt)", figsize=(5, 8), color=(0.1, 0.7, 0.1, 0.7),
    #                 edgecolor='black')
    # plt.xlabel("Protein count")
    # plt.savefig(os.path.join(ptm_search_path, "PTMs_count_of_potential_prots_ENG.png"), dpi=200, bbox_inches='tight')
    # plt.close()
    #
    # ''' 3 '''
    # with open(fasta_path, 'r') as fasta_file:
    #     make_fasta_file_for_searches(list(ptm_info_df['accession'].unique()), list_of_grouped_prots_by_ptms, fasta_file, ptm_search_path, config)
    #
    # ''' 4 '''
    # ptm_list = list(ptm_info_df['PTM'].unique())
    # for mode in [1, 2]:  # малые / большие базы
    #     parse_config_file(ptm_list, ptm_search_path, config, mode)
    #
    # ''' 5 '''
    # # Создание mgf-файлов из не идентифицированных спектров для PTM-поиска
    # if all(['_for_PTM.mgf' not in file for file in os.listdir(f'{base_path}/')]):
    #     make_mgfs_for_ptm(base_path, experiment_name)

    print(ptm_search_path)
    print(st_search_path)
    print(mgf_path)
    print('prepare_ptm_search -- connected')
