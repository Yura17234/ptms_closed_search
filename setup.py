from setuptools import setup, find_packages

setup(
    name="ptm_search",
    version="0.1.0",
    packages=find_packages(),
    include_package_data=True,
    package_data={
        "ptm_search": [
            "data/*.json",
            "data/*.json.gz",
        ],
    },
    install_requires=[
        "pandas>=1.3",
        "numpy>=1.21",
        "scipy>=1.7",
        "matplotlib>=3.4",
        "seaborn>=0.11",
        "tqdm>=4.62",
        "lxml>=4.6",
        "pyteomics>=4.5",
        "ConfigUpdater>=3.1",
        "configparser>=5.0",
        "unipressed>=0.3.1",
        "identipy @ git+https://github.com/levitsky/identipy.git#egg=identipy",
        "setuptools>=65.0.0",
        "scikit-learn>=1.0",
        "openpyxl>=3.0",
    ],
    entry_points={
        'console_scripts': [
            'run_prepare_ptm_search = ptm_search.run_prepare_ptm_search:main',
            'run_multiple_search = ptm_search.run_multiple_search:main',
            'run_aggregate_results = ptm_search.run_aggregate_results:main',
        ],
    },
)
