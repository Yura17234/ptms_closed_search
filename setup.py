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
        "pandas>=1.0",
    ],
    entry_points={
        'console_scripts': [
            'run_prepare_ptm_search = ptm_search.run_prepare_ptm_search:main',
            'run_multiple_search = ptm_search.run_multiple_search:main',
            'run_aggregate_results = ptm_search.run_aggregate_results:main',
        ],
    },
)
