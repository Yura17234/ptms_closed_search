from setuptools import setup, find_packages

setup(
    name="test_module",
    version="0.1.0",
    packages=find_packages(),
    install_requires=[],  # сюда добавь зависимости при необходимости
    entry_points={
        'console_scripts': [
            'run_prepare_ptm_search = test_module.run_prepare_ptm_search:main',
            'run_multiple_search = test_module.run_multiple_search:main',
            'run_aggregate_results = test_module.run_aggregate_results:main',
        ],
    },
)
