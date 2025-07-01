# PTM Pipeline

**PTM Pipeline** is a module for running a complete post-translational modification (PTM) analysis of proteins based on mass spectrometry results, using IdentiPy and UniProt / dbPTM databases.

---

## Features

- Parsing standard protein search results.
- Matching identified proteins to PTM data from UniProt and dbPTM.
- Generating FASTA and configuration files for PTM search.
- Running batch PTM searches with IdentiPy using various configurations.
- Aggregating and filtering PTM search results.
- Visualizing results via plots and summary tables.

---
## Installation

```bash
pip install git+https://github.com/Yura17234/ptm_searching.git
```

---
## PTM search workflow
The PTM search is performed in three steps using the following commands. Each command reads all required parameters and file paths from a single configuration file like this one:

```ini
[experiment]
experiment_name = HEK293
analysis_index = all
# can be 01, 02, 03 or all
search_mode = fast_search
# can be fast_search or full_search
fdr_strategy = transferred_fdr
# can be transferred_fdr or separate_fdr

[paths]
work_dir = identipy_work/HEK293
uniprot_query_path = Query_uniprotkb_proteome_UP000005640_AND_prot_2024_07_27.txt
fasta_path = UP000005640_9606_reverse.fasta
base_config_path = identipy_work/edited_default.cfg
```

```bash
run_prepare_ptm_search --config parameters.cfg
```

```bash
run_multiple_search --config parameters.cfg
```

```bash
run_aggregate_results --config parameters.cfg
```