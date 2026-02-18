# ptms_closed_search

**ptms_closed_search** is a module for running a complete post-translational modification (PTM) analysis of proteins based on mass spectrometry results, using IdentiPy and info about human proteins from UniProt and dbPTM databases.

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
pip install git+https://github.com/Yura17234/ptms_closed_search.git
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
work_dir = work_dir
uniprot_query_path = Query_uniprotkb_proteome_UP000005640_AND_prot_xxxx_xx_xx.txt
fasta_path = uniprotkb_proteome_UP000005640_AND_revi_xxxx_xx_xx_reverse.fasta
base_config_path = work_dir/edited_default.cfg
additional_lists_path = work_dir/additional_lists.txt
# additional_lists.txt or / if no additional proteins for ptm search
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
---
# Citation

**ptms_closed_search** is described in this paper:
Strogov, Y. Y., Spirin, S. A., Ivanov, M. V., Kulebyakina, M. A., Efimenko, A. Y., & Klychnikov, O. I. (2026).
PTMs_Closed_Search: Multiple Post-Translational Modification Closed Search Using Reduced Search Space and Transferred FDR.
Proteomes, 14(1), 7. https://doi.org/10.3390/proteomes14010007

Please cite it when using ptms_closed_search or its parts.
