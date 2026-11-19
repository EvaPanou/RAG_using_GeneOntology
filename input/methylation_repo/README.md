# Input folder: sites from the methylation project

These files come from a separate repository, `Predicting_age_from_DNA_methylation`, which predicts age from blood DNA methylation with an elastic net and a small multilayer perceptron (study GSE179325, platform GPL21145, Illumina Infinium MethylationEPIC version 1). They were created by Step 10.4 of the notebook of that repository on 2026-11-18 and copied here without changes. This repository only reads them (Step 11 of `02_Gene_Ontology_analysis.ipynb` and `01_missmethyl_check.R`).

## Files

| File | What it holds | In Git? |
|---|---|---|
| `elastic_net_sites_seed_0.csv`, `elastic_net_sites_seed_1.csv`, `elastic_net_sites_seed_42.csv` | One file per random seed. Columns: `site` (the array probe name), `weight` (the weight the elastic net gave the site in that seed; `0.0` or `-0.0` means the site was not used) and `weight_is_not_zero` (`True` or `False`). The seed 0 file has 2,000 rows, 126 of them with a weight other than zero | yes |
| `bridge_information.json` | Study `GSE179325`, platform `GPL21145`, seeds `0, 1, 42`, `sites_in_study` 768,067, date 2026-11-18 | yes |
| `all_sites_in_study.txt` | One site name per line, all 768,067 sites of the study (765,569 start with `cg`, 2,498 are `ch.` sites). Used as the background of the tests | **no** (listed in `.gitignore`; copy it from the methylation repository, or recreate it with Step 10.4 of its notebook) |

## Notes

- Site names such as `cg16867657` are probe names of the array, not genes. The table that links them to genes is `data/site_to_gene_table.csv`.