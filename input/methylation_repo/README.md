# input

Files made by the separate methylation project & can be located there(`Predicting_age_from_DNA_methylation`, study GSE179325, platform GPL21145, Illumina MethylationEPIC). They are the starting point for linking age-related methylation sites to Gene Ontology terms. **The notebook of this repository does not use them yet.**

Folder `methylation_repo/`:

| File | What it is |
|---|---|
| `elastic_net_sites_seed_0.csv`, `_1.csv`, `_42.csv` | The 2,000 sites offered to the elastic net for each seed, with columns `site`, `weight` (years per standard deviation of the site value) and `weight_is_not_zero`. Weights other than zero: 126 (seed 0), 755 (seed 1), 186 (seed 42). The fit of seed 1 did not converge, so its list is less reliable |
| `bridge_information.json` | Study, platform, seeds, number of sites (768,067) and date |
| `all_sites_in_study.txt` | The names of all 768,067 sites, one per line. Not committed (about 8 megabytes); it is made again by Step 10.4 of the methylation notebook |

Sites named `ch.…` are not CpG sites and have no gene annotation; they are removed before any gene lookup.