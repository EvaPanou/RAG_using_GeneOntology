# Results folder

Files written by the notebook (`02_Gene_Ontology_analysis.ipynb`) and by `01_missmethyl_check.R`. They are the evidence for the numbers in the main README. In the methylation files, `<list>` is `seed_0`, `seed_1`, `seed_42` or `all_seeds`.

## Retrieval, generation and gene list test (Steps 7 to 10)

| File | What it holds |
|---|---|
| `summary.txt` | Every reported line of the notebook, in order, including the Step 11 lines |
| `retrieval_per_question.csv` | The terms retrieved for each test question |
| `comparison_per_question.csv` | Embedding search against word search, per question |
| `answers_with_and_without_retrieval.txt` | The answers of the local model with and without retrieval |
| `gene_list_per_seed.csv` | The gene list test for seeds 0, 1 and 42 (significant terms, rank of the seed term, random control, retrieval overlap, identifiers cited) |

## Methylation lists (Step 11 and the R script)

| File | What it holds | Written by |
|---|---|---|
| `missmethyl_<list>.csv` | The `missMethyl` result (`gometh`, array type EPIC, all study `cg` sites as background, correction for the number of sites per gene). First column: GO identifier. Then `ONTOLOGY`, `TERM`, `N` (genes in the term), `DE` (genes of the list in the term), `P.DE` and `FDR`. For the top 20 terms two columns added by the script: `genes_in_list` (gene names) and `genes_listed_count` | `01_missmethyl_check.R` |
| `python_enrichment_<list>.csv` | The over-representation test of the notebook inside the restricted background: `identifier`, `name`, `genes_in_list`, `genes_in_term`, `probability_value`, `adjusted_probability_value`. For the top 20 terms: `genes_in_list_names` and `permutation_probability_value` (share of 500 random site lists with an equal or smaller value for that term) | Step 11 |
| `bridge_summary.csv` | One row per list: sites, genes with a matching symbol, significant terms in this notebook and in `missMethyl`, smallest probability value, share of 500 random lists with an equal or smaller smallest value, number of the ten lowest terms shared with `missMethyl`, rank correlation (Spearman), number of terms compared | Step 11 |

## How to read them

- `genes_listed_count` in a `missmethyl_` file can be smaller than `DE`: the script looks the gene names up itself with another annotation table. The reason was not tested.
- The two analyses use different gene-to-term tables, so their results are not expected to agree exactly.
- In `python_enrichment_<list>.csv`, a probability value is only meaningful together with the permutation values and the share in `bridge_summary.csv`: the simple test does not correct for the number of sites per gene, and the permutation test shows that its small values are common in random lists. No list shows reliable enrichment.
- The smallest value that a permutation can give is 1 / 501 = 0.002.