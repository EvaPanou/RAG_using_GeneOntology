# results

Files written by the notebook `02_Gene_Ontology_analysis.ipynb`. Only the local model (Qwen2.5-0.5B-Instruct) was run.

| File | Step | What it contains |
|---|---|---|
| `summary.txt` | 11 | Every result line printed by the notebook (66 lines) |
| `retrieval_per_question.csv` | 8.1 | For each question with a known correct term: whether the correct or a related term was among the five found by embedding search and by word search |
| `comparison_per_question.csv` | 8.2 | For each question, backend and condition (with or without retrieval): identifiers cited, how many exist, how many are inside the retrieved five, and whether the strict citation check passed |
| `answers_with_and_without_retrieval.txt` | 8.2 | The full text of every answer, with its citation check |
| `gene_list_per_seed.csv` | 10 | For seeds 0, 1 and 42: significant terms, rank of the planted term, significant terms in a random list, and overlap of the retrieval with the answer key |

`01_missmethyl_check.R` writes `missmethyl_*.csv` files here once it has been run.