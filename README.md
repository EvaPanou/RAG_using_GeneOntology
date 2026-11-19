# Gene Ontology retrieval-augmented generation
### *Retrieval of Gene Ontology terms with sentence embeddings, answers from a small local language model, and a gene list test against an enrichment analysis*

## Topic
This project asks whether retrieving the five closest Gene Ontology terms for a question<sup>1,2</sup> makes the answer of a language model checkable, and whether the retrieval agrees with a standard enrichment analysis on a human gene list. A last step uses real lists of methylation sites from a separate project and tests them for Gene Ontology enrichment with a permutation test and with the R package missMethyl<sup>9,10</sup>. Retrieval uses sentence embeddings<sup>4,5</sup> and one matrix multiplication in PyTorch<sup>7</sup>; answers were written by a small local model, Qwen2.5-0.5B-Instruct<sup>8</sup>. The Claude backend is built into the notebook but was **not run**. It is a learning project, not a scientific or clinical tool.

## Repository structure and data sources
| Path | What it is | Source and how it was obtained |
|---|---|---|
| `README.md` | Description, setup, notes, results and references | |
| `00_requirements.txt` | Pinned package versions | Written by Step 1 of the notebook |
| `01_missmethyl_check.R` | Reference test of the age-related methylation sites with the Bioconductor package missMethyl | |
| `02_Gene_Ontology_analysis.ipynb` | The analysis notebook, saved with its outputs | |
| `04_full_execution.sh` | Terminal commands for the whole project in order (Git Bash) | |
| `input/` | The sites chosen by the Elastic Net in my methylation project repo `Predicting_age_from_DNA_methylation`. |
| `data/` | The ontology, the human annotation file and the cache files; only `data/README.md` is committed. See `data/README.md` | Downloaded and created by the notebook (Steps 2, 4 and 7) |
| `results/` | The result tables and text files of the run. See `results/README.md` | Created by the notebook |

## Setup
1. Create a virtual environment in the project folder, with Python 3.12:

   ```text
   python -m venv .venv
   ```

2. Activate it:

   ```text
   # Windows (PowerShell) - the one I used
   .venv\Scripts\Activate.ps1

   # macOS and Linux
   source .venv/bin/activate
   ```

3. Install the packages:

   ```text
   pip install -r 00_requirements.txt
   ```

4. Open the notebook and select the `.venv` environment as its kernel:
   - **Visual Studio Code *(used for this project)*:** click *Select Kernel*, and choose `.venv`.

5. Run `02_Gene_Ontology_analysis.ipynb` from top to bottom. Step 2 downloads the two data files into `data/`. No graphics card is needed (the embeddings of the 38,092 terms are computed on the processor and cached in `data/`). The first answer of the local model also downloads the model.

6. A box at the top of the Visual Studio Code window asks for a Claude key in Step 7. **Press Enter on the empty box to use only the local model** (what was done here). A key is never written to a file; it can also be given in the environment variable `ANTHROPIC_API_KEY`.

7. Optional, for the terminal: `bash 04_full_execution.sh` in Git Bash. It creates the environment, runs `01_missmethyl_check.R`, and executes the notebook in place.
   - R (4.6.1 was used) must be installed and `Rscript` must be on the path of Git Bash (for example `export PATH="$PATH:/c/Program Files/R/R-4.6.1/bin"` in `~/.bashrc`). The R script installs missing Bioconductor packages by itself on its first run.
   - The R step is skipped when `results/missmethyl_all_seeds.csv` and `data/site_to_gene_table.csv` already exist; `RERUN_R=1 bash 04_full_execution.sh` forces it. Inside the R script, lists whose result file exists are skipped as well.
   - The script installs `nbconvert` (it runs the notebook from the terminal); it is not part of `00_requirements.txt`.
   - In the terminal, no box can ask for the Claude key: the Claude backend is used only if `ANTHROPIC_API_KEY` is set, otherwise the notebook prints a message and uses the local model only.

## Intro
The Gene Ontology is a controlled vocabulary of terms that describe what genes do. Each term has an identifier (written like `GO:0006914`), a name and a definition<sup>1,2</sup>. A language model that answers from memory can invent such identifiers. **Retrieval-augmented generation**<sup>3</sup> first searches the real ontology for the terms closest to the question and gives only those to the model, so every cited identifier can be checked.

Each term (name plus definition) is turned into 384 numbers with the sentence-embedding model all-MiniLM-L6-v2<sup>4,5</sup>. Cosine similarity between a question and all terms is one matrix multiplication in PyTorch<sup>7</sup>, and the five highest are kept. The model is not trained here: the retrieval step has 0 learnable numbers and 22,713,216 frozen numbers. The citation check is a pattern match for `GO:` followed by seven digits, compared with the retrieved list. For the gene list test, the answer key is an over-representation analysis: a hypergeometric test for each biological process term, with the Benjamini-Hochberg correction<sup>6</sup>, after annotations were propagated to the ancestors of each term. For the methylation lists (Step 11) the same test is used inside a restricted background, and a permutation test checks it against the known bias of gene-set tests on methylation arrays<sup>10</sup>.

## Notes about the data
- **Obsolete terms:** 48,340 terms were read and 10,248 obsolete terms removed, leaving 38,092 (23,974 biological process, 10,041 molecular function, 4,077 cellular component). No kept term lacks a definition, and no name or definition appears twice.
- **Annotations:** 910,476 lines were read; 1,497 with a negative relation (`NOT`) were skipped. 38,810 annotated entries form the background, and 21,102 terms have at least one entry after propagation.
- **Entries that are not genes:** the UniProt-based annotation file includes RNAcentral identifiers (names starting `URS`), protein complexes (names ending `_human`) and microRNAs (for example `hsa-miR-650`). They are part of the background and appear in the planted gene lists (for seed 0, 8 of 40 entries are `URS` identifiers). Their effect on the results was not tested here.
- **Cut-off texts:** 30 terms have a text longer than the 256 word pieces the embedding model reads; the end of these texts is cut off.
- **Two questions skipped:** the names `aging` and `DNA methylation` were not found in this ontology release, so 8 of the 10 test questions had a known correct term. Why was not tested here.
- **Download:** the first automatic download of the two data files was refused with an HTTP 403 error. It worked after the download cell was given a browser-like `User-Agent` header (three lines at the top of the Step 2 cell). Later runs of Step 2 report `already on disk`.
- **Planted gene list:** the gene lists are built from the genes of `autophagy` plus 10 random genes, using the same annotation file as the test (circular, see Limits).

## Results and conclusions
All numbers come from the saved notebook. Only the local model (Qwen2.5-0.5B-Instruct) was run, so nothing here says how a stronger model behaves. The retrieval pipeline and the language models were not run on the methylation lists.

| Test | Result |
|---|---|
| Correct term among the five retrieved, embedding search | 5 of 8 questions (first place in 3) |
| Correct term among the five retrieved, naive word search (baseline) | 2 of 8 questions |
| Chance level for one question | 0.00013 |
| Local model with retrieval: answers that cited identifiers | 1 of 8 (5 identifiers, all real and among the retrieved five) |
| Local model without retrieval: identifiers cited | 0 in 8 answers |
| Gene list test, seeds 0, 1 and 42: rank of the planted term | 1, 2 and 1 |
| Gene list test: overlap of the retrieved terms with the answer key | 0 in every seed |
| Random list of 40 genes: significant terms | 4, 1 and 0 (expected: none) |
| Methylation lists (seeds 0, 1, 42 with 126, 753 and 186 sites, and the 36 sites shared by all seeds): significant terms in the test of this notebook | 7, 25, 13 and 3 |
| Same lists: significant terms in `missMethyl` (corrects for the number of sites per gene) | 0 in every list |
| Same lists: share of 500 random site lists whose best term is as extreme as the observed one | 0.232, 0.160, 0.026 and 0.399 |

- **Retrieval.** Embedding search found the correct term among the five in more questions than word search (5 against 2) and far more than chance. With only 8 questions, a difference of three is too small to call reliable.
- **Generation.** When the local model cited an identifier, it was real and among the retrieved five. But 7 of 8 answers with retrieval cited none, and without retrieval it cited none either. This test therefore **could not show whether retrieval prevents invented identifiers**.
- **Gene list.** The enrichment analysis recovered the planted term, but the list was built from the same annotation file, so this only shows that the code works. The terms retrieved for the gene list matched none of the answer-key terms, and the random control gave significant terms in two of three seeds. The reasons were not tested here.
- **Methylation lists.** No list shows reliable enrichment. The simple test of this notebook finds significant terms, but `missMethyl` finds none, and random site lists of the same size give a best term as extreme as the observed one in 16 to 40% of draws for three of the four lists. Seed 42 (2.6%) is below 0.05 but not below 0.0125 (0.05 divided by four lists), and the lists share sites, so at most it is a hint. The small probability values come mainly from gene clusters: 14 of 21 and 16 of 17 genes behind `homophilic cell-cell adhesion` in seeds 1 and 42 are protocadherin genes, and a site inside a cluster is named for many genes at once (not tested directly). The two analyses rank terms alike (rank correlation 0.76 to 0.83) but share between 0 and 8 of their ten lowest terms.
- Two of the ten test questions (`aging` and `DNA methylation`) were skipped because their exact names were not found in this ontology release.

### Environment used
| Item | Version or detail |
|---|---|
| Operating system | Windows 11, version 10.0.26300 |
| Python | 3.12.0 |
| Visual Studio Code | `1.147` |
| R | R 4.6.1 {Bioconductor 3.23, missMethyl 1.46.0} |
| Computing | Processor only, no graphics card |
| Date of the run | `18/11/2026` |

### Package versions
The exact versions are pinned in `00_requirements.txt` and were printed by Step 1 of the notebook.

| Package | Version |
|---|---|
| torch | 2.15.0 |
| numpy | 2.5.5 |
| pandas | 3.1.1 |
| scipy | 1.18.1 |
| sentence-transformers | 6.3.0 |
| transformers | 5.23.0 |
| huggingface_hub | 2.4.0 |
| anthropic | 2.10.0 |

## Link to the methylation project
The sites chosen by the elastic net in my other repository [Predicting_age_from_DNA_methylation](https://github.com/EvaPanou/Predicting_age_from_DNA_methylation) are stored in `input/methylation_repo/`. TStep 11 of the notebook and `01_missmethyl_check.R` read them, map them to genes and test them for Gene Ontology terms (see Results). The result is that **none of the four lists shows reliable enrichment**. This is not a claim that the sites are unrelated to biology: a list chosen to predict age is not a list of sites known to change with age. Next steps (see the notebook): count each gene cluster once or weight genes by their number of sites, map genes through Entrez identifiers instead of symbols, and run the retrieval pipeline on the lists only if a list gives significant terms to serve as an answer key.

## References
1. Ashburner M, Ball CA, Blake JA, et al. Gene Ontology: tool for the unification of biology. *Nat Genet*. 2000;25(1):25-29. doi:10.1038/75556
2. The Gene Ontology Consortium; Aleksander SA, Balhoff J, Carbon S, et al. The Gene Ontology knowledgebase in 2023. *Genetics*. 2023;224(1):iyad031. doi:10.1093/genetics/iyad031
3. Lewis P, Perez E, Piktus A, et al. Retrieval-augmented generation for knowledge-intensive NLP tasks. In: *Advances in Neural Information Processing Systems 33 (NeurIPS 2020)*. 2020. Page range not verified and left out.
4. Reimers N, Gurevych I. Sentence-BERT: sentence embeddings using Siamese BERT-networks. In: *Proceedings of the 2019 Conference on Empirical Methods in Natural Language Processing and the 9th International Joint Conference on Natural Language Processing (EMNLP-IJCNLP)*. Association for Computational Linguistics; 2019:3982-3992. doi:10.18653/v1/D19-1410
5. sentence-transformers/all-MiniLM-L6-v2. Hugging Face. Accessed November 18, 2026. https://huggingface.co/sentence-transformers/all-MiniLM-L6-v2
6. Benjamini Y, Hochberg Y. Controlling the false discovery rate: a practical and powerful approach to multiple testing. *J R Stat Soc Series B Methodol*. 1995;57(1):289-300. doi:10.1111/j.2517-6161.1995.tb02031.x
7. Paszke A, Gross S, Massa F, et al. PyTorch: an imperative style, high-performance deep learning library. In: Wallach H, Larochelle H, Beygelzimer A, d'Alché-Buc F, Fox E, Garnett R, eds. *Advances in Neural Information Processing Systems 32*. Curran Associates, Inc.; 2019:8024-8035.
8. Qwen, Yang A, Yang B, et al. Qwen2.5 technical report. arXiv. Preprint posted online December 19, 2024. arXiv:2412.15115. Model used: Qwen/Qwen2.5-0.5B-Instruct, Hugging Face. Accessed November 18, 2026. https://huggingface.co/Qwen/Qwen2.5-0.5B-Instruct
9. Phipson B, Maksimovic J, Oshlack A. missMethyl: an R package for analyzing data from Illumina's HumanMethylation450 platform. *Bioinformatics*. 2016;32(2):286-288. doi:10.1093/bioinformatics/btv560
10. Geeleher P, Hartnett L, Egan LJ, Golden A, Raja Ali RA, Seoighe C. Gene-set analysis is severely biased when applied to genome-wide methylation data. *Bioinformatics*. 2013;29(15):1851-1857. doi:10.1093/bioinformatics/btt311