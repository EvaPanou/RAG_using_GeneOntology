# Gene Ontology retrieval-augmented generation
### *Retrieval of Gene Ontology terms with sentence embeddings, answers from a small local language model, and a gene list test against an enrichment analysis*

## Topic
This project asks whether retrieving the five closest Gene Ontology terms for a question<sup>1,2</sup> makes the answer of a language model checkable, and whether the retrieval agrees with a standard enrichment analysis on a human gene list. Retrieval uses sentence embeddings<sup>4,5</sup> and one matrix multiplication in PyTorch<sup>7</sup>; answers were written by a small local model, Qwen2.5-0.5B-Instruct<sup>8</sup>. The Claude backend is built into the notebook but was **not run**. It is a learning project, not a scientific or clinical tool.

## Repository structure and data sources
| Path | What it is | Source and how it was obtained |
|---|---|---|
| `README.md` | Description, setup, notes, results and references | |
| `00_requirements.txt` | Pinned package versions | Written by Step 1 of the notebook |
| `01_missmethyl_check.R` | Reference test of the age-related methylation sites with the Bioconductor package missMethyl | Written, **not yet run** when this README was written; no result of it is reported here |
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

5. Run `02_Gene_Ontology_analysis.ipynb` from top to bottom. Step 2 downloads the two data files into `data/`. The embedding step (38,092 terms) took 16 minutes and 17 seconds on the processor; no graphics card is needed. The first answer of the local model also downloads the model.
6. A box at the top of the Visual Studio Code window asks for a Claude key in Step 7. **Press Enter on the empty box to use only the local model** (what was done here). A key is never written to a file; it can also be given in the environment variable `ANTHROPIC_API_KEY`.
7. Optional, for the terminal: `bash 04_full_execution.sh` in Git Bash. This also needs R for the missMethyl step.

## Intro
The Gene Ontology is a controlled vocabulary of terms that describe what genes do. Each term has an identifier (written like `GO:0006914`), a name and a definition<sup>1,2</sup>. A language model that answers from memory can invent such identifiers. **Retrieval-augmented generation**<sup>3</sup> first searches the real ontology for the terms closest to the question and gives only those to the model, so every cited identifier can be checked.

Each term (name plus definition) is turned into 384 numbers with the sentence-embedding model all-MiniLM-L6-v2<sup>4,5</sup>. Cosine similarity between a question and all terms is one matrix multiplication in PyTorch<sup>7</sup>, and the five highest are kept. The model is not trained here: the retrieval step has 0 learnable numbers and 22,713,216 frozen numbers. The citation check is a pattern match for `GO:` followed by seven digits, compared with the retrieved list. For the gene list test, the answer key is an over-representation analysis: a hypergeometric test for each biological process term, with the Benjamini-Hochberg correction<sup>6</sup>, after annotations were propagated to the ancestors of each term.

## Notes about the data
- **Obsolete terms:** 48,340 terms were read and 10,248 obsolete terms removed, leaving 38,092 (23,974 biological process, 10,041 molecular function, 4,077 cellular component). No kept term lacks a definition, and no name or definition appears twice.
- **Annotations:** 910,476 lines were read; 1,497 with a negative relation (`NOT`) were skipped. 38,810 annotated entries form the background, and 21,102 terms have at least one entry after propagation.
- **Entries that are not genes:** the UniProt-based annotation file includes RNAcentral identifiers (names starting `URS`), protein complexes (names ending `_human`) and microRNAs (for example `hsa-miR-650`). They are part of the background and appear in the planted gene lists (for seed 0, 8 of 40 entries are `URS` identifiers). Their effect on the results was not tested here.
- **Cut-off texts:** 30 terms have a text longer than the 256 word pieces the embedding model reads; the end of these texts is cut off.
- **Two questions skipped:** the names `aging` and `DNA methylation` were not found in this ontology release, so 8 of the 10 test questions had a known correct term. Why was not tested here.
- **Download:** the first automatic download of the two data files was refused with an HTTP 403 error. It worked after the download cell was given a browser-like `User-Agent` header (three lines at the top of the Step 2 cell). Later runs of Step 2 report `already on disk`.
- **Planted gene list:** the gene lists are built from the genes of `autophagy` plus 10 random genes, using the same annotation file as the test (circular, see Limits).

## Results and conclusions
All numbers come from the saved notebook. Only the local model (Qwen2.5-0.5B-Instruct) was run, so nothing here says how a stronger model behaves.

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

- **Retrieval.** Embedding search found the correct term among the five in more questions than word search (5 against 2) and far more than chance. With only 8 questions, a difference of three is too small to call reliable.
- **Generation.** When the local model cited an identifier, it was real and among the retrieved five. But 7 of 8 answers with retrieval cited none, and without retrieval it cited none either. This test therefore **could not show whether retrieval prevents invented identifiers**.
- **Gene list.** The enrichment analysis recovered the planted term, but the list was built from the same annotation file, so this only shows that the code works. The terms retrieved for the gene list matched none of the answer-key terms, and the random control gave significant terms in two of three seeds. The reasons were not tested here.
- Two of the ten test questions (`aging` and `DNA methylation`) were skipped because their exact names were not found in this ontology release.

### Environment used
| Item | Version or detail |
|---|---|
| Operating system | Windows 11, version 10.0.26300 |
| Python | 3.12.0 |
| Visual Studio Code | `1.146` |
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
The sites chosen by the elastic net in my other repository [Predicting_age_from_DNA_methylation](https://github.com/EvaPanou/Predicting_age_from_DNA_methylation) are stored in `input/methylation_repo/`. The step that maps these sites to genes and tests their Gene Ontology terms has **not been built yet**, so no biological claim about those sites is made here.

## References
1. Ashburner M, Ball CA, Blake JA, et al. Gene Ontology: tool for the unification of biology. *Nat Genet*. 2000;25(1):25-29. doi:10.1038/75556
2. The Gene Ontology Consortium; Aleksander SA, Balhoff J, Carbon S, et al. The Gene Ontology knowledgebase in 2023. *Genetics*. 2023;224(1):iyad031. doi:10.1093/genetics/iyad031
3. Lewis P, Perez E, Piktus A, et al. Retrieval-augmented generation for knowledge-intensive NLP tasks. In: *Advances in Neural Information Processing Systems 33 (NeurIPS 2020)*. 2020. Page range not verified and left out.
4. Reimers N, Gurevych I. Sentence-BERT: sentence embeddings using Siamese BERT-networks. In: *Proceedings of the 2019 Conference on Empirical Methods in Natural Language Processing and the 9th International Joint Conference on Natural Language Processing (EMNLP-IJCNLP)*. Association for Computational Linguistics; 2019:3982-3992. doi:10.18653/v1/D19-1410
5. sentence-transformers/all-MiniLM-L6-v2. Hugging Face. Accessed November 18, 2026. https://huggingface.co/sentence-transformers/all-MiniLM-L6-v2
6. Benjamini Y, Hochberg Y. Controlling the false discovery rate: a practical and powerful approach to multiple testing. *J R Stat Soc Series B Methodol*. 1995;57(1):289-300. doi:10.1111/j.2517-6161.1995.tb02031.x
7. Paszke A, Gross S, Massa F, et al. PyTorch: an imperative style, high-performance deep learning library. In: Wallach H, Larochelle H, Beygelzimer A, d'Alché-Buc F, Fox E, Garnett R, eds. *Advances in Neural Information Processing Systems 32*. Curran Associates, Inc.; 2019:8024-8035.
8. Qwen, Yang A, Yang B, et al. Qwen2.5 technical report. arXiv. Preprint posted online December 19, 2024. arXiv:2412.15115. Model used: Qwen/Qwen2.5-0.5B-Instruct, Hugging Face. Accessed November 18, 2026. https://huggingface.co/Qwen/Qwen2.5-0.5B-Instruct