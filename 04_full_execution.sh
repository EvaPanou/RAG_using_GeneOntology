#!/usr/bin/env bash
# 04_full_execution.sh - every terminal step of the project, in order.
# Run from the repository folder.
#   Windows: open "Git Bash" (installed together with Git for Windows) and type:  bash 04_full_execution.sh
#   macOS / Linux: bash 04_full_execution.sh
# The R part is skipped when its results already exist in results/. To run it again anyway, type:
#   RERUN_R=1 bash 04_full_execution.sh
set -e    # stop at the first error

# ---------- Part 1: environment ----------
if [ ! -d ".venv" ]; then
    python -m venv .venv || python3 -m venv .venv
fi
if [ -f ".venv/Scripts/activate" ]; then
    source .venv/Scripts/activate      # Windows (Git Bash)
else
    source .venv/bin/activate          # macOS / Linux
fi
python -m pip install --upgrade pip
python -m pip install -r 00_requirements.txt
python -m pip install nbconvert ipykernel      # needed only to run the notebook from the terminal

# ---------- Part 2: R reference test with missMethyl ----------
# This part takes long (the first run also installs Bioconductor packages). It is skipped when its result already exists.
if [ -f "results/missmethyl_all_seeds.csv" ] && [ "$RERUN_R" != "1" ]; then
    echo "The missMethyl results already exist in results/. Skipping the R part."
else
    for required_file in \
        "input/methylation_repo/bridge_information.json" \
        "input/methylation_repo/elastic_net_sites_seed_0.csv" \
        "input/methylation_repo/elastic_net_sites_seed_1.csv" \
        "input/methylation_repo/elastic_net_sites_seed_42.csv" \
        "input/methylation_repo/all_sites_in_study.txt"; do
        if [ ! -f "$required_file" ]; then
            echo "Missing file: $required_file"
            exit 1
        fi
    done
    echo "All methylation input files found."
    if ! command -v Rscript > /dev/null 2>&1; then
        echo "Rscript was not found. Install R from https://cloud.r-project.org and open a new Git Bash window, then run this file again."
        exit 1
    fi
    Rscript 01_missmethyl_check.R
fi

# ---------- Part 3: run the notebook from top to bottom ----------
# Optional: to use the Claude backend without typing the key, set it before running (never write it in a file):
#   export ANTHROPIC_API_KEY="..."
python -m jupyter nbconvert --to notebook --execute --inplace --ExecutePreprocessor.timeout=-1 02_Gene_Ontology_analysis.ipynb
echo "Notebook finished. Open it in Visual Studio Code to read the outputs."