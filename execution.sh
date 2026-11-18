#!/usr/bin/env bash
# run_all.sh - every terminal step of the project, in order.
# Run from the repository folder.
#   Windows: open "Git Bash" (installed together with Git for Windows) and type:  bash run_all.sh
#   macOS / Linux: bash run_all.sh
# Each part is marked. Parts that are not needed yet are commented out and will be switched on as the project grows.
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
python -m pip install -r requirements.txt
python -m pip install nbconvert ipykernel      # needed only to run the notebook from the terminal

# ---------- Part 2: check that the methylation input files are present ----------
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

# ---------- Part 3: run the notebook from top to bottom ----------
# Optional: to use the Claude backend without typing the key, set it before running (never write it in a file):
#   export ANTHROPIC_API_KEY="..."
python -m jupyter nbconvert --to notebook --execute --inplace --ExecutePreprocessor.timeout=-1 Gene_Ontology_analysis.ipynb
echo "Notebook finished. Open it in Visual Studio Code to read the outputs."

# ---------- Part 4 (later): optional R cross-check ----------
# Rscript r/missmethyl_check.R
