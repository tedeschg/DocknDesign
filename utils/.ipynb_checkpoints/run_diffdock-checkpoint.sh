#!/bin/bash
# Exit on error
set -e

# Base directory (directory da cui lanci lo script)
BASE_DIR="$(pwd)"

# Attiva l’environment conda diffdock
source "$HOME/miniforge3/etc/profile.d/conda.sh"
conda activate diffdock

# Output directory
OUTDIR="${BASE_DIR}/diffdock_results/test_multi"
mkdir -p "$OUTDIR"

# Vai nella cartella DiffDock
cd /home/tedeschg/prj/dockndesign/pkgs/DiffDock

# Esegui inference
python -m inference \
  --config default_inference_args.yaml \
  --protein_ligand_csv "${BASE_DIR}/pdk1_multi.csv" \
  --out_dir "$OUTDIR"
