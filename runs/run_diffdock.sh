#!/bin/bash
# Exit on error
set -e

# Attiva l’environment conda diffdock
# (modifica il path a conda.sh se necessario)
source "$HOME/miniforge3/etc/profile.d/conda.sh"
conda activate diffdock

# Vai nella cartella DiffDock
cd /home/tedeschg/prj/dockndesign/pkgs/DiffDock

# Esegui inference
python -m inference \
  --config default_inference_args.yaml \
  --protein_ligand_csv /home/tedeschg/prj/dockndesign/pdbs/hiv/4LL3.csv \
  --out_dir /home/tedeschg/prj/dockndesign/hiv
