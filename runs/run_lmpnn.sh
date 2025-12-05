#!/bin/bash
# Exit on error
set -e

# Attiva environment conda ligandmpnn_env
source "$HOME/miniforge3/etc/profile.d/conda.sh"
conda activate ligandmpnn_env

# Vai nella cartella LigandMPNN
cd /home/tedeschg/prj/dockndesign/pkgs/LigandMPNN

# Leggi i residui dal file pocket_residues.txt
POCKET_RESIDUES=$(cat /home/tedeschg/prj/dockndesign/pdbs/hiv/pocket_residues.txt)

# Esegui LigandMPNN
python run.py \
    --model_type "ligand_mpnn" \
    --checkpoint_ligand_mpnn "model_params/ligandmpnn_v_32_010_25.pt" \
    --pdb_path "/home/tedeschg/prj/dockndesign/hiv/complex.pdb" \
    --out_folder "/home/tedeschg/prj/dockndesign/hiv/lmpnn" \
    --temperature 0.1 \
    --seed 37 \
    --redesigned_residues "$POCKET_RESIDUES"
