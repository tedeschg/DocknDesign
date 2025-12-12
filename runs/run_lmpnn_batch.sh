#!/bin/bash
# Exit on error
set -e

# Numero di complessi
N=20

# Base directory: current directory where you run the script
BASE_DIR="$(pwd)"

# Attiva environment conda ligandmpnn_env
source "$HOME/miniforge3/etc/profile.d/conda.sh"
conda activate ligandmpnn_env

# Vai nella cartella LigandMPNN
cd /home/tedeschg/prj/dockndesign/pkgs/LigandMPNN

# Leggi residui pocket (from current directory)
POCKET_RESIDUES="$(cat "${BASE_DIR}/pocket_residues.txt")"

OUTDIR="${BASE_DIR}/lmpnn_results"
mkdir -p "$OUTDIR"

# Loop su complex1.pdb ... complex20.pdb
for i in $(seq 1 $N); do

    IN_PDB="${BASE_DIR}/complexes/complex${i}.pdb"
    OUT_DIR="${OUTDIR}/lmpnn${i}"

    mkdir -p "$OUT_DIR"

    echo "==============================="
    echo " Running LigandMPNN for complex${i}.pdb"
    echo " Input  → ${IN_PDB}"
    echo " Output → ${OUT_DIR}"
    echo "==============================="

    python run.py \
        --model_type "ligand_mpnn" \
        --checkpoint_ligand_mpnn "model_params/ligandmpnn_v_32_010_25.pt" \
        --pdb_path "$IN_PDB" \
        --out_folder "$OUT_DIR" \
        --temperature 0.1 \
        --seed 37 \
        --redesigned_residues "$POCKET_RESIDUES"

done
