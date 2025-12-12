#!/bin/bash
# Do not exit immediately on error for the entire script
# set -e removed here

N=20

source "$HOME/miniforge3/etc/profile.d/conda.sh"
conda activate ligandmpnn_env

cd /home/tedeschg/prj/dockndesign/pkgs/LigandMPNN

POCKET_RESIDUES=$(cat /home/tedeschg/prj/dockndesign/test_score/test_multi_gnina/pocket_residues.txt)

for i in $(seq 1 $N); do

    IN_PDB="/home/tedeschg/prj/dockndesign/test_score/test_multi_gnina/complex${i}.pdb"
    OUT_DIR="/home/tedeschg/prj/dockndesign/test_score/test_multi_gnina/lmpnn${i}"

    if [ ! -f "$IN_PDB" ]; then
        echo "SKIP: $IN_PDB not found"
        continue
    fi

    echo "==============================="
    echo " Running LigandMPNN for complex${i}.pdb"
    echo " Output -> ${OUT_DIR}"
    echo "==============================="

    # Catch any errors from the python command
    if python run.py \
        --model_type "ligand_mpnn" \
        --checkpoint_ligand_mpnn "model_params/ligandmpnn_v_32_010_25.pt" \
        --pdb_path "$IN_PDB" \
        --out_folder "$OUT_DIR" \
        --temperature 0.1 \
        --seed 111 \
        --redesigned_residues "$POCKET_RESIDUES"; then
        echo "Successfully completed complex${i}"
    else
        echo "ERROR during processing of complex${i}"
    fi

done

echo "==============================="
echo " Processing completed"
echo "==============================="