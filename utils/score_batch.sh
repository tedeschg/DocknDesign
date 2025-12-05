#!/bin/bash
set -e

N=20

SCORE_SCRIPT="/home/tedeschg/prj/dockndesign/utils/lmpnn_score.py"
PDB="/home/tedeschg/prj/dockndesign/pdbs/pdk1/pdk1.pdb"
POCKET="pocket_residues.txt"

OUTFILE="scores_summary.txt"
echo "Complex    Score" > "$OUTFILE"
echo "-----------------" >> "$OUTFILE"

for i in $(seq 1 $N); do

    FOLDER="lmpnn${i}"
    FASTA="${FOLDER}/seqs/complex${i}.fa"

    echo "======================================"
    echo " Scoring LMPNN output: ${FASTA}"
    echo "======================================"

    # Esegui comando e cattura l'output
    OUTPUT=$(python "$SCORE_SCRIPT" \
        -f "$FASTA" \
        -p "$PDB" \
        -r "$POCKET")

    # Estrai la riga "FINAL SCORE : xxxx"
    SCORE=$(echo "$OUTPUT" | grep "FINAL SCORE" | awk '{print $4}')

    # Scrivi nel file
    echo "complex${i}    ${SCORE}" >> "$OUTFILE"

    # Stampa anche a schermo
    echo "complex${i} score = ${SCORE}"

done

echo ""
echo "Tutti gli score salvati in: $OUTFILE"
