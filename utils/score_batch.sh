#!/bin/bash
set -e

N=20
BASE_DIR="$(pwd)"

SCORE_SCRIPT="/home/tedeschg/prj/dockndesign/utils/lmpnn_score.py"
PDB="/home/tedeschg/prj/dockndesign/pdbs/pdk1/pdk1.pdb"
POCKET="${BASE_DIR}/pocket_residues.txt"

OUTFILE="${BASE_DIR}/scores_summary.txt"
LOGFILE="${BASE_DIR}/score_with_mutations.log"

echo "Complex    Score" > "$OUTFILE"
echo "-----------------" >> "$OUTFILE"
: > "$LOGFILE"   # svuota il log a inizio run

for i in $(seq 1 $N); do
    FOLDER="${BASE_DIR}/lmpnn_results/lmpnn${i}"
    FASTA="${FOLDER}/seqs/complex${i}.fa"

    echo "Scoring: $FASTA"

    if [[ ! -f "$FASTA" ]]; then
        echo "WARNING: $FASTA not found, skipping."
        continue
    fi

    OUTPUT=$(python "$SCORE_SCRIPT" \
        -f "$FASTA" \
        -p "$PDB" \
        -r "$POCKET" \
        --show-mutations
    )

    {
      echo "===== complex${i} ====="
      echo "$OUTPUT"
      echo
    } >> "$LOGFILE"

    SCORE=$(echo "$OUTPUT" | awk '/FINAL SCORE/ {print $NF}')
    echo "complex${i}    ${SCORE}" >> "$OUTFILE"
    echo "complex${i} score = ${SCORE}"
done
