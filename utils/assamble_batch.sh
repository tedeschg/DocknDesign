#!/bin/bash
set -e

# Numero totale di cartelle
N=20

# Base directory (directory da cui lanci lo script)
BASE_DIR="$(pwd)"

# Percorsi fissi (tool e proteina restano assoluti)
PY_SCRIPT="/home/tedeschg/prj/dockndesign/utils/assemble_complex.py"
PROTEIN="/home/tedeschg/prj/dockndesign/pdbs/pdk1/pdk1.pdb"

# Output directory (relativa alla current directory)
OUTDIR="${BASE_DIR}/complexes"
mkdir -p "$OUTDIR"

for i in $(seq 1 $N); do
    DIR="${BASE_DIR}/test_multi_${i}"
    LIG="${DIR}/rank1.sdf"
    OUT="${OUTDIR}/complex${i}.pdb"

    echo "==============================="
    echo " Processing $LIG"
    echo " Output → $OUT"
    echo "==============================="

    # Check di sicurezza (consigliato)
    if [[ ! -f "$LIG" ]]; then
        echo "WARNING: $LIG not found, skipping."
        continue
    fi

    python "$PY_SCRIPT" \
        -p "$PROTEIN" \
        -l "$LIG" \
        -o "$OUT"
done
