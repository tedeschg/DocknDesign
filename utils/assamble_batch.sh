#!/bin/bash

# Numero totale di cartelle
N=20

# Percorsi fissi
PY_SCRIPT="/home/tedeschg/prj/dockndesign/utils/assemble_complex.py"
PROTEIN="/home/tedeschg/prj/dockndesign/pdbs/pdk1/pdk1.pdb"

for i in $(seq 1 $N); do
    DIR="test_multi_${i}"
    LIG="${DIR}/rank1.sdf"
    OUT="complex${i}.pdb"

    echo "Processing $LIG → $OUT"

    python "$PY_SCRIPT" \
        -p "$PROTEIN" \
        -l "$LIG" \
        -o "$OUT"
done
