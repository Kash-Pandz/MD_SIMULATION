#!/bin/bash

TPR=$1
LIG=$2
shift 2

for XTC in "$@"; do
    BASENAME=$(basename "$XTC" .xtc)
    
    if [[ -z "$LIG" ]]; then
        echo "Protein System" | gmx trjconv -s "$TPR" -f "$XTC" \
            -o "${BASENAME}_pbc.xtc" \
            -pbc nojump -pbc mol -center -ur compact -quiet
    else
        gmx make_ndx -f "$TPR" -o index.ndx << EOF > /dev/null 2>&1
r $LIG | Protein
name 3 Prot_Lig
q
EOF
        echo "Prot_Lig System" | gmx trjconv -s "$TPR" -f "$XTC" \
            -o "${BASENAME}_pbc.xtc" \
            -pbc nojump -pbc mol -center -ur compact -n index.ndx -quiet
        rm -f index.ndx
    fi
done
