#!/bin/bash

# Activate AMBER conda environment
eval "$(conda shell.bash hook)"
conda activate ambertools_build

# Optional: set CUDA if needed
export CUDA_VISIBLE_DEVICES=0
export CUDA_HOME=/usr/local/cuda

# Loop over all directories
for dir in pos_*; do
    echo "Processing directory: $dir"
    cd "$dir" || continue

    pdb="model_0.pdb"
    if [ ! -f "$pdb" ]; then
        echo "PDB not found, skipping $dir"
        cd ..
        continue
    fi

    # Use awk to split the PDB file. This is the final fix.
    # It defines the receptor as any line where the chain ID (column 5) starts with 'A'.
    awk '$5 ~ /^A/ {print}' "$pdb" > receptor.pdb
    # It defines the ligand as any line where the chain ID (column 5) starts with 'B'.
    awk '$5 ~ /^B/ {print}' "$pdb" > ligand.pdb

    # Run tleap with inline input
    tleap -I ~/Downloads/amber24_src/dat/leap/cmd \
          -I ~/Downloads/amber24_src/dat/leap/parm \
          -I ~/Downloads/amber24_src/dat/leap/lib \
          -f - <<EOF
source leaprc.protein.ff14SB
receptor = loadpdb receptor.pdb
saveamberparm receptor receptor.prmtop receptor.inpcrd
ligand = loadpdb ligand.pdb
saveamberparm ligand ligand.prmtop ligand.inpcrd
complex = loadpdb $pdb
saveamberparm complex complex.prmtop complex.inpcrd
quit
EOF

    # Check LEaP success
    if [ ! -f complex.inpcrd ]; then
        echo "LEaP failed in $dir" > Static_AMBER.txt
        cd ..
        continue
    fi

    # Minimization input
    cat > min.in <<EOF
Minimization
&cntrl
imin=1,
maxcyc=500,
ncyc=250,
cut=10.0,
ntb=0,
/
EOF

    sander -O -i min.in -p complex.prmtop -c complex.inpcrd -o min.out -r min.rst -ref complex.inpcrd

    # MMPBSA input
    cat > mmpbsa.in <<EOF
&general
startframe=1,
endframe=1,
interval=1,
verbose=1,
keep_files=0,
/
&gb
igb=5,
saltcon=0.150
/
EOF

    # Run MMPBSA.py from conda environment
    MMPBSA.py -O \
        -i mmpbsa.in \
        -cp complex.prmtop \
        -rp receptor.prmtop \
        -lp ligand.prmtop \
        -y min.rst > Static_AMBER.log 2>&1

    # Extract result
    if [ -f FINAL_RESULTS_MMPBSA.dat ]; then
        grep 'DELTA TOTAL' FINAL_RESULTS_MMPBSA.dat | awk '{print $3}' > Static_AMBER.txt
    else
        echo "MMPBSA failed in $dir" > Static_AMBER.txt
    fi

    cd ..
done

echo "All evaluations complete."
