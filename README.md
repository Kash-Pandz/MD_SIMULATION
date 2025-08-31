# MD_SIMULATION

## Overview

Running MD simulations on protein-only and protein-ligand systems using a combination of OpenMM and Gromacs.


### protein-only system

The `prepare_protein_only_system.py` script allows a protein structure to be corrected, minimised and converted to gromacs format. An AMBER FF was used but this can be replaced with a CHARMM FF.

`python prepare_protein_only_system.py --input_pdb "INPUT.pdb" --output_pdb "OUTPUT.pdb", --output_prefix "MD_SIM" --forcefield "amber14-all.xml", --ph 7.0, --water "amber14/tip3pfb.xml", --padding 1.0, --ionic_strength 0.15, --temperature 300.0, --timestep 2.0 --verbose`
