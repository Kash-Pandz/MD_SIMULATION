# MD_SIMULATION

## Overview

Running MD simulations on protein-only and protein-ligand systems using a combination of OpenMM and Gromacs.


### protein-only system

The `prepare_protein_only_system.py` script allows a protein structure to be corrected, minimised and converted to gromacs format. An AMBER FF was used but this can be replaced with a CHARMM FF.

The following command can be used: 

`python prepare_protein_only_system.py --input_pdb "INPUT.pdb" --output_pdb "OUTPUT.pdb", --output_prefix "MD_SIM" --forcefield "amber14-all.xml", --ph 7.0, --water "amber14/tip3pfb.xml", --padding 1.0, --ionic_strength 0.15, --temperature 300.0, --timestep 2.0 --verbose`


### protein-ligand system

The `python prepare_protein_ligand_complex.py` script allows a protein-ligand system to be setup, minimised and converted to gromacs format. The ligand AMBER FF parameters are derived using antechamber.

The following command can be used:

`python prepare_protein_ligand_complex.py --input_pdb "INPUT_COMPLEX.pdb", --ligand_resname "LIG", --ligand_charge 0, --do_minimise, --gmx_gro "system_EM.gro", --gmx_top "system_EM.top", --min_pdb "EM.pdb" --verbose`


