# MD_SIMULATION

## Overview

Running MD simulations on protein-only and protein-ligand systems using a combination of OpenMM and Gromacs.


### protein-only system

`python prepare_protein_only_system.py --input_pdb {protein.pdb}'



    parser.add_argument("--input_pdb", required=True, help="Input .pdb file")
    parser.add_argument("--output_pdb", required=True, help="Output .pdb file")
    parser.add_argument("--output_prefix", required=True, help="Output prefix for Gromacs files")
    parser.add_argument("--forcefield", default="amber14-all.xml", help="OpenMM force field")
    parser.add_argument("--ph", type=float, default=7.0, help="System pH")
    parser.add_argument("--water", default="amber14/tip3pfb.xml", help="OpenMM water model")
    parser.add_argument("--padding", type=float, default=1.0, help="Box padding (nm)")
    parser.add_argument("--ionic_strength", type=float, default=0.15, help="Ionic strength (mM)")
    parser.add_argument("--temperature", type=float, default=300.0, help="Temperature (K)")
    parser.add_argument("--timestep", type=float, default=2.0, help="Time Step (fs)")
    parser.add_argument("--verbose", action="store_true", help="Enable verbose logging")
