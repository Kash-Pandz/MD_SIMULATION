# MD_SIMULATION

## Overview

This repo contains a MD simulation code for protein-only and protein-ligand systems using OpenMM and Gromacs.

## Installation
### OpenMM
1. Install Conda
2. Clone this repository
```bash
git clone https://github.com/Kash-Pandz/MD_SIMULATION.git
cd MD_SIMULATION
```
3. Create the environment
```bash
conda env create -f openmm.yaml
```
4. Activate the environment
```bash
conda activate openmm_sim
```
### GROMACS

GROMACS from source can be compiled inside the `openmm_sim` conda environment.

1. Install required system dependence
```bash
sudo apt-get update
sudo apt-get install -y build-essential cmake git \
    libfftw3-dev libgsl-dev mpi-default-bin mpi-default-dev
```
- `cmake`, `build-essential` for compilation tools
- `libfftw3-dev` for Fourier transforms (required for GROMACS)
- `libgsl-dev` for scientific libraries
- `mpi-default-bin` and `mpi-default-dev` for MPI for parallel runs
2. Download GROMACS from source
```bash
wget ftp://ftp.gromacs.org/gromacs/gromacs-2025.3.tar.gz
tar xfz gromacs-2025.3.tar.gz
cd gromacs-2025.3
mkdir build && cd build
```
3. Configure build with GPU + MPI
```bash
cmake .. \
  -DGMX_BUILD_OWN_FFTW=ON \
  -DGMX_GPU=CUDA \
  -DGMX_MPI=ON \
  -DCMAKE_INSTALL_PREFIX=/usr/local/gromacs \
  -DREGRESSIONTEST_DOWNLOAD=ON
```
4. Compile and install
```bash
make -j$(nproc)
sudo make install
```
5. Load GROMACS environment
```bash
source /usr/local/gromacs/bin/GMXRC
```
6. Test Installation
```bash
gmx --version
```

## Procedure

### 1. protein-only system preparation 

The `prepare_protein_only_system.py` script allows a protein structure to be corrected, minimised and converted to gromacs format. An AMBER FF was used but this can be replaced with a CHARMM FF.

The following command can be used: 

`python prepare_protein_only_system.py --input_pdb INPUT.pdb, --output_pdb OUTPUT.pdb, --output_prefix MD_SIM, --forcefield amber14-all.xml, --ph 7.0, --water amber14/tip3pfb.xml, --padding 1.0, --ionic_strength 0.15, --temperature 300.0, --timestep 2.0, --verbose`

### OR

### 2. protein-ligand system preparation

The `python prepare_protein_ligand_complex.py` script allows a protein-ligand system to be setup, minimised and converted to gromacs format. The ligand AMBER FF parameters are derived using antechamber.

The following command can be used to obtain the AMBER prmtop and inpcrd files:

`python prepare_protein_ligand_complex.py --input_pdb "INPUT_COMPLEX.pdb", --ligand_resname "LIG", --ligand_charge 0, --verbose`

The AMBER prmtop and inpcrd files can then be converted into GROMACS using the acpype command:

`acpype -p SYSTEM.prmtop -x SYSTEM.inpcrd -o gmx`

This will generate SYSTEM.gro, SYSTEM.top, LIG.itp files. 

Alternatively, use `amber_to_gromacs.py` to convert the .prmtop and .inpcrd files to GROMACS format using parmed. N.B Manual inspection and manipulation is required for generated .gro and .top files.


### 3. MD simulation run

MD simulations can be run using GROMACS (`gmx_md_job.sh`) or OpenMM (`openmm_sim.py`). Both scripts follow the same routine 1) energy minimisation 2) NVT equilibration 3) NPT equilibration 4) production MD.

### 4. Analysis

##### Fix PBC effects on md trajectories

To fix PBC effects on the MD trajectories from GROMACS, use the ```fix_pbc.sh``` script.

- protein-only system
``` bash 
fix_gmx_pbc.sh topol.tpr "" traj1.xtc traj2.xtc...
```
- protein-ligand system
```bash
fix_gmx_pbc.sh topol.tpr LIG traj1.xtc traj2.xtc...
```


