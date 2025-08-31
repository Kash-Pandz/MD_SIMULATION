#!/bin/bash
#SBATCH --job-name=gmx_md
#SBATCH --output=gmx_md_%j.out
#SBATCH --error=gmx_md_%j.err
#SBATCH --time=48:00:00
#SBATCH --partition=hc24rs   
#SBATCH --nodes=2      
#SBATCH --ntasks-per-node=44 
#SBATCH --cpus-per-task=1
#SBATCH --exclusive

# Load modules 
module purge
module load gcc mpi
module load gromacs/2025.3

GMX="srun gmx_mpi"

echo "Job started on $(date)"
echo "Running on $SLURM_NNODES nodes with $SLURM_NTASKS tasks"

# NVT equilibration 
gmx grompp -f nvt.mdp -c SYSTEM_EM.gro -r SYSTEM_EM.gro -p SYSTEM.top -o nvt.tpr
$GMX mdrun -deffnm nvt

# NPT equilibration 
gmx grompp -f npt.mdp -c nvt.gro -r nvt.gro -t nvt.cpt -p SYSTEM.top -o npt.tpr
$GMX mdrun -deffnm nvt

# Production MD
gmx grompp -f md.mdp -c npt.gro -r npt.gro -t npt.cpt -p SYSTEM.top -o md_rep.tpr
$GMX mdrun -deffnm md_rep -cpi -append

echo "Job finished on $(date)"
