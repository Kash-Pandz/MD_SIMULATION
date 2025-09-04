#!/bin/bash
#SBATCH --job-name=gmx_md
#SBATCH --output=gmx_md_%j.out
#SBATCH --error=gmx_md_%j.err
#SBATCH --time=48:00:00
#SBATCH --partition=nc40ads 
#SBATCH --nodes=1      
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=40
#SBATCH --gpus=1


# Load modules 
module purge
module load gcc cuda/12.4 openmpi
module load gromacs/2025.3-cuda

GMX="srun gmx_mpi"

echo "Job started on $(date)"
echo "Running on $SLURM_NNODES node(s), $SLURM_NTASKS task(s), $SLURM_CPUS_PER_TASK CPU threads, $SLURM_GPUS GPU(s)"

# Energy minimisation
gmx grompp -f em.mdp -c SYSTEM.gro -p SYSTEM.top -o em.tpr
$GMX mdrun -deffnm em -ntmpi 1 -ntomp $SLURM_CPUS_PER_TASK -pin on -gpu_id 0

# NVT equilibration 
gmx grompp -f nvt.mdp -c SYSTEM_EM.gro -r SYSTEM_EM.gro -p SYSTEM.top -o nvt.tpr
$GMX mdrun -deffnm nvt -cpi -append -ntmpi 1 -ntomp $SLURM_CPUS_PER_TASK -pin on -gpu_id 0

# NPT equilibration 
gmx grompp -f npt.mdp -c nvt.gro -r nvt.gro -t nvt.cpt -p SYSTEM.top -o npt.tpr
$GMX mdrun -deffnm npt -cpi -append -ntmpi 1 -ntomp $SLURM_CPUS_PER_TASK -pin on -gpu_id 0

# Production MD
gmx grompp -f md.mdp -c npt.gro -r npt.gro -t npt.cpt -p SYSTEM.top -o md_rep.tpr
$GMX mdrun -deffnm md_rep -cpi -append -ntmpi 1 -ntomp $SLURM_CPUS_PER_TASK -pin on -gpu_id 0

echo "Job finished on $(date)"
