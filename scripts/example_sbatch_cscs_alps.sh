#!/bin/bash
#SBATCH --job-name=psdf
#SBATCH --account=c40
#SBATCH --nodes=4
#SBATCH --ntasks-per-node=4
#SBATCH --cpus-per-task=72
#SBATCH --time=00:20:00
#SBATCH --output=slurm-psdf-%j.out
#SBATCH --error=slurm-psdf-%j.err
#SBATCH --exclusive
#SBATCH --partition=normal

set -euo pipefail

export MPICH_GPU_SUPPORT_ENABLED=1
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export OMP_PROC_BIND=true

echo "#---------------#"
date
echo "#---------------#"

# Compile with MPI and OpenMP for hybrid parallelism
srun ./bin/psdf_exe pump_surface_1/ pump_1 out.raw

echo "#---------------#"
date
echo "#---------------#"