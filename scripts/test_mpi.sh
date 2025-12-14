#!/bin/bash
#SBATCH --job-name=test_mpi
#SBATCH --account=c40
#SBATCH --nodes=4
#SBATCH --ntasks-per-node=72
#SBATCH --cpus-per-task=1
#SBATCH --time=00:05:00
#SBATCH --output=slurm-%j.out
#SBATCH --error=slurm-%j.err
#SBATCH --exclusive
#SBATCH --hint=nomultithread
#SBATCH --partition=debug

set -euo pipefail

workdir="$PWD/temp_test_mpi_"
mkdir -p "$workdir"
cd "$workdir"

# Write the source
cat > test_mpi.c <<'EOF'
#include <mpi.h>
#include <stdio.h>

int main(int argc, char *argv[]) {
  MPI_Init(&argc, &argv);

  int rank, size;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  for (int r = 0; r < size; r++) {
    if (r == rank) {
      printf("%d/%d\n", rank + 1, size);
      fflush(stdout);
    }
    MPI_Barrier(MPI_COMM_WORLD);
  }

  return MPI_Finalize();
}
EOF

export MPICH_GPU_SUPPORT_ENABLED=1

# Compile (use the MPI wrapper)
mpicc -O2 -std=c11 -Wall -Wextra -o test_mpi test_mpi.c

# Run
srun --cpu-bind=cores ./test_mpi

rm -r $workdir
