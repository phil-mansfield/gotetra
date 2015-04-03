#!/bin/sh
#SBATCH --job-name=tetra.Density
#SBATCH --output=L0125.out
#SBATCH --error=L0125.err
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=14
#SBATCH --time=2:00:00
#SBATCH --mem=4GB
#SBATCH --account=pi-kravtsov

module load go/1.3

# Add this line to turn off bounds checks:
# -gcflags=-B
go build -gcflags=-B && time ./main -Density density.txt L0063/L0063_h5_slice_bounds.txt
