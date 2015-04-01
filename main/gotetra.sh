#!/bin/sh
#SBATCH --job-name=tetra.Density
#SBATCH --output=L0250_S4.out
#SBATCH --error=L0250_S4.err
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=14
#SBATCH --time=1:00:00
#SBATCH --mem=20GB
#SBATCH --account=pi-kravtsov

module load go/1.3

# Add this line to turn off bounds checks:
# -gcflags=-B
go build && time ./main -Density L0125/density.txt
