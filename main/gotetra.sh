#!/bin/sh
#SBATCH --job-name=tetra.Density
#SBATCH --output=Density.out
#SBATCH --error=Density.err
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=14
#SBATCH --time=1:00:00
#SBATCH --mem=8GB
#SBATCH --account=pi-kravtsov

module load go/1.3

#go build && time ./gotetra -CreateSheet -Cells 8 /project/surph/diemer/Box_L0125_N1024_CBol/Snaps/snapdir_100/snapshot_100.* /project/surph/mansfield/data/sheet_segments/

# Add this line to turn off bounds checks:
# -gcflags=-B
go build -gcflags=-B && time ./main -BoundedDensity -Points 1000 -Cells 5000 -MinSheet 0 -MaxSheet 511 -BoundsFile bounds.txt /project/surph/mansfield/data/tetra_fields/sheet_segments/Box_L0125_N1024_G0008_CBol/snapdir_100/ .
