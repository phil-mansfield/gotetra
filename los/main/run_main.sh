#!/bin/sh
#SBATCH --job-name=tetra.halo_gen2
#SBATCH --output=old_gtet.out
#SBATCH --error=old_gtet.err
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=6
#SBATCH --time=6:00:00
#SBATCH --mem=16GB
#SBATCH --account=pi-kravtsov
[ -d plots ] || mkdir plots
[ -d text ] || mkdir text
[ -d save ] || mkdir save
# For the love of God, remember to clear save/ if you change sim boxes.
go build main.go && ./main /project/surph/mansfield/data/sheet_segments/Box_L0063_N1024_G0008_CBol/snapdir_100/ /project/surph/diemer/Box_L0063_N1024_CBol/Rockstar/hlists/hlist_1.00000.list plots text save
