#!/bin/bash
#SBATCH --job-name=iqtree
#SBATCH --partition=colella
#SBATCH --time=100:00:00
#SBATCH --mem=100G
#SBATCH --nodes=1
#SBATCH --ntasks=12
#SBATCH --output=output_iqtree.%j.txt

module load iqtree

# remove extra line that stacks adds to phylip file
sed '/^#/ d' 01.whitelist.rbv/populations.all.phylip > rbv.phy

iqtree2 -s rbv.phy -m MFP -bb 1000 -nt AUTO
#-s specifies the input sequence data
#-m MFP specifies to perform model testing and use the best model of sequence evolution
#-bb specifies performing 1000 ultrafast bootstraps to assess support
#-nt AUTO allows the program to use the optimal number of threads (15 specified here)