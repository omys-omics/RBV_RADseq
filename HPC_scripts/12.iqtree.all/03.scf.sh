#!/bin/bash
#SBATCH --job-name=scf
#SBATCH --partition=colella
#SBATCH --time=100:00:00
#SBATCH --mem=200G
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --output=output_scf.%j.txt

module load iqtree

iqtree2 -te rbv.phy.treefile -s rbv.phy --scfl 1000 -nt 8 -redo
#-s specifies the input sequence data
#--scfl specifies the number of quartets (randomly sampled around each internal branch) for computing sCF
#  recommended at least 100
#-te input species tree
#-nt AUTO allows the program to use the optimal number of threads (8 specified here)