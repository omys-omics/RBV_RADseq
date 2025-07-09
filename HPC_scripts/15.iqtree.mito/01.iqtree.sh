#!/bin/bash
#SBATCH --job-name=iqtree
#SBATCH --partition=sixhour
#SBATCH --time=6:00:00
#SBATCH --mem=60G
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --output=output_iqtree.%j.txt

module load iqtree

iqtree2 -s mito.SNPs.both.txt -m GTR+ASC -bb 10000 -nt AUTO -o UAM23201
#-s specifies the input sequence data
#-m GTR+ASC specifies to use GTR model and account for ascertain bias (because this is SNP data that does not contain constant sites)
#-bb specifies performing 1000 ultrafast bootstraps to assess support
#-nt AUTO allows the program to use the optimal number of threads (15 specified here)
#-o specify outgroup