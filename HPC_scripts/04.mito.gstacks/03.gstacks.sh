#!/bin/bash
#SBATCH --job-name=gstacks
#SBATCH --partition=colella
#SBATCH --time=6:00:00
#SBATCH --mem=120G
#SBATCH --nodes=1
#SBATCH --ntasks=2
#SBATCH --output=out.%j.gstacks.txt

module load stacks

gstacks --threads 32 -I ./01.mito.bams -M popmap.txt -S .mito.bam -O ./
# -I directory containing input bam files
# -M path to pop map
# -S suffix expected for the bam files
# -O output directory
