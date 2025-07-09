#!/bin/bash
#SBATCH --job-name=index
#SBATCH --partition=colella
#SBATCH --time=6:00:00
#SBATCH --mem=8G
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --output=index.output/out.%j.index.txt
#SBATCH --array=1-265

# set up
module load samtools
bams=`ls /kuhpc/work/colella/ben/004.RBV.rad/09.plate.1-5.analyses/01.align.reads/02.gstacks/01.simlink.bams/*.bam`
mkdir index.output

# use array numbers to get each unique fasta file
bam=$(echo $bams | cut -d" " -f "$SLURM_ARRAY_TASK_ID")
echo $bam

samtools index $bam


