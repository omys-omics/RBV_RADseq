#!/bin/bash
#SBATCH --job-name=mito
#SBATCH --partition=colella
#SBATCH --time=6:00:00
#SBATCH --mem=8G
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --output=mito.output/out.%j.mito.txt
#SBATCH --array=1-265

# set up
module load samtools
bams=`ls /kuhpc/work/colella/ben/004.RBV.rad/09.plate.1-5.analyses/01.align.reads/02.gstacks/01.simlink.bams/*.bam`
mkdir 01.mito.bams
mkdir mito.output

# use array numbers to get each unique fasta file
bam=$(echo $bams | cut -d" " -f "$SLURM_ARRAY_TASK_ID")

# unique name
name=`basename $bam | cut -f1 -d"."`

# subset bam file for reads aligning to mitochondrion
echo $name
samtools view -b $bam NC_024538.1 | samtools sort > 01.mito.bams/$name.mito.bam
#NC_024538.1 is the mitochondrion
# -b output in bam format





