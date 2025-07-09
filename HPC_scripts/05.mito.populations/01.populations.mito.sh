#!/bin/bash
#SBATCH --job-name=populations
#SBATCH --partition=colella
#SBATCH --time=6:00:00
#SBATCH --mem=60G
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --output=out.%j.populations.txt

module load stacks

mkdir 01.all.loci

popmap=popmap.txt
stacks_dir=../04.mito.gstacks
out_dir=01.all.loci
log_file=$out_dir/log.txt


min_samples=0.0 #-R minimum percentage of individuals across populations required to process a locus
min_pops=1 # -p minimum number of populations a locus must be present in to process a locus
min_mac=1 # --min_mac minimum minor allele count required to process a SNP (applied to metapopulation)

populations -t 32 -P $stacks_dir -O $out_dir --popmap $popmap -R $min_samples -p $min_pops --min-mac $min_mac --vcf --phylip-var --phylip-var-all --fasta-loci --fasta-samples &> $log_file
# -t number of threads
# -P path to directory containing Stacks files
# -O output directory
# --popmap path to population map
# --fasta-loci — output locus consensus sequences in FASTA format







