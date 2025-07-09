#!/bin/bash
#SBATCH --job-name=nuc
#SBATCH --partition=sixhour
#SBATCH --time=6:00:00
#SBATCH --mem=60G
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --output=out.%j.nuc.txt

module load vcftools

vcftools --vcf 01.all.loci/populations.snps.vcf --not-chr NC_024538.1 --recode --recode-INFO-all --out all.nuc

# --vcf input vcf file
# --not-chr exclude sites with identifiers matching <chromosome>
# --out prefix of output file
# --recode generate a new vcf file with the filtering options
# --recode-INFO-all keep INFO fields




