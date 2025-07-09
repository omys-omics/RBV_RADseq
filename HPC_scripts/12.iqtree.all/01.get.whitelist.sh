#!/bin/bash
#SBATCH --job-name=whitelist
#SBATCH --partition=sixhour
#SBATCH --time=6:00:00
#SBATCH --mem=250G
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --output=%j.output_whitelist.txt

module load stacks

popmap=iqtree.popmap.rbv.txt
stacks_dir=../02.gstacks
out_dir=01.whitelist.rbv
log_file=$out_dir/log.txt

mkdir $out_dir

min_samples=0.0 #-R minimum percentage of individuals across populations required to process a locus
min_pops=1 # -p minimum number of populations a locus must be present in to process a locus
min_mac=1 # --min_mac minimum minor allele count required to process a SNP (applied to metapopulation)

populations -t 12 -P $stacks_dir -O $out_dir --popmap $popmap --whitelist whitelist.rbv.txt -R $min_samples -p $min_pops --min-mac $min_mac --phylip-var-all &> $log_file
# -t number of threads
# -P path to the directory containing the Stacks files
# -O output directory
# --popmap path to population map
# --phylip-var-all flag indicates to output the phylip including invariant sites (best for phylogenetic reconstruction)











