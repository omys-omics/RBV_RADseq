#!/bin/bash
#SBATCH --job-name=dsuite
#SBATCH --partition=sixhour
#SBATCH --time=1:00:00
#SBATCH --mem=90G
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --output=out.%j.dsuite.txt


sites=`ls dsuite_SEAK/*.vcf | cut -d"." -f2`

for s in $sites
do
  sed "s/SITE/$s/g" guide.tree.txt > dsuite_SEAK/dsuite.$s.guide.tree.txt
  /kuhpc/work/colella/lab_software/Dsuite/Build/Dsuite Dtrios -k 1000 -o dsuite_SEAK/ -n $s -t dsuite_SEAK/dsuite.$s.guide.tree.txt dsuite_SEAK/dsuite.$s.vcf dsuite_SEAK/dsuite.$s.popmap.txt


done

# -k number of jackknife blocks to divide the dataset into
# -n run name
# -t guide tree (newick format)
# <INPUT_VCF>
# <INPUT POPMAP>
