#!/bin/bash
#SBATCH --job-name=dsuite
#SBATCH --partition=sixhour
#SBATCH --time=1:00:00
#SBATCH --mem=90G
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --output=out.%j.dsuite.txt


sites=`ls dsuite/*.vcf | cut -d"." -f2`

for s in $sites
do
  sed "s/SITE/$s/g" guide.tree.txt > dsuite/dsuite.$s.guide.tree.txt
  /kuhpc/work/colella/lab_software/Dsuite/Build/Dsuite Dtrios -k 1000 -o dsuite/ -n $s -t dsuite/dsuite.$s.guide.tree.txt dsuite/dsuite.$s.vcf dsuite/dsuite.$s.popmap.txt


done



#/kuhpc/work/colella/lab_software/Dsuite/Build/Dsuite Dtrios -k 1000 -n test -t guide.tree.txt dsuite.vcf dsuite.popmap.txt

# -k number of jackknife blocks to divide the dataset into
# -n run name
# -t guide tree (newick format)
# <INPUT_VCF>
# <INPUT POPMAP>
