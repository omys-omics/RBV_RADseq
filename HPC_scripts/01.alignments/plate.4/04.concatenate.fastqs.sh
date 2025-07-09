#!/bin/bash
#SBATCH --job-name=concatenate
#SBATCH --partition=sixhour
#SBATCH --time=6:00:00
#SBATCH --mem=250G
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --output=out.%j.cat.txt

# Get names of fastq files
cd BJW-MSG-NdeI-P1-80-UDI-01_S1_R1_001
fastqs=`ls *.fq`
cd ..

# Echo how many there are
num=`echo $fastqs | wc -w`
echo "Processing" $num "fastq files"

# Make output directory for concatenated fastqs
mkdir 04.concatenated.fastqs
cd 04.concatenated.fastqs

# Concatenate matching fastq files
for f in $fastqs
do
  cat ../BJW-MSG-NdeI-P1-80-UDI-01_S1_R1_001/$f ../BJW-MSG-NdeI-P1-80-UDI-02_S2_R1_001/$f ../BJW-MSG-NdeI-P1-80-UDI-03_S3_R1_001/$f> $f
  echo "Processed:" $f
done
