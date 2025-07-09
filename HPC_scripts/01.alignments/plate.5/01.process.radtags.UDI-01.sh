#!/bin/bash
#SBATCH --job-name=process_radtags
#SBATCH --partition=sixhour
#SBATCH --time=6:00:00
#SBATCH --mem=186G
#SBATCH --nodes=1
#SBATCH --ntasks=6
#SBATCH --output=out.%j.process_radtags.txt


module load stacks/2.41

mkdir BJW-MSG-NdeI-P1-96-UDI-01_S1_R1_001


process_radtags \
-f /panfs/pfs.local/work/colella/ben/004.RBV.rad/01.raw.reads/plate.5/BJW-MSG-NdeI-P1-96-UDI-01_S1_R1_001.fastq.gz \
-o BJW-MSG-NdeI-P1-96-UDI-01_S1_R1_001/ \
-b plate.5.barcodes.txt \
-e ndeI -r -c -q -y fastq
#-f path to the input file if processing single-end sequences
#-o path to output the processed files
#-b path to a file containing barcodes for this run
#-e specifies the enzyme used, so that process_radtags can look for the correct RAD cutsite
#-c,--clean: clean data, remove any read with an uncalled base
#-q,--quality: discard reads with low quality scores
#-r,--rescue: rescue barcodes and RAD-Tag cut sites
#-y,--out-type: output type, either 'fastq', 'gzfastq', 'fasta', or 'gzfasta' (default: match input type)

# This is the default barcode option, should be the case here
#--inline-null: barcode is inline with sequence, occurs only on single-end read (default).
