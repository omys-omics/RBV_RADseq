#!/bin/bash
#SBATCH --job-name=bwa
#SBATCH --partition=colella
#SBATCH --time=47:00:00
#SBATCH --mem=500G
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --output=out.%j.bwa.txt


module load bwa
module load samtools

#Get fastqs from fourth plate
fastqs=`ls /panfs/pfs.local/work/colella/ben/004.RBV.rad/08.plate.5.analyses/01.process.radtags/04.concatenated.fastqs/*.fq`

names=`basename -a $fastqs | cut -f1 -d"."`
num=`echo $names | wc -w`
echo $names
echo "Processing" $num "fastq files"

for fastq in $fastqs;
do

#Get sample name
    name=`basename -a $fastq | cut -f1 -d"."`

#Make sure to start in correct directory
    cd /panfs/pfs.local/work/colella/ben/004.RBV.rad/08.plate.5.analyses/02.align.reads/01.alignments/

#Make directory for alignment
    mkdir alignment_$name
    cd alignment_$name

#Align raw reads to reference genome and convert to .bam
    bwa mem -t 32 /panfs/pfs.local/work/colella/ben/004.RBV.rad/00.myodes.glareolus.ref.genome/indexed/GCF_902806735.1_Bank_vole1_10x_genomic $fastq | samtools view -b -@32 -o $name.raw.bam
      # -t number of threads to use
      # -b output in bam format
      # -@ number of threads to use
      # -o output file name
    echo "aligned raw reads for $name"
    date
    
#Sort alignments by leftmost coordinates
    samtools sort -@32 $name.raw.bam -o $name.sorted.bam
    # -o output file
    echo "sorted bam for $name"
    date
    wait
done

