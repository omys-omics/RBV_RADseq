#!/bin/bash
#SBATCH --job-name=bwa
#SBATCH --partition=colella
#SBATCH --time=23:00:00
#SBATCH --mem=250G
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --output=out.%j.bwa.txt


module load bwa
module load samtools

#Get fastqs
fastqs=`ls /kuhpc/work/colella/ben/004.RBV.rad/09.plate.1-5.analyses/01.align.reads/06.mito.rutilus.gstacks/01.simlink.fqs/*.fq`

names=`basename -a $fastqs | cut -f1 -d"."`
num=`echo $names | wc -w`
echo $names
echo "Processing" $num "fastq files"

for fastq in $fastqs;
do
    
#Make sure to start in correct directory
    cd /kuhpc/work/colella/ben/004.RBV.rad/09.plate.1-5.analyses/01.align.reads/06.mito.rutilus.gstacks/02.alignments

#Get sample name
    name=`basename -a $fastq | cut -f1 -d"."`

#Only do identified rutilus mitos
    rutilus=`grep $name ../rutilus.mt.txt | wc -l`
    if [[ $rutilus -eq 1 ]]
    then
        echo $name >> ../02.aligned.txt

        #Make directory for alignment
        mkdir mito_alignment_$name
        cd mito_alignment_$name

        #Align raw reads to reference genome and convert to .bam
        bwa mem -t 32 /kuhpc/work/colella/ben/004.RBV.rad/09.plate.1-5.analyses/01.align.reads/06.mito.rutilus.gstacks/rutilus_mitogenome/rutilus_mitogenome $fastq | samtools view -b -@32 -o $name.raw.bam
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
    
    else
        echo $name >> ../02.not.aligned.txt
    fi
    
done

