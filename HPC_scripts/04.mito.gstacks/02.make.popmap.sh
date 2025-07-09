#!/bin/bash

module load samtools

bams=`ls /kuhpc/work/colella/ben/004.RBV.rad/09.plate.1-5.analyses/01.align.reads/04.mito.gstacks/01.mito.bams/*.bam`

# make a popmap of only samples that have sufficient data (measured by the number of lines in the sam file)
for bam in $bams;
do
  num_lines=`samtools view $bam | wc -l`
  if [[ $num_lines -gt 31 ]]
  then
    name=`basename -a $bam | cut -f1 -d"."`
    printf "$name\tpop1\n" >> popmap.txt
  fi
done

