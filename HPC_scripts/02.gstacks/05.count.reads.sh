#!/bin/bash


module load samtools


bams=`ls 04.simlink.bams/`

for bam in $bams
do
  
  id=`basename $bam | cut -f1 -d"."`
  
  reads=`samtools flagstat 04.simlink.bams/$bam | head -1 | cut -f1 -d" "`

  printf "$id\t$reads\n" >> reads.txt
done





