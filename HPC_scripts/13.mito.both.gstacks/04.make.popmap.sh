#!/bin/bash

bams=`ls ./03.simlink.bams/`

for bam in $bams;
do
  name=`echo $bam | cut -f1 -d"."`
  printf "$name\tpop1\n" >> popmap.txt
done

