#!/bin/bash

bams=`ls ./01.simlink.bams/`

names=`basename -a $bams | cut -f1 -d"."`

for name in $names;
do
  printf "$name\tpop1\n" >> popmap.txt
done

