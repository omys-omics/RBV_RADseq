#!/bin/bash

fqs=`ls ./02.alignments/`

for fq in $fqs;
do
  name=`echo $fq | cut -f3 -d"_"`
  printf "$name\tpop1\n" >> popmap.txt
done

