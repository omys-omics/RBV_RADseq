#!/bin/bash

mkdir 01.simlink.bams

# simlink only the good samples from the first plate
for cat_num in FN413 FN414 FN421 FN422 FN428 FN442 FN446 FN450 FN499 FN501 FN534 FN540 FN547 FN550 FN556 FN560 FN572 FN577 FN605 FN608 FN740 MSB147123 MSB147137 MSB147459 MSB147576 MSB147685 MSB147723 MSB149263 MSB149264 MSB149265 MSB157962 MSB198875 MSB199211 UAM20694 UAM20722 UAM20794 UAM23474 UAM48857 UAM50297 UAM50974 UAM59916 WK1514 WK1521;
do
  ln -s /panfs/pfs.local/work/colella/ben/004.RBV.rad/02.plate.1.analyses/03.align.reads/01.basic.fastqs/01.alignments/alignment_$cat_num/$cat_num.sorted.bam ./01.simlink.bams
done

# simlink all the samples from the third plate
ln -s /panfs/pfs.local/work/colella/ben/004.RBV.rad/04.plate.3.analyses/02.align.reads/01.alignments/alignment_*/*sorted.bam ./01.simlink.bams

# simlink only the RBVs from the fourth plate
for cat_num in MSB147138 MSB147139 MSB147140 MSB147141 MSB147142 MSB147449 MSB147451 MSB147455 MSB147457 MSB147463 MSB147466 MSB155686 MSB155687 MSB155694 MSB155710 MSB155711 MSB155712 MSB155730 MSB155731 MSB155751 MSB155752 MSB155753 MSB155862 MSB266678 MSB266686 MSB266688 MSB266704 UAM100098 UAM100102 UAM20791 UAM23201 UAM34216 UAM41893 UAM47614 UAM47615 UAM50299 UAM50300 UAM50301 UAM50302 UAM50303 UAM50304 UAM50305 UAM50483 UAM50489 UAM52520 UAM52534 UAM68435 UAM68436 UAM68437 UAM68438 UAM68442 UAM68443 UAM68444 UAM68445 UAM68446 UAM68452 UAM68454;
do
  ln -s /panfs/pfs.local/work/colella/ben/004.RBV.rad/06.plate.4.analyses/02.align.reads/01.alignments/alignment_$cat_num/$cat_num.sorted.bam ./01.simlink.bams
done

# simlink all the samples from the fifth plate
ln -s /kuhpc/work/colella/ben/004.RBV.rad/08.plate.5.analyses/02.align.reads/01.alignments/alignment_*/*sorted.bam ./01.simlink.bams