#!/bin/bash
#SBATCH --job-name=eems
#SBATCH --partition=bi
#SBATCH --time=24:00:00
#SBATCH --mem=120G
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --output=output_eems_rutilus.%j.txt

module load eems

mkdir rep100.rutilus
mkdir rep200.rutilus
mkdir rep400.rutilus

runeems_snps --params /kuhpc/work/colella/ben/004.RBV.rad/09.plate.1-5.analyses/01.align.reads/18.rutilus.eems/params.ndemes100.rutilus.txt
runeems_snps --params /kuhpc/work/colella/ben/004.RBV.rad/09.plate.1-5.analyses/01.align.reads/18.rutilus.eems/params.ndemes200.rutilus.txt
runeems_snps --params /kuhpc/work/colella/ben/004.RBV.rad/09.plate.1-5.analyses/01.align.reads/18.rutilus.eems/params.ndemes400.rutilus.txt