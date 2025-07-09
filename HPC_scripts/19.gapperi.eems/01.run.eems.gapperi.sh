#!/bin/bash
#SBATCH --job-name=eems
#SBATCH --partition=colella
#SBATCH --time=24:00:00
#SBATCH --mem=200G
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --output=output_eems_gapperi.%j.txt

module load eems

mkdir rep100.gapperi
mkdir rep200.gapperi
mkdir rep400.gapperi

runeems_snps --params /kuhpc/work/colella/ben/004.RBV.rad/09.plate.1-5.analyses/01.align.reads/19.gapperi.eems/params.ndemes100.gapperi.txt
runeems_snps --params /kuhpc/work/colella/ben/004.RBV.rad/09.plate.1-5.analyses/01.align.reads/19.gapperi.eems/params.ndemes200.gapperi.txt
runeems_snps --params /kuhpc/work/colella/ben/004.RBV.rad/09.plate.1-5.analyses/01.align.reads/19.gapperi.eems/params.ndemes400.gapperi.txt