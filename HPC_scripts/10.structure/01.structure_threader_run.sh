#!/bin/bash
#SBATCH --job-name=structure
#SBATCH --partition=colella
#SBATCH --time=400:00:00
#SBATCH --mem=400G
#SBATCH --nodes=1
#SBATCH --ntasks=24
#SBATCH --output=output_structure.%j.txt

out_dir=01.output.1M.K2_8
mkdir $out_dir
dos2unix *.str
structure_threader run -t 24 -i rbv.for.structure.str -o $out_dir --params mainparams -Klist 2 3 4 5 6 7 8 -R 10 -st ~/.local/bin/structure --log TRUE
# -i infile
# -o output directory
# --params File with run parameters
# -t threads
# -Klist '2 4 6' ['2 4 6' ...] List of Ks to calculate
# -R replicates
# -st path to structure
# --log enable logging