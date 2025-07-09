#!/bin/bash
#SBATCH --job-name=bwa
#SBATCH --partition=colella
#SBATCH --time=47:00:00
#SBATCH --mem=500G
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --output=out.%j.bwa.txt


module load bwa
module load samtools

#Get fastqs from fourth plate
fastqs=`ls /panfs/pfs.local/work/colella/ben/004.RBV.rad/06.plate.4.analyses/01.process.radtags/04.concatenated.fastqs/*.fq`

#only align the RBVs on the fourth plate, not the shrews
RBVs=`echo MSB147138 MSB147139 MSB147140 MSB147141 MSB147142 MSB147449 MSB147451 MSB147455 MSB147457 MSB147463 MSB147466 MSB155686 MSB155687 MSB155694 MSB155710 MSB155711 MSB155712 MSB155730 MSB155731 MSB155751 MSB155752 MSB155753 MSB155862 MSB266678 MSB266686 MSB266688 MSB266704 UAM100098 UAM100102 UAM20791 UAM23201 UAM34216 UAM41893 UAM47614 UAM47615 UAM50299 UAM50300 UAM50301 UAM50302 UAM50303 UAM50304 UAM50305 UAM50483 UAM50489 UAM52520 UAM52534 UAM68435 UAM68436 UAM68437 UAM68438 UAM68442 UAM68443 UAM68444 UAM68445 UAM68446 UAM68452 UAM68454`


names=`basename -a $fastqs | cut -f1 -d"."`
num=`echo $names | wc -w`
echo $names
echo "Processing" $num "fastq files"

for fastq in $fastqs;
do

#Get sample name
    name=`basename -a $fastq | cut -f1 -d"."`

#Check that sample is RBV, if not then skip
    if [[ $RBVs =~ $name ]]
    then 
      echo "$name is an RBV!"
    else
      echo "$name is not an RBV :("
      continue
    fi

#Make sure to start in correct directory
    cd /panfs/pfs.local/work/colella/ben/004.RBV.rad/06.plate.4.analyses/02.align.reads/01.alignments/

#Make directory for alignment
    mkdir alignment_$name
    cd alignment_$name

#Align raw reads to reference genome and convert to .bam
    bwa mem -t 32 /panfs/pfs.local/work/colella/ben/004.RBV.rad/00.myodes.glareolus.ref.genome/indexed/GCF_902806735.1_Bank_vole1_10x_genomic $fastq | samtools view -b -@32 -o $name.raw.bam
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
done

