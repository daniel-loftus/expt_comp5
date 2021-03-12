#!/bin/bash
#SBATCH -n 1 # Number of cores
#SBATCH -N 1 # Ensure that all cores are on one machine
#SBATCH -t 0-03:00 # Runtime in D-HH:MM
#SBATCH -p whipple # Partition to submit to
#SBATCH --mem=16G # Memory pool for all cores (see also --mem-per-cpu)
#SBATCH -o %x%j.out # File to which STDOUT will be written
#SBATCH -e %x%j.err # File to which STDERR will be written
#SBATCH --mail-type=ALL
#SBATCH --mail-user="danielloftus@g.harvard.edu"

#install necessary packages 
module load samtools/1.5-fasrc02
module load jdk/10.0.1-fasrc01 
PATH=$PATH:/n/whipple_lab/Users/dloftus/software/juicer/juicer-master/SLURM/scripts



for file in C1_S1_R1_2_001.hicup.G1_all.bam C1_S1_R1_2_001.hicup.G2_all.bam; do



#convert to sam file 
samtools view -o ${file%.bam}.sam $file



#generate medium file 
mysam=${file%.bam}.sam
awk '{if(NR%2==0) print $1, $3, $4, $5}' $mysam > evens_temp.txt
awk '{if(NR%2!=0) print $1, $3, $4, $5}' $mysam > odds_temp.txt

awk 'BEGIN{OFS="\t"} FNR==NR { a[FNR""] = $0; next } { print a[FNR""], $0 }' evens_temp.txt odds_temp.txt > combined_temp.txt 

awk 'BEGIN{OFS="\t"} {print $1, 0, "chr"$2, $3, 0, 0, "chr"$6, $7, 1, $4, $8}' combined_temp.txt > ${mysam%.sam}.medium.txt

sort -k3,3 -k4,4n ${mysam%.sam}.medium.txt > ${mysam%.sam}.medium.sorted.txt

rm evens_temp.txt
rm odds_temp.txt
rm combined_temp.txt 


#convert to .hic
input_file=${mysam%.sam}.medium.sorted.txt
out_name=${file%.bam}.hic
juicer_tools pre -d $input_file $out_name mm10


done 
