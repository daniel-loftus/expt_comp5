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

for file in C1_S1_R1_2_001.hicup.G1_all.medium.sorted.txt C1_S1_R1_2_001.hicup.G2_all.medium.sorted.txt; do

awk '{ if($3=="chr15" && $4>72460000 && $4<73200000 && $7=="chr15" && $8>72460000 && $8<73200000) print}' $file > whole_locus.temp 

awk '{ if( ($4<72963827 || $4>72980478) && ($8<72963827 || $8>72980478) ) print}' whole_locus.temp > ${file%.txt}.filtered.txt

rm whole_locus.temp 

done 

