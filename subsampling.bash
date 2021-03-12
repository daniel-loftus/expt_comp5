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

module load Anaconda3/5.1.0

source activate py36

for file in C1_S1_R1_2_001.hicup.G1_all.medium.sorted.filtered.txt ; do 

python subsampling.py -m $file \
-s 0.8642188 \
-o ${file%.txt}.subsampled.txt

done 
