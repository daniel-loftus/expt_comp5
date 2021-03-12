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


module load samtools/0.1.19-fasrc02
module load bedtools2/2.26.0-fasrc01


input_bam=C1_S1_R1_2_001.hicup.bam

echo The first count is:
samtools view $input_bam | wc -l

region_bed=/n/whipple_lab/Users/dloftus/capture_hic/coords_mm10_peg13_locus.bed

bedtools intersect -u -bed -a $input_bam -b $region_bed > ${input_bam}_on_target.bed 

echo On target count is: 
cut -f4 ${input_bam}_on_target.bed | cut -f1 -d"/" | sort | uniq | wc -l