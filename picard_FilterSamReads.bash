#!/bin/bash
#SBATCH -n 1 # Number of cores
#SBATCH -N 1 # Ensure that all cores are on one machine
#SBATCH -t 0-01:00 # Runtime in D-HH:MM
#SBATCH -p whipple # Partition to submit to
#SBATCH --mem=4G # Memory pool for all cores (see also --mem-per-cpu)
#SBATCH -o %x%j.out # File to which STDOUT will be written
#SBATCH -e %x%j.err # File to which STDERR will be written
#SBATCH --mail-type=ALL
#SBATCH --mail-user="danielloftus@g.harvard.edu"

module load samtools/1.5-fasrc02


module load Java/1.8

input_bam_picard=C1_S1_R1_2_001.hicup.G1_all.region.sorted.bam

java -jar /n/whipple_lab/Users/dloftus/software/picard/picard.jar FilterSamReads \
--FILTER includeReadList \
--READ_LIST_FILE region_of_no_snps_ids.txt \
--INPUT $input_bam_picard \
--OUTPUT ${input_bam_picard%.bam}.filtered_for_no_snps_reads.bam 


input_bam_samtools2=${input_bam%.bam}.filtered_for_no_snps_reads.bam 

samtools sort $input_bam_samtools2 -o ${input_bam%.bam}.filtered_for_no_snps_reads.sorted.bam

samtools index ${input_bam%.bam}.filtered_for_no_snps_reads.sorted.bam ${input_bam%.bam}.filtered_for_no_snps_reads.sorted.bam.bai

