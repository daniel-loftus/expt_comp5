#!/bin/bash
#SBATCH -n 16 # Number of cores
#SBATCH -N 1 # Ensure that all cores are on one machine
#SBATCH -t 0-03:00 # Runtime in D-HH:MM
#SBATCH -p whipple # Partition to submit to
#SBATCH --mem=16G # Memory pool for all cores (see also --mem-per-cpu)
#SBATCH -o %x%j.out # File to which STDOUT will be written
#SBATCH -e %x%j.err # File to which STDERR will be written
#SBATCH --mail-type=ALL
#SBATCH --mail-user="danielloftus@g.harvard.edu"

module load samtools/1.5-fasrc02

input_bam=C1_S1_R1_2_001.hicup.bam
output_dir=/n/whipple_lab/Users/dloftus/capture_hic/results_run1/hicup/WASP/alignment_output_N-masked_digestion/bam_to_fastq_for_wasp_single_end_n_masked_digestion

#the -n tag ensures that the names stay the same (otherwise a 1 and 2 are added to individual reads in the pair)
samtools fastq --threads 16 -n -1 ${output_dir}/${input_bam%.bam}.read1.fq -2 ${output_dir}/${input_bam%.bam}.read2.fq -s ${output_dir}/${input_bam%.bam}.singletons.fq $input_bam

#merge the output fastq file s 
cat ${output_dir}/${input_bam%.bam}.read1.fq ${output_dir}/${input_bam%.bam}.read2.fq ${output_dir}/${input_bam%.bam}.singletons.fq > ${output_dir}/${input_bam%.bam}.all_reads.fq