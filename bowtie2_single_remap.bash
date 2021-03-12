#!/bin/bash
#SBATCH -n 32 # Number of cores
#SBATCH -N 1 # Ensure that all cores are on one machine
#SBATCH -t 1-00:00 # Runtime in D-HH:MM
#SBATCH -p whipple # Partition to submit to
#SBATCH --mem=32G # Memory pool for all cores (see also --mem-per-cpu)
#SBATCH -o %x%j.out # File to which STDOUT will be written
#SBATCH -e %x%j.err # File to which STDERR will be written
#SBATCH --mail-type=ALL
#SBATCH --mail-user="danielloftus@g.harvard.edu"

module load bowtie2/2.3.2-fasrc02

#remap reads crossing a SNP using the same bowtie2 options as in the original mapping (see HiCUP_mapper line 332)

my_index=/n/whipple_lab/Users/dloftus/genome_references_indices/bowtie2_mm10/mm10_index

output_dir=/n/whipple_lab/Users/dloftus/capture_hic/results_run1/hicup/WASP/alignment_output/bam_to_fastq_for_wasp_single_end/single_end_realignment

for fastq in C1_S1_R1_2_001.hicup.read1.fq C1_S1_R1_2_001.hicup.read2.fq ; do

bowtie2 --very-sensitive -x $my_index --no-unal --threads 32 -U $fastq --reorder -S ${output_dir}/${fastq%.gz}.sam 

module load samtools/1.5-fasrc02

samtools view -bS ${output_dir}/${fastq%.gz}.sam > ${output_dir}/${fastq%.gz}.bam
samtools sort ${output_dir}/${fastq%.gz}.bam -o ${output_dir}/${fastq%.gz}.sorted.bam

done 