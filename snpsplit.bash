#!/bin/bash
#SBATCH -n 1 # Number of cores
#SBATCH -N 1 # Ensure that all cores are on one machine
#SBATCH -t 0-15:00 # Runtime in D-HH:MM
#SBATCH -p whipple # Partition to submit to
#SBATCH --mem=32G # Memory pool for all cores (see also --mem-per-cpu)
#SBATCH -o %x%j.out # File to which STDOUT will be written
#SBATCH -e %x%j.err # File to which STDERR will be written
#SBATCH --mail-type=ALL
#SBATCH --mail-user="danielloftus@g.harvard.edu"


module load GCCcore/8.2.0 Perl/5.28.0
module load samtools/1.5-fasrc02

snp_file=/n/whipple_lab/Users/dloftus/genome_references_indices/mm10_N-masked_b6_cast_snpsplit/all_SNPs_CAST_EiJ_GRCm38.txt 
input_bam=/n/whipple_lab/Users/dloftus/capture_hic/results_run1/hicup/WASP/alignment_output/bam_to_fastq_for_wasp_single_end/single_end_realignment/find_intersecting_snps_output_PASS_ONLY/realignment_output/reads_to_keep_bams/n_masked_alignment/C1_S1_R1_2_001.hicup.wasp.n_masked.sorted_by_name.filtered.bam 
output_dir=/n/whipple_lab/Users/dloftus/capture_hic/results_run1/hicup/WASP/alignment_output/bam_to_fastq_for_wasp_single_end/single_end_realignment/find_intersecting_snps_output_PASS_ONLY/realignment_output/reads_to_keep_bams/n_masked_alignment/final_snpsplit

perl SNPsplit --hic -o $output_dir --snp_file $snp_file $input_bam

