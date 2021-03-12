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

module load samtools/1.5-fasrc02

input_bam=/n/whipple_lab/Users/dloftus/capture_hic/results_run1/hicup/WASP/alignment_output/bam_to_fastq_for_wasp_single_end/single_end_realignment/find_intersecting_snps_output_PASS_ONLY/realignment_output/reads_to_keep_bams/n_masked_alignment/C1_S1_R1_2_001.hicup.wasp.n_masked.sorted.bam 
output_dir=/n/whipple_lab/Users/dloftus/capture_hic/results_run1/hicup/WASP/alignment_output/bam_to_fastq_for_wasp_single_end/single_end_realignment/find_intersecting_snps_output_PASS_ONLY/realignment_output/reads_to_keep_bams/n_masked_alignment

#samtools sort -n --threads 32 $input_bam -o ${output_dir}/C1_S1_R1_2_001.hicup.wasp.n_masked.sorted_by_name.bam



#samtools view C1_S1_R1_2_001.hicup.wasp.n_masked.sorted_by_name.bam | cut -f1 | uniq -c > counts_per_qname_space_delim.txt
#awk 'BEGIN{OFS="\t"} { if($1==1) print $2}' counts_per_qname_space_delim.txt > qnames_to_remove.txt



module load Java/1.8

input_bam_picard=C1_S1_R1_2_001.hicup.wasp.n_masked.sorted.bam

java -jar /n/whipple_lab/Users/dloftus/software/picard/picard.jar FilterSamReads \
--FILTER excludeReadList \
--READ_LIST_FILE qnames_to_remove.txt \
--SORT_ORDER queryname \
--INPUT $input_bam_picard \
--OUTPUT ${input_bam_picard%.bam}.filtered.bam 