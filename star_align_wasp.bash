#!/bin/bash
#SBATCH -n 32 # Number of cores
#SBATCH -N 1 # Ensure that all cores are on one machine
#SBATCH -t 0-04:00 # Runtime in D-HH:MM
#SBATCH -p whipple # Partition to submit to
#SBATCH --mem=64G # Memory pool for all cores (see also --mem-per-cpu)
#SBATCH -o %x%j.out # File to which STDOUT will be written
#SBATCH -e %x%j.err # File to which STDERR will be written
#SBATCH --mail-type=ALL
#SBATCH --mail-user="danielloftus@g.harvard.edu"

module load parallel/20180522-fasrc01
module load GCC/7.3.0-2.30 OpenMPI/3.1.1 STAR/2.7.1a
module load samtools/1.5-fasrc02

echo MODULES LOADED 
date

genome=/n/whipple_lab/Users/dloftus/genome_references_indices/star_mm10_with_annotations_no_chr

#align reads using STAR 

vcfFile=/n/whipple_lab/Users/dloftus/vcf_files/CAST_EiJ.mgp.v5.snps.dbSNP142.vcf
output_dir=/n/whipple_lab/Users/dloftus/capture_hic/results_run1/hicup/WASP/alignment_output/bam_to_fastq_for_wasp_single_end/star_wasp_alignment

STAR \
--genomeDir $genome \
--runThreadN 32 \
--readFilesIn C1_S1_R1_2_001.hicup.all_reads.fq \
--outFileNamePrefix ${output_dir}/C1_S1_R1_2_001.hicup.all_reads.wasp \
--outSAMtype BAM SortedByCoordinate \
--waspOutputMode SAMtag \
--varVCFfile $vcfFile \
--outSAMattributes vW \
--alignIntronMax 1 

echo COMPLETED
date 


