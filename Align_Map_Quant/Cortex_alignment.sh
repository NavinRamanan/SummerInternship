#! /bin/bash -l
#$ -P adomics
#$ -l h_rt=96:00:00
#$ -pe omp 24
#$ -m e
#$ -N STAR_alignment_cortex
#$ -j y
#$ -o STAR_alignment_cortex.qlog
module load star
for i in /projectnb/adomics/navin_analysis/STAR/Cortex/*_R1.fastq.gz
do
STAR --runMode alignReads --runThreadN 24 --genomeDir /projectnb/adomics/rna_sequencing_data/STAR/star_mouse_genome_index --readFilesCommand zcat --readFilesIn $i ${i%_R1.fastq.gz}_R2.fastq.gz --outSAMtype BAM Unsorted SortedByCoordinate --outFileNamePrefix /projectnb/adomics/navin_analysis/STAR/Cortex_STAR_output/${i%_R1.fastq.gz} --quantMode GeneCounts TranscriptomeSAM
done 