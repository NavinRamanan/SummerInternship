#!/bin/bash -l
#$ -P adomics
#$ -l h_rt=48:00:00
#$ -pe omp 24
source /etc/bashrc

module load fastqc

for i in /projectnb/adomics/rna_sequencing_data/Batch_of_192/FASTQ/Hippocampus/*.fastq.gz
do
  	fastqc -o /projectnb/adomics/navin_analysis/fastqc/hc_fastqc $i
done
