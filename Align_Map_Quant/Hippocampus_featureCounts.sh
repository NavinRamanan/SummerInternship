#! /bin/bash -l
#$ -P adomics
#$ -l h_rt=24:00:00
#$ -pe omp 24
#$ -m e
#$ -N featureCounts_Hippocampus_new
#$ -j y
#$ -o featureCounts_Hippocampus_new.qlog


module load subread
featureCounts -p -a /project/adomics/shared/Thomas_STAR/Mouse_FASTA_Genome_and_GTF_Annotation_From_ENSEMBL/Mus_musculus.GRCm39.105.gtf -o /projectnb/adomics/navin_analysis/featureCounts/Hippocampus_gene_name_counts_matrix.txt /projectnb/adomics/navin_analysis/STAR/Hippocampus_STAR_output/projectnb/adomics/navin_analysis/STAR/Hippocampus/*Aligned.sortedByCoord.out.bam
