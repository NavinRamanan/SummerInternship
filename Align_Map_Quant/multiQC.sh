#! /bin/bash -l
#$ -P adomics
#$ -l h_rt=48:00:00
#$ -pe omp 24
#$ -m e
#$ -N multiQC_Hippocampus
#$ -j y
#$ -o multiQC_Hippocampus.qlog
#$ -cwd


module load python3/3.7.9 
module load multiqc/1.10.1

multiqc /projectnb/adomics/navin_analysis/fastqc/hc_fastqc /projectnb/adomics/navin_analysis/STAR/Cortex_STAR_output/projectnb/adomics/navin_analysis/STAR/Cortex /projectnb/adomics/navin_analysis/featureCounts