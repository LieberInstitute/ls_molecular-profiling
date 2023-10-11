#!/bin/bash
#SBATCH --array=1-3

SAMPLE=$(awk "NR==${SLURM_ARRAY_TASK_ID}" raw-mkdir.txt)

mkdir -p ./logs/
mkdir -p /dcs04/lieber/marmaypag/ls_molecular-profiling_LIBD1070/ls_molecular-profiling/raw-data/FASTQ/${SAMPLE}/

mv slurm-*.out ./logs