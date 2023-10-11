#!/bin/bash
#SBATCH --mem=80G
#SBATCH --job-name=lscellranger
#SBATCH -o logs/cellranger-2310.txt
#SBATCH --array=1-3

# Files to process
# 5c_dACC_MRV 

echo "**** Job starts ****"
date

echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${SLURM_JOBID}"
echo "Job name: ${SLURM_JOB_NAME}"
echo "Hostname: ${SLURM_NODENAME}"
echo "Task id: ${SLURM_ARRAY_TASK_ID}"

## load CellRanger
module load cellranger/7.2.0

## List current modules for reproducibility
module list

## Locate file
SAMPLE=$(awk 'BEGIN {FS="\t"} {print $1}' samples.txt | awk "NR==${SLURM_ARRAY_TASK_ID}")
SAMPLEID=$(awk 'BEGIN {FS="\t"} {print $2}' samples.txt | awk "NR==${SLURM_ARRAY_TASK_ID}")
## SAMPLE=$(awk "NR==${SGE_TASK_ID}" 03_cellranger.txt)
echo "Processing sample ${SAMPLE}"
date

## Run CellRanger
cellranger count --id=${SAMPLE} \
    --transcriptome=/dcs04/lieber/lcolladotor/annotationFiles_LIBD001/10x/refdata-gex-GRCh38-2020-A \
    --fastqs=/dcs04/lieber/marmaypag/ls_molecular-profiling_LIBD1070/ls_molecular-profiling/raw-data/FASTQ/${SAMPLE} \
    --sample=${SAMPLEID} \
    --jobmode=local \
    --localcores=8 \
    --localmem=64

## Move output
echo "Moving data to new location"
date
mkdir -p /dcs04/lieber/marmaypag/ls_molecular-profiling_LIBD1070/ls_molecular-profiling/processed-data/01_cellranger/${SAMPLE}
mv ${SAMPLE} /dcs04/lieber/marmaypag/ls_molecular-profiling_LIBD1070/ls_molecular-profiling/processed-data/01_cellranger/

echo "**** Job ends ****"
date
