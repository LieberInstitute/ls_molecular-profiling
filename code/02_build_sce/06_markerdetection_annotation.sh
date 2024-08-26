#!/bin/bash
#
#SBATCH --job-name=06_anno_DEG
#SBATCH --output=logs/06_anno_DEG.log
#SBATCH --error=logs/06_anno_DEG.log
#
# Number of CPUs allocated to each task.
#SBATCH --cpus-per-task=2
#
# Mimimum memory required per allocated  CPU
#SBATCH --mem-per-cpu=5G
#
# Send mail to the email address when the job fails
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=robert.phillips@libd.org


echo "********* Job Starts *********"
date

#load R
module load conda_R/4.4

#list modules for reproducibility purposes
module list

#run the Rjob
Rscript 06_markerdetection_annotation.R 

echo "********* Job Ends *********"
date
