#!/bin/bash
#
#SBATCH --job-name=human_LS_broad
#SBATCH --output=logs/human_LS_broad.log
#SBATCH --error=logs/human_LS_broad.log
#
# Number of CPUs allocated to each task.
#SBATCH --cpus-per-task=2
#
# Mimimum memory required per allocated  CPU
#SBATCH --mem-per-cpu=10G
#
# Send mail to the email address when the job fails
#SBATCH --mail-type=END
#SBATCH --mail-user=robert.phillips@libd.org


echo "********* Job Starts *********"
date

#load R
module load conda_R/4.4

#list modules for reproducibility purposes
module list

#run the Rjob
Rscript human_LS_broad_DEGs.R

echo "********* Job Ends *********"
date
