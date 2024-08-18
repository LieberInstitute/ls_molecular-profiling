#!/bin/bash
#
#SBATCH --job-name=01_build_sce
#SBATCH --output=logs/01_build_sce.log
#SBATCH --error=logs/01_build_sce.log
#
# Number of CPUs allocated to each task.
#SBATCH --cpus-per-task=1
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
Rscript 01_build_basic_sce.R

echo "********* Job Ends *********"
date
