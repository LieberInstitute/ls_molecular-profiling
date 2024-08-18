#!/bin/bash
#
#SBATCH --job-name=05_clustering
#SBATCH --output=logs/05_clustering.log
#SBATCH --error=logs/05_clustering.log
#
# Number of CPUs allocated to each task.
#SBATCH --cpus-per-task=1
#
# Mimimum memory required per allocated  CPU
#SBATCH --mem-per-cpu=8G
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
Rscript 05_clustering.R 

echo "********* Job Ends *********"
date
