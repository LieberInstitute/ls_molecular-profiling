#!/bin/bash
#
#SBATCH --job-name=Fig2_Plots
#SBATCH --output=logs/Fig2_Plots.log
#SBATCH --error=logs/Fig2_Plots.log
#
# Number of CPUs allocated to each task.
#SBATCH --cpus-per-task=1
#
# Mimimum memory required per allocated  CPU
#SBATCH --mem-per-cpu=12G
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
Rscript Figure2_Plots.R

echo "********* Job Ends *********"
date
