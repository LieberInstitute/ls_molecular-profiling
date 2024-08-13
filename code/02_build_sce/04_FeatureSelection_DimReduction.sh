#!/bin/bash
#
#SBATCH --job-name=04_FeatSlxn_DimRed
#SBATCH --output=logs/04_FeatSlxn_DimRed.log
#SBATCH --error=logs/04_FeatSlxn_DimRed.log
#
# Number of CPUs allocated to each task.
#SBATCH --cpus-per-task=2
#
# Mimimum memory required per allocated  CPU
#SBATCH --mem-per-cpu=10G
#
# Send mail to the email address when the job fails
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=robert.phillips@libd.org


echo "********* Job Starts *********"
date

#load R
module load conda_R/4.3

#list modules for reproducibility purposes
module list

#run the Rjob
Rscript 04_FeatureSelection_DimReduction.R 

echo "********* Job Ends *********"
date
