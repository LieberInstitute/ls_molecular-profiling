#!/bin/bash
#
#SBATCH --job-name=emptyDrops
#SBATCH --output=emptyDrops_dropletQC.out
#SBATCH --error=emptyDrops_dropletQC.err
#
# Number of CPUs allocated to each task.
#SBATCH --cpus-per-task=5
#
# Mimimum memory required per allocated  CPU
#SBATCH --mem-per-cpu=20G
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
Rscript 02_calculate_droplet_scores.R

echo "********* Job Ends *********"
date
