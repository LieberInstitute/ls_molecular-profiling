#!/bin/bash
#
#SBATCH --job-name=RegionSpec_FeaturePlots
#SBATCH --output=logs/RegionSpec_FeaturePlots.log
#SBATCH --error=logs/RegionSpec_FeaturePlots.log
#
# Number of CPUs allocated to each task.
#SBATCH --cpus-per-task=1
#
# Mimimum memory required per allocated  CPU
#SBATCH --mem-per-cpu=4G
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
Rscript RegionSpecific_FeaturePlots.R

echo "********* Job Ends *********"
date
