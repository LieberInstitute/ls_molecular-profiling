#!/bin/bash
#SBATCH --job-name="Human stat shuffling"
#SBATCH --mem=10GB
#SBATCH --output=human_shuffle.out
#SBATCH --error=human_shuffle.err

echo "`date`: loading modules"
module use /jhpce/shared/libd/modulefiles
module load conda_R/4.3

cd /dcs04/lieber/marmaypag/ls_molecular-profiling_LIBD1070/ls_molecular-profiling/

echo "`date`: Running R Job"
R CMD BATCH shuffle_human_t_stats.R

echo "`date`: R job complete"


 
 