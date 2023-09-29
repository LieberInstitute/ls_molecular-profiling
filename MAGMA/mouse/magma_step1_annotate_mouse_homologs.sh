#!/bin/bash
#
#SBATCH --job-name=magma_mouse_homs
#SBATCH --output=magma_mouse_homs.out
#SBATCH --error=magma_mouse_homs.err
#
# Number of CPUs allocated to each task.
#SBATCH --cpus-per-task=5
#
# Mimimum memory required per allocated  CPU  in  MegaBytes.
#SBATCH --mem-per-cpu=5G
#
# Send mail to the email address when the job fails
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=robert.phillips@libd.org


echo "**** Job starts ****"
date

#Load the magma module
module load magma/1.10


#Create
ANNO=/dcs04/lieber/marmaypag/ls_molecular-profiling_LIBD1070/ls_molecular-profiling/mouse/mouse_expressing_hg19_homs.gene.loc
SNPDATA=/dcs04/lieber/marmaypag/ls_molecular-profiling_LIBD1070/ls_molecular-profiling/MAGMA/SNP_Data

## Step 1 - Annotation (SNP : gene mapping)
#ADHD
magma --annotate window=35,10 --snp-loc $SNPDATA/ADHD_PGC2018.snploc --gene-loc $ANNO --out $SNP_Data/mouse/ADHD_PGC2018_mouse_homologs

#Alzheimers
magma --annotate window=35,10 --snp-loc $SNPDATA/AD_PGC-IGAP-ADSP-UKB_2019.snploc --gene-loc $ANNO --out $SNP_Data/mouse/AD_PGC-IGAP-ADSP-UKB_2019_mouse_homologs

#ASD
magma --annotate window=35,10 --snp-loc $SNPDATA/iPSYCH-PGC_ASD_Nov2017.snploc --gene-loc $ANNO --out $SNP_Data/mouse/iPSYCH-PGC_ASD_Nov2017_mouse_homologs

#PTSD
magma --annotate window=35,10 --snp-loc $SNPDATA/PTSD_Nievergelt2019.snploc --gene-loc $ANNO --out $SNP_Data/mouse/PTSD_Nievergelt2019_mouse_homologs

#Alcohol/Tobacco
magma --annotate window=35,10 --snp-loc $SNPDATA/AgeofInitiation.snploc --gene-loc $ANNO --out $SNP_Data/mouse/AgeofInitiation_mouse_homologs
magma --annotate window=35,10 --snp-loc $SNPDATA/CigarettesPerDay.snploc --gene-loc $ANNO --out $SNP_Data/mouse/CigarettesPerDay_mouse_homologs
magma --annotate window=35,10 --snp-loc $SNPDATA/DrinksPerWeek.snploc --gene-loc $ANNO --out $SNP_Data/mouse/DrinksPerWeek_mouse_homologs
magma --annotate window=35,10 --snp-loc $SNPDATA/SmokingCessation.snploc --gene-loc $ANNO --out $SNP_Data/mouse/SmokingCessation_mouse_homologs

#OpioidUseDisorder
magma --annotate window=35,10 --snp-loc $SNPDATA/ODvsEx_AA_mapped.snploc --gene-loc $ANNO --out $SNP_Data/mouse/ODvsEx_AA_mouse_homologs
magma --annotate window=35,10 --snp-loc $SNPDATA/ODvsEx_EA_mapped.snploc --gene-loc $ANNO --out $SNP_Data/mouse/ODvsEx_EA_mapped_homologs
magma --annotate window=35,10 --snp-loc $SNPDATA/ODvsEx_mapped.snploc --gene-loc $ANNO --out $SNP_Data/mouse/ODvsEx_mapped_homolocgs
magma --annotate window=35,10 --snp-loc $SNPDATA/ODvsUn_AA_mapped.snploc --gene-loc $ANNO --out $SNP_Data/mouse/ODvsUn_AA_mapped_homologs
magma --annotate window=35,10 --snp-loc $SNPDATA/ODvsUn_EA_mapped.snploc --gene-loc $ANNO --out $SNP_Data/mouse/ODvsUn_EA_mapped_homologs
magma --annotate window=35,10 --snp-loc $SNPDATA/ODvsUn_mapped.snploc --gene-loc $ANNO --out $SNP_Data/mouse/ODvsUn_mapped_homologs

#Coronary Artery Disease - negative control. 
magma --annotate window=35,10 --snp-loc $SNPDATA/CoronaryArteryDisease.snploc --gene-loc $ANNO --out $SNP_Data/mouse/CoronaryArteryDisease_mouse_homologs

echo "**** Job ends ****"
date