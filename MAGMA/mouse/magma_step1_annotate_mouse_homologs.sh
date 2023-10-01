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

#Variables for specific paths. 
ANNO=/dcs04/lieber/marmaypag/ls_molecular-profiling_LIBD1070/ls_molecular-profiling/mouse/mouse_expressing_hg19_homs.gene.loc
SNPDATA=/dcs04/lieber/marmaypag/ls_molecular-profiling_LIBD1070/ls_molecular-profiling/MAGMA/SNP_Data
MOUSE=/dcs04/lieber/marmaypag/ls_molecular-profiling_LIBD1070/ls_molecular-profiling/MAGMA/SNP_Data/mouse

#Variables for .snploc files.
AD=/dcs04/lieber/marmaypag/ls_molecular-profiling_LIBD1070/ls_molecular-profiling/MAGMA/SNP_Data/GWAS_Tables/AD
ADHD=/dcs04/lieber/marmaypag/ls_molecular-profiling_LIBD1070/ls_molecular-profiling/MAGMA/SNP_Data/GWAS_Tables/ADHD
ASD=/dcs04/lieber/marmaypag/ls_molecular-profiling_LIBD1070/ls_molecular-profiling/MAGMA/SNP_Data/GWAS_Tables/ASD
BIP=/dcs04/lieber/marmaypag/ls_molecular-profiling_LIBD1070/ls_molecular-profiling/MAGMA/SNP_Data/GWAS_Tables/BIP
CUD=/dcs04/lieber/marmaypag/ls_molecular-profiling_LIBD1070/ls_molecular-profiling/MAGMA/SNP_Data/GWAS_Tables/CUD
OUD=/dcs04/lieber/marmaypag/ls_molecular-profiling_LIBD1070/ls_molecular-profiling/MAGMA/SNP_Data/GWAS_Tables/OUD
PD=/dcs04/lieber/marmaypag/ls_molecular-profiling_LIBD1070/ls_molecular-profiling/MAGMA/SNP_Data/GWAS_Tables/PD
SCZ=/dcs04/lieber/marmaypag/ls_molecular-profiling_LIBD1070/ls_molecular-profiling/MAGMA/SNP_Data/GWAS_Tables/SCZ

## Step 1 - Annotation (SNP : gene mapping)
#ADHD
magma --annotate window=35,10 --snp-loc $ADHD/adhd_jul2017.snploc --gene-loc $ANNO --out $MOUSE/adhd_jul2017_mouse_homologs

#Alzheimers
magma --annotate window=35,10 --snp-loc $AD/AD_sumstats_Jansenetal_2019sept.snploc --gene-loc $ANNO --out $MOUSE/AD_Jansenetal_2019sept_mouse_homologs

#ASD
magma --annotate window=35,10 --snp-loc $ASD/iPSYCH-PGC_ASD_Nov2017.snploc --gene-loc $ANNO --out $MOUSE/iPSYCH-PGC_ASD_Nov2017_mouse_homologs

#BIP
magma --annotate window=35,10 --snp-loc $BIP/PGC_BIP32b.snploc --gene-loc $ANNO --out $MOUSE/BIP_PGC_2018_mouse_homologs

#CUD
magma --annotate window=35,10 --snp-loc $CUD/CUD_GWAS_iPSYCH_June2019.snploc --gene-loc $ANNO --out $MOUSE/CUD_iPSYCH_June2019_mouse_homologs

#PD
magma --annotate window=35,10 --snp-loc $PD/nallsetal2019_PD.snploc --gene-loc $ANNO --out $MOUSE/PD_nallsetal2019_mouse_homologs

#SCZ
magma --annotate window=35,10 --snp-loc $SCZ/PGC_SCZ.snploc --gene-loc $ANNO --out $MOUSE/PGC_SCZ_mouse_homologs

#OpioidUseDisorder
magma --annotate window=35,10 --snp-loc $OUD/ODvsEx_AA_mapped.snploc --gene-loc $ANNO --out $MOUSE/ODvsEx_AA_mouse_homologs
magma --annotate window=35,10 --snp-loc $OUD/ODvsEx_EA_mapped.snploc --gene-loc $ANNO --out $MOUSE/ODvsEx_EA_mapped_homologs
magma --annotate window=35,10 --snp-loc $OUD/ODvsEx_mapped.snploc --gene-loc $ANNO --out $MOUSE/ODvsEx_mapped_homologs
magma --annotate window=35,10 --snp-loc $OUD/ODvsUn_AA_mapped.snploc --gene-loc $ANNO --out $MOUSE/ODvsUn_AA_mapped_homologs
magma --annotate window=35,10 --snp-loc $OUD/ODvsUn_EA_mapped.snploc --gene-loc $ANNO --out $MOUSE/ODvsUn_EA_mapped_homologs
magma --annotate window=35,10 --snp-loc $OUD/ODvsUn_mapped.snploc --gene-loc $ANNO --out $MOUSE/ODvsUn_mapped_homologs


echo "**** Job ends ****"
date
