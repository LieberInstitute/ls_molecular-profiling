#!/bin/bash
#
#SBATCH --job-name=magma_step2_mouse
#SBATCH --output=magma_step2_mouse.out
#SBATCH --error=magma_step2_mouse.err
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

#Variables for paths and such
BFILE=/dcs04/lieber/marmaypag/ls_molecular-profiling_LIBD1070/ls_molecular-profiling/MAGMA/g1000_eur/g1000_eur
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

## Step 2 - Gene analysis (from SNP-wise summary stats)
#AD
magma --bfile $BFILE --gene-annot $MOUSE/AD_Jansenetal_2019sept_mouse_homologs.genes.annot --pval $AD/AD_sumstats_Jansenetal_2019sept.txt use=SNP,P ncol=Nsum --out $AD/AD

#ADHD
magma --bfile $BFILE --gene-annot $MOUSE/adhd_jul2017_mouse_homologs.genes.annot --pval $ADHD/adhd_jul2017.txt use=SNP,P N=55374 --out $ADHD/adhd

#ASD
magma --bfile $BFILE --gene-annot $MOUSE/iPSYCH-PGC_ASD_Nov2017_mouse_homologs.genes.annot --pval $ASD/iPSYCH-PGC_ASD_Nov2017.txt use=SNP,P N=46351 --out $ASD/ASD

#BIP
magma --bfile $BFILE --gene-annot $MOUSE/BIP_PGC_2018_mouse_homologs.genes.annot --pval $BIP/PGC_BIP32b_withN.txt use=SNP,P ncol=Nca_Plus_Nco --out $BIP/BIP

#CUD
magma --bfile $BFILE --gene-annot $MOUSE/CUD_iPSYCH_June2019_mouse_homologs.genes.annot --pval $CUD/CUD_GWAS_iPSYCH_June2019.txt use=SNP,P N=51372 --out $CUD/CUD

#OUD
magma --bfile $BFILE --gene-annot $MOUSE/ODvsEx_AA_mouse_homologs.genes.annot --pval $OUD/OD_cases_vs._opioid-exposed_controls_in_African-ancestry_cohorts.txt use=rsID,P-Value ncol=Total_N  --out $OUD/ODvsEX_AA
magma --bfile $BFILE --gene-annot $MOUSE/ODvsEx_EA_mapped_homologs.genes.annot --pval $OUD/OD_cases_vs._opioid-exposed_controls_in_European-ancestry_cohorts.txt use=rsID,P-Value ncol=Total_N  --out $OUD/ODvsEX_EA
magma --bfile $BFILE --gene-annot $MOUSE/ODvsEx_mapped_homologs.genes.annot --pval $OUD/OD_cases_vs._opioid-exposed_controls_in_the_trans-ancestry_meta-analysis.txt use=rsID,P-Value ncol=Total_N  --out $OUD/ODvsEX_Trans
magma --bfile $BFILE --gene-annot $MOUSE/ODvsUn_AA_mapped_homologs.genes.annot --pval $OUD/OD_cases_vs._opioid-unexposed_controls_in_African-ancestry_cohorts.txt  use=rsID,P-Value ncol=Total_N  --out $OUD/ODvsUN_AA
magma --bfile $BFILE --gene-annot $MOUSE/ODvsUn_EA_mapped_homologs.genes.annot --pval $OUD/OD_cases_vs._opioid-unexposed_controls_in_European-ancestry_cohorts.txt   use=rsID,P-Value ncol=Total_N  --out $OUD/ODvsUN_EA
magma --bfile $BFILE --gene-annot $MOUSE/ODvsUn_mapped_homologs.genes.annot --pval $OUD/OD_cases_vs._opioid-unexposed_controls_in_the_trans-ancestry_meta-analysis.txt  use=rsID,P-Value ncol=Total_N  --out $OUD/ODvsUN_Trans

#PD
magma --bfile $BFILE --gene-annot $MOUSE/PD_nallsetal2019_mouse_homologs.genes.annot --pval $PD/nallsetal2019_mapped.txt use=ID,p ncol=N --out $PD/PD

#SCZ
magma --bfile $BFILE --gene-annot $MOUSE/PGC_SCZ_mouse_homologs.genes.annot --pval $SCZ/rall.txt use=snpid,p N=150064 --out $SCZ/SCZ

echo "**** Job ends ****"
date
