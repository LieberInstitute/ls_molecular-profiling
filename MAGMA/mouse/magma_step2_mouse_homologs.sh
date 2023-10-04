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
MOUSE=/dcs04/lieber/marmaypag/ls_molecular-profiling_LIBD1070/ls_molecular-profiling/MAGMA/mouse_analysis/magma_output
model="snp-wise"

#Variables for disorder/disease specific files containing all SNPs. 
AD=/dcs04/lieber/marmaypag/ls_molecular-profiling_LIBD1070/ls_molecular-profiling/MAGMA/SNP_Data/GWAS_Tables/AD
ADHD=/dcs04/lieber/marmaypag/ls_molecular-profiling_LIBD1070/ls_molecular-profiling/MAGMA/SNP_Data/GWAS_Tables/ADHD
ASD=/dcs04/lieber/marmaypag/ls_molecular-profiling_LIBD1070/ls_molecular-profiling/MAGMA/SNP_Data/GWAS_Tables/ASD
BIP=/dcs04/lieber/marmaypag/ls_molecular-profiling_LIBD1070/ls_molecular-profiling/MAGMA/SNP_Data/GWAS_Tables/BIP
CUD=/dcs04/lieber/marmaypag/ls_molecular-profiling_LIBD1070/ls_molecular-profiling/MAGMA/SNP_Data/GWAS_Tables/CUD
SCZ=/dcs04/lieber/marmaypag/ls_molecular-profiling_LIBD1070/ls_molecular-profiling/MAGMA/SNP_Data/GWAS_Tables/SCZ
OUD=/dcs04/lieber/marmaypag/ls_molecular-profiling_LIBD1070/ls_molecular-profiling/MAGMA/SNP_Data/GWAS_Tables/OUD
MDD=/dcs04/lieber/marmaypag/ls_molecular-profiling_LIBD1070/ls_molecular-profiling/MAGMA/SNP_Data/GWAS_Tables/MDD

## Step 2 - Gene analysis (from SNP-wise summary stats) 
#AD
magma --bfile $BFILE --gene-annot $MOUSE/AD_Jansenetal_2019sept_mouse_homologs.genes.annot --pval $AD/AD_sumstats_Jansenetal_2019sept.txt use=SNP,P ncol=Nsum --gene-model ${model} --out $AD/AD

#ADHD
magma --bfile $BFILE --gene-annot $MOUSE/adhd_jul2017_mouse_homologs.genes.annot --pval $ADHD/adhd_jul2017.txt use=SNP,P N=55374 --gene-model ${model} --out $ADHD/adhd

#ASD
magma --bfile $BFILE --gene-annot $MOUSE/iPSYCH-PGC_ASD_Nov2017_mouse_homologs.genes.annot --pval $ASD/iPSYCH-PGC_ASD_Nov2017.txt use=SNP,P N=46350 --gene-model ${model} --out $ASD/ASD

#BIP
magma --bfile $BFILE --gene-annot $MOUSE/BIP_PGC_2018_mouse_homologs.genes.annot --pval $BIP/PGC_BIP32b_withN.txt use=SNP,P ncol=Nca_Plus_Nco --gene-model ${model} --out $BIP/BIP

#CUD
magma --bfile $BFILE --gene-annot $MOUSE/CUD_iPSYCH_June2019_mouse_homologs.genes.annot --pval $CUD/CUD_GWAS_iPSYCH_June2019.txt use=SNP,P N=51372 --gene-model ${model} --out $CUD/CUD

#SCZ
magma --bfile $BFILE --gene-annot $MOUSE/PGC_SCZ_mouse_homologs.genes.annot --pval $SCZ/PGC3_SCZ_wave3.european.autosome.public.v3.pval ncol=N --gene-model ${model} --out $SCZ/SCZ

#MDD
magma --bfile $BFILE --gene-annot $MOUSE/MDD_pha005122_mouse_homologs.genes.annot --pval $MDD/MDD.phs001672.pha005122.pval N=1154267 --gene-model ${model} --out $MDD/MDD

#OUD
magma --bfile $BFILE --gene-annot $MOUSE/OUD_pha004954_mouse_homologs.genes.annot --pval $OUD/OUD.phs001672.pha004954.txt N=639063 --gene-model ${model} --out $OUD/OUD

echo "**** Job ends ****"
date
