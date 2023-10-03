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
ANNO_Hg19=/dcs04/lieber/marmaypag/ls_molecular-profiling_LIBD1070/ls_molecular-profiling/MAGMA/mouse_analysis/mouse_expressing_hg19_homs.gene.loc
ANNO_Hg38=/dcs04/lieber/marmaypag/ls_molecular-profiling_LIBD1070/ls_molecular-profiling/MAGMA/mouse_analysis/mouse_expressing_hg38_homs.gene.loc
MOUSE=/dcs04/lieber/marmaypag/ls_molecular-profiling_LIBD1070/ls_molecular-profiling/MAGMA/mouse_analysis/magma_output

#Variables for .snploc files.
AD=/dcs04/lieber/marmaypag/ls_molecular-profiling_LIBD1070/ls_molecular-profiling/MAGMA/SNP_Data/GWAS_Tables/AD
ADHD=/dcs04/lieber/marmaypag/ls_molecular-profiling_LIBD1070/ls_molecular-profiling/MAGMA/SNP_Data/GWAS_Tables/ADHD
ASD=/dcs04/lieber/marmaypag/ls_molecular-profiling_LIBD1070/ls_molecular-profiling/MAGMA/SNP_Data/GWAS_Tables/ASD
BIP=/dcs04/lieber/marmaypag/ls_molecular-profiling_LIBD1070/ls_molecular-profiling/MAGMA/SNP_Data/GWAS_Tables/BIP
CUD=/dcs04/lieber/marmaypag/ls_molecular-profiling_LIBD1070/ls_molecular-profiling/MAGMA/SNP_Data/GWAS_Tables/CUD
SCZ=/dcs04/lieber/marmaypag/ls_molecular-profiling_LIBD1070/ls_molecular-profiling/MAGMA/SNP_Data/GWAS_Tables/SCZ
OUD=/dcs04/lieber/marmaypag/ls_molecular-profiling_LIBD1070/ls_molecular-profiling/MAGMA/SNP_Data/GWAS_Tables/OUD
MDD=/dcs04/lieber/marmaypag/ls_molecular-profiling_LIBD1070/ls_molecular-profiling/MAGMA/SNP_Data/GWAS_Tables/MDD

## Step 1 - Annotation (SNP : gene mapping)
#Alzheimers
magma --annotate window=35,10 --snp-loc $AD/AD_sumstats_Jansenetal_2019sept.snploc --gene-loc $ANNO_Hg19 --out $MOUSE/AD_Jansenetal_2019sept_mouse_homologs

#ADHD
magma --annotate window=35,10 --snp-loc $ADHD/adhd_jul2017.snploc --gene-loc $ANNO_Hg19 --out $MOUSE/adhd_jul2017_mouse_homologs

#ASD
magma --annotate window=35,10 --snp-loc $ASD/iPSYCH-PGC_ASD_Nov2017.snploc --gene-loc $ANNO_Hg19 --out $MOUSE/iPSYCH-PGC_ASD_Nov2017_mouse_homologs

#BIP
magma --annotate window=35,10 --snp-loc $BIP/PGC_BIP32b.snploc --gene-loc $ANNO_Hg19 --out $MOUSE/BIP_PGC_2018_mouse_homologs

#CUD
magma --annotate window=35,10 --snp-loc $CUD/CUD_GWAS_iPSYCH_June2019.snploc --gene-loc $ANNO_Hg19 --out $MOUSE/CUD_iPSYCH_June2019_mouse_homologs

#SCZ
magma --annotate window=35,10 --snp-loc $SCZ/PGC3_SCZ_wave3.european.autosome.public.v3.snploc --gene-loc $ANNO_Hg19 --out $MOUSE/PGC_SCZ_mouse_homologs

#MDD 
magma --annotate window=35,10 --snp-loc $MDD/MDD.phs001672.pha005122.snploc --gene-loc $ANNO_Hg38 --out $MOUSE/MDD_pha005122_mouse_homologs

#OUD 
magma --annotate window=35,10 --snp-loc $OUD/OUD.phs001672.pha004954.snploc --gene-loc $ANNO_Hg38 --out $MOUSE/OUD_pha004954_mouse_homologs

echo "**** Job ends ****"
date
