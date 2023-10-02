#!/bin/bash
#
#SBATCH --job-name=magma_step3_mouse
#SBATCH --output=magma_step3_mouse.out
#SBATCH --error=magma_step3_mouse.err
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
setcol=1
genecol=2
model="snp-wise"
BFILE=/dcs04/lieber/marmaypag/ls_molecular-profiling_LIBD1070/ls_molecular-profiling/MAGMA/g1000_eur/g1000_eur
ANNO=/dcs04/lieber/marmaypag/ls_molecular-profiling_LIBD1070/ls_molecular-profiling/mouse/mouse_expressing_hg19_homs.gene.loc
SNPDATA=/dcs04/lieber/marmaypag/ls_molecular-profiling_LIBD1070/ls_molecular-profiling/MAGMA/SNP_Data/mouse
MOUSE=/dcs04/lieber/marmaypag/ls_molecular-profiling_LIBD1070/ls_molecular-profiling/mouse/mouse_markers_human_homologs_1e-2.txt

#Variables for genes.raw
AD=/dcs04/lieber/marmaypag/ls_molecular-profiling_LIBD1070/ls_molecular-profiling/MAGMA/SNP_Data/GWAS_Tables/AD
ADHD=/dcs04/lieber/marmaypag/ls_molecular-profiling_LIBD1070/ls_molecular-profiling/MAGMA/SNP_Data/GWAS_Tables/ADHD
ASD=/dcs04/lieber/marmaypag/ls_molecular-profiling_LIBD1070/ls_molecular-profiling/MAGMA/SNP_Data/GWAS_Tables/ASD
BIP=/dcs04/lieber/marmaypag/ls_molecular-profiling_LIBD1070/ls_molecular-profiling/MAGMA/SNP_Data/GWAS_Tables/BIP
CUD=/dcs04/lieber/marmaypag/ls_molecular-profiling_LIBD1070/ls_molecular-profiling/MAGMA/SNP_Data/GWAS_Tables/CUD
OUD=/dcs04/lieber/marmaypag/ls_molecular-profiling_LIBD1070/ls_molecular-profiling/MAGMA/SNP_Data/GWAS_Tables/OUD
PD=/dcs04/lieber/marmaypag/ls_molecular-profiling_LIBD1070/ls_molecular-profiling/MAGMA/SNP_Data/GWAS_Tables/PD
SCZ=/dcs04/lieber/marmaypag/ls_molecular-profiling_LIBD1070/ls_molecular-profiling/MAGMA/SNP_Data/GWAS_Tables/SCZ

## Step 3 
#AD
magma --gene-results $AD/AD.genes.raw --set-annot $MOUSE gene-col=${genecol} set-col=${setcol} --out $AD/AD

#ADHD
magma --gene-results $ADHD/adhd.genes.raw --set-annot $MOUSE gene-col=${genecol} set-col=${setcol} --out $ADHD/ADHD

#ASD
magma --gene-results $ASD/ASD.genes.raw --set-annot $MOUSE gene-col=${genecol} set-col=${setcol} --out $ASD/ASD

#BIP
magma --gene-results $BIP/BIP.genes.raw --set-annot $MOUSE gene-col=${genecol} set-col=${setcol} --out $BIP/BIP

#CUD
magma --gene-results $CUD/CUD.genes.raw --set-annot $MOUSE gene-col=${genecol} set-col=${setcol} --out $CUD/CUD

#OUD
magma --gene-results $OUD/ODvsEX_AA.genes.raw --set-annot $MOUSE gene-col=${genecol} set-col=${setcol} --out $OUD/ODvsEx_AA
magma --gene-results $OUD/ODvsEX_EA.genes.raw --set-annot $MOUSE gene-col=${genecol} set-col=${setcol} --out $OUD/ODvsEx_EA
magma --gene-results $OUD/ODvsEX_Trans.genes.raw --set-annot $MOUSE gene-col=${genecol} set-col=${setcol} --out $OUD/ODvsEx_Trans

#PD
magma --gene-results $PD/PD.genes.raw --set-annot $MOUSE gene-col=${genecol} set-col=${setcol} --out $PD/PD

#SCZ
magma --gene-results $SCZ/SCZ.genes.raw --set-annot $MOUSE gene-col=${genecol} set-col=${setcol} --out $SCZ/SCZ


echo "**** Job ends ****"
date
