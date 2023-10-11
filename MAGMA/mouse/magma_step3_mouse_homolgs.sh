#!/bin/bash
#
#SBATCH --job-name=magma_step3_mouse
#SBATCH --output=magma_step3_mouse.out
#SBATCH --error=magma_step3_mouse.err
#
# Number of CPUs allocated to each task.
#SBATCH --cpus-per-task=1
#
# Mimimum memory required per allocated  CPU  in  MegaBytes.
#SBATCH --mem-per-cpu=10G
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
MOUSE=/dcs04/lieber/marmaypag/ls_molecular-profiling_LIBD1070/ls_molecular-profiling/MAGMA/mouse_analysis/mouse_markers_human_homologs_1e-2.txt
MAGMA_OUTPUT=/dcs04/lieber/marmaypag/ls_molecular-profiling_LIBD1070/ls_molecular-profiling/MAGMA/mouse_analysis/magma_output

## Step 3 
#AD
magma --gene-results $MAGMA_OUTPUT/AD.genes.raw --set-annot $MOUSE gene-col=${genecol} set-col=${setcol} --out $MAGMA_OUTPUT/AD

#ADHD
magma --gene-results $MAGMA_OUTPUT/adhd.genes.raw --set-annot $MOUSE gene-col=${genecol} set-col=${setcol} --out $MAGMA_OUTPUT/ADHD

#ASD
magma --gene-results $MAGMA_OUTPUT/ASD.genes.raw --set-annot $MOUSE gene-col=${genecol} set-col=${setcol} --out $MAGMA_OUTPUT/ASD

#BIP
magma --gene-results $MAGMA_OUTPUT/BIP.genes.raw --set-annot $MOUSE gene-col=${genecol} set-col=${setcol} --out $MAGMA_OUTPUT/BIP

#CUD
magma --gene-results $MAGMA_OUTPUT/CUD.genes.raw --set-annot $MOUSE gene-col=${genecol} set-col=${setcol} --out $MAGMA_OUTPUT/CUD

#MDD
magma --gene-results $MAGMA_OUTPUT/MDD.genes.raw --set-annot $MOUSE gene-col=${genecol} set-col=${setcol} --out $MAGMA_OUTPUT/MDD

#OUD
magma --gene-results $MAGMA_OUTPUT/OUD.genes.raw --set-annot $MOUSE gene-col=${genecol} set-col=${setcol} --out $MAGMA_OUTPUT/OUD

#SCZ
magma --gene-results $MAGMA_OUTPUT/SCZ.genes.raw --set-annot $MOUSE gene-col=${genecol} set-col=${setcol} --out $MAGMA_OUTPUT/SCZ


echo "**** Job ends ****"
date
