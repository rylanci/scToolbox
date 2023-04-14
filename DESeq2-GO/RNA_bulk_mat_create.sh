#!/bin/sh
#PBS -q pdafm
#PBS -m abe
#PBS -M rlancione@ucsd.edu
#PBS -V
#PBS -A epigen-group
#PBS -e /home/rlancione/jeo/err/LungBPD/annoRNA_MC.R.e
#PBS -o /home/rlancione/jeo/out/LungBPD/annoRNA_MC.R.o
#PBS -l walltime=12:00:00
#PBS -l nodes=1:ppn=4


# Navigate to oasis scratch
cd /oasis/tscc/scratch/rlancione/
# activate cond env
source ~/miniconda3/etc/profile.d/conda.sh
conda activate renv4


# Input merged object path 
IDIR="/home/rlancione/ps-epigen/users/rlan/Lung_BPD/DESeq/healedControl_subD352/DE_Objects/Lung_BPD_hcontrol_subD352_batch.RDS"
# Dataset column name in object metadata. E.G. "orig.ident"
DSET_COL="patient"
COND_COL="condition"
CT_COL="celltype_manual_anno_broad" 
#Output path
ODIR="/home/rlancione/ps-epigen/users/rlan/Lung_BPD/DESeq/healedControl_subD352/bulk_matrices/"
# Project: To name files
PROJ="LungBPD_HC_subD352"
# Extra covariate columns for meta data. 
CV1="sex"
CV2="donor_batch"
# Path to script
echo "run Rscript"
Script='/home/rlancione/ps-epigen/users/rlan/Lung_BPD/DESeq/healedControl_subD352/RNA_mCreate.R'
# args 1. input object path 2. output path 3. project name  
Rscript $Script -sobj $IDIR -dcol $DSET_COL -ccol $COND_COL -ctcol $CT_COL --project_name $PROJ -o $ODIR \
	-cov1 $CV1 -cov2 $CV2 


