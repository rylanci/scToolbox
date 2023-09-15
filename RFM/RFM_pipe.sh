#!/bin/bash
#PBS -q hotel
#PBS -m ae
#PBS -M rlancione@ucsd.edu
#PBS -V
#PBS -A epigen-group
#PBS -e /home/rlancione/jeo/err/rmf_test.R.e
#PBS -o /home/rlancione/jeo/out/rfm_test.R.o
#PBS -l walltime=8:00:00
#PBS -l nodes=1:ppn=10

# Navigate to oasis scratch
cd /oasis/tscc/scratch/rlancione/
# activate cond env
source ~/miniconda3/etc/profile.d/conda.sh
conda activate renv4

# 1. Input merged object path
SEU_OBJ="/projects/ps-epigen/users/rlan/Lung_BPD/DESeq/healedControl_subD352/DE_Objects/Lung_BPD_hcontrol_subD352_batch.RDS"

# 5. path to output directory
OUTDIR='/projects/ps-epigen/users/rlan/Lung_BPD/DESeq/healedControl_subD352/DE_CP_UNV/Results_n10/'
# 6. project name (must be the same from matrix creation script)
PROJ='LungBPD_HC_subD352'
#7. Cell type column. Default seurat_clusters
CT_COL="celltype_manual_anno_broad"
#8. Path to tsv containing comparisons between conditions to run DESeq on. If left null (commed out) all conditions will be run pairwise.
COMP_PATH="/projects/ps-epigen/users/rlan/Lung_BPD/DESeq/healedControl_subD352/conditions.tsv"

# Path to script
Script='/home/rlancione/scToolbox/RFM/RFM_pipe.R'
######## Run R Script #########
Rscript $Script -sobj $SEU_OBJ -ctcol $CT_COL -mpath $MATRIX_PATH -ap $PADJ_CUTOFF -lfc $LOG2FC_CUTOFF \
	-o $OUTDIR -proj $PROJ -comp $COMP_PATH -form "~ donor_batch + condition" -nf 10
