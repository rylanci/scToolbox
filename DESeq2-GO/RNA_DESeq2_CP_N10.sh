#!/bin/bash
#PBS -q hotel
#PBS -m ae
#PBS -M rlancione@ucsd.edu
#PBS -V
#PBS -A epigen-group
#PBS -e /home/rlancione/jeo/err/LungBPD/annoRNA_DEseq_enrcher_select_condi.R.e
#PBS -o /home/rlancione/jeo/out/LungBPD/annoRNA_DEseq_enrichr_select_condi.R.o
#PBS -l walltime=8:00:00
#PBS -l nodes=1:ppn=20

# Navigate to oasis scratch
cd /oasis/tscc/scratch/rlancione/
# activate cond env
source ~/miniconda3/etc/profile.d/conda.sh
conda activate renv4

# 1. Input merged object path
SEU_OBJ="/projects/ps-epigen/users/rlan/Lung_BPD/DESeq/healedControl_subD352/DE_Objects/Lung_BPD_hcontrol_subD352_batch.RDS"
#SEU_OBJ="/projects/ps-epigen/users/rlan/Lung_BPD/DESeq/healedControl_subD352/DE_Objects/Lung_BPD_cap1_utest_batch.RDS"
#SEU_OBJ="/projects/ps-epigen/users/rlan/Lung_BPD/DESeq/healedControl_subD352/DE_Objects/Lung_BPD_ciliated_utest_batch.RDS"

# Record Input
# 2. pathway to matrices and metadata created from MatrixCreation.R script
MATRIX_PATH="/home/rlancione/ps-epigen/users/rlan/Lung_BPD/DESeq/healedControl_subD352/bulk_matrices/"
# 3. adjusted p-value cutoff for EnrichR
PADJ_CUTOFF=0.1
# 4. log2FC threshold cutoff for EnrichR
# .58=log2(1.5)  1=log2(2)
LOG2FC_CUTOFF=0.5849625
# 5. path to output directory
OUTDIR='/projects/ps-epigen/users/rlan/Lung_BPD/DESeq/healedControl_subD352/DE_CP_UNV/Results_n10/'
# 6. project name (must be the same from matrix creation script)
PROJ='LungBPD_HC_subD352'
#7. Cell type column. Default seurat_clusters
CT_COL="celltype_manual_anno_broad"
#8. Path to tsv containing comparisons between conditions to run DESeq on. If left null (commed out) all conditions will be run pairwise.
COMP_PATH="/projects/ps-epigen/users/rlan/Lung_BPD/DESeq/healedControl_subD352/conditions.tsv"

# Path to script
Script='/projects/ps-epigen/users/rlan/Lung_BPD/DESeq/healedControl_subD352/DE_CP_UNV/RNA_pwDE_CP_unv.R'
######## Run R Script #########
Rscript $Script -sobj $SEU_OBJ -ctcol $CT_COL -mpath $MATRIX_PATH -ap $PADJ_CUTOFF -lfc $LOG2FC_CUTOFF \
	-o $OUTDIR -proj $PROJ -comp $COMP_PATH -form "~ donor_batch + condition" -nf 10
