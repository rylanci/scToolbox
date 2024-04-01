#!/bin/bash
#SBATCH -J LungBPD_DESeq 
#SBATCH -N 1                            #Nodes
#SBATCH -c 3                            #Cores
#SBATCH --mem=20G                       #Memory
#SBATCH -t 24:00:00
#SBATCH -o /tscc/nfs/home/rlancione/jeo/lungbpd_deseq_ac.sh.o
#SBATCH -e /tscc/nfs/home/rlancione/jeo/lungbpd_deseq_ac.sh.e
#SBATCH -p platinum
#SBATCH -q hcp-csd772
#SBATCH -A csd772
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=rlancione@ucsd.edu

# Navigate to oasis scratch
cd /tscc/lustre/scratch/rlancione/
# activate cond env
source ~/miniconda3/etc/profile.d/conda.sh
conda activate renv43

# 1. Input merged object path
SEU_OBJ="/tscc/projects/ps-epigen/users/rlan/Lung_BPD/DESeq/healedControl_subD352/Iter12/Object/231010_02_lung_bpd_clust_anno_mitorm.RDS"

# Record Input
# 2. pathway to matrices and metadata created from MatrixCreation.R script
MATRIX_PATH="/tscc/projects/ps-epigen/users/rlan/Lung_BPD/DESeq/healedControl_subD352/Iter12/matrices/"
# 3. adjusted p-value cutoff for EnrichR
PADJ_CUTOFF=0.05
# 4. log2FC threshold cutoff for EnrichR
# .58=log2(1.5)  1=log2(2)
LOG2FC_CUTOFF=0.5849625
# 5. path to output directory
OUTDIR='/tscc/projects/ps-epigen/users/rlan/Lung_BPD/DESeq/healedControl_subD352/Iter12/Results_n10_p05_AC/'
# 6. project name (must be the same from matrix creation script)
PROJ='LungBPD'
#7. Cell type column. Default seurat_clusters
CT_COL="celltype_manual_anno_broad_v02"
#8. Path to tsv containing comparisons between conditions to run DESeq on. If left null (commed out) all conditions will be run pairwise.
COMP_PATH="/tscc/projects/ps-epigen/users/rlan/Lung_BPD/DESeq/healedControl_subD352/Iter12/comparisons_AC.tsv"


# Path to script
#Script='/tscc/projects/ps-epigen/users/rlan/Lung_BPD/DESeq/healedControl_subD352/Iter12/RNA_DEseq2_i12.R'
######## Run R Script #########
#Rscript $Script -sobj $SEU_OBJ -ctcol $CT_COL -mpath $MATRIX_PATH -ap $PADJ_CUTOFF -lfc $LOG2FC_CUTOFF \
#        -o $OUTDIR -proj $PROJ -comp $COMP_PATH -nf 10


conda deactivate 
conda activate CProfiler
#### Now call report script 
Report='/tscc/projects/ps-epigen/users/rlan/Lung_BPD/DESeq/healedControl_subD352/Iter12/RNA_DE_report_12.R'
######## Run R Script #########
Rscript $Report -rp $OUTDIR -mp $MATRIX_PATH -ap $PADJ_CUTOFF -lfc $LOG2FC_CUTOFF -org "human" \
	-pn $PROJ -sp $SEU_OBJ -ct $CT_COL



