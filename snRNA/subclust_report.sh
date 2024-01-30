#!/bin/bash
#PBS -q hotel
#PBS -m ae
#PBS -M rlancione@ucsd.edu
#PBS -V
#PBS -A epigen-group
#PBS -e /home/rlancione/jeo/err/LungBPD/subclust_CAP1.2.R.e
#PBS -o /home/rlancione/jeo/out/LungBPD/subclust_CAP1.2.R.o
#PBS -l walltime=48:00:00
#PBS -l nodes=1:ppn=20

# Navigate to oasis scratch
cd /oasis/tscc/scratch/rlancione/
# activate cond env
source ~/miniconda3/etc/profile.d/conda.sh
conda activate renv4

### INPUTS 
# object path
SEU_OBJ="/projects/ps-epigen/users/rlan/Space_Organoid/Reports/Objects/so_merge_alldata_sobj.RDS"
# output path
OUT_PATH="/projects/ps-epigen/users/rlan/Space_Organoid/Reports/Output/"
# celltype
CELLTYPE="SO"
# assay to test
ASSAY="RNA"
# test to use 
TEST="MAST"
# Project name 
NAME="merge_alldata_mast"
# column to correct by
HARM_COL="orig.ident"
# Comparison table
#COMP_PATH="/projects/ps-epigen/users/rlan/Lung_BPD/Subclustering/comparisons.txt"

#### Processing & QC script
QC_Report='/projects/ps-epigen/users/rlan/Space_Organoid/Reports/subclust_proc_qc.R'
Rscript $QC_Report -sobj $SEU_OBJ -c $CELLTYPE -n $NAME -o $OUT_PATH \
	-hg $HARM_COL

#### FindMarker report script 
FM_Report='/projects/ps-epigen/users/rlan/Lung_BPD/Subclustering/scripts/subclust_fm.R'
Rscript $FM_Report -o $OUT_PATH -a $ASSAY -t $TEST -n $NAME -w 18

#### scCustomize dotplot script 
conda deactivate 
conda activate seurat
scDP="/projects/ps-epigen/users/rlan/Lung_BPD/Subclustering/scripts/clust_dotplot.R"
Rscript $scDP -o $OUT_PATH -n $NAME 


