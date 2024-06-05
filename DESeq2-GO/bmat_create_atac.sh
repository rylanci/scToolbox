#!/bin/sh
#SBATCH -J bulk_mat_create 
#SBATCH -N 1                           
#SBATCH -c 6                            
#SBATCH --mem=30G                       
#SBATCH -t 24:00:00
#SBATCH -o /tscc/nfs/home/rlancione/jeo/eye_macs3_bmatcreate.sh.o
#SBATCH -e /tscc/nfs/home/rlancione/jeo/eye_macs3_bmatcreate.sh.e
#SBATCH -p gold
#SBATCH -q hcg-csd772
#SBATCH -A csd772
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=rlancione@ucsd.edu
 
## Navigate to oasis scratch
#cd /tscc/lustre/scratch/rlancione/
# activate cond env
source ~/miniconda3/etc/profile.d/conda.sh
conda activate renv43


# Input merged object path 
SOBJ="/tscc/projects/ps-epigen/users/rlan/Liver_FNIH/Callpeaks/Objects/LFNIH_ct_edit.RDS"
# Path to peaks 
PP="/tscc/projects/ps-epigen/users/rlan/Liver_FNIH//Callpeaks/merge_out/"
# Dataset column name in object metadata. E.G. "orig.ident"
DSET_COL="Patient_identifier"
COND_COL="Final_Decision"
CT_COL="celltype_wnn_v03" 
#Output path
ODIR="/tscc/nfs/home/rlancione/myprojects/sandbox/ATAC_bmat_create/test_out/"
# Project: To name files
PROJ="LFNIH"
# n Cores for parrallel
NC=5

# Extra covariate columns for meta data. 
CV1="Gender"
CV2="exp_batch"


# Path to script
echo "run Rscript"
Script='/tscc/nfs/home/rlancione/scToolbox/DESeq2-GO/bmat_create_atac.R'
# args 1. input object path 2. output path 3. project name  
Rscript $Script -sobj $SOBJ -pp $PP -nc $NC -dcol $DSET_COL -ccol $COND_COL -ctcol $CT_COL --proj $PROJ -o $ODIR \
	-cov1 $CV1 -cov2 $CV2 


