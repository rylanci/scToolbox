#!/bin/bash
#SBATCH -J bmat_create
#SBATCH -N 1                            #Nodes
#SBATCH -c 2                           #Cores
#SBATCH --mem=30G                       #Memory
#SBATCH -t 08:00:00
#SBATCH -o /tscc/nfs/home/rlancione/jeo/bmat_create_heartfnih.sh.o
#SBATCH -e /tscc/nfs/home/rlancione/jeo/bmat_create_heartfnih.sh.e
#SBATCH -p gold
#SBATCH -q hcg-csd772
#SBATCH -A csd772
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=rlancione@ucsd.edu

# Navigate to oasis scratch
cd /tscc/lustre/scratch/rlancione/
# activate cond env
source ~/miniconda3/etc/profile.d/conda.sh
conda activate renv43


# Input merged object path 
IDIR="/tscc/nfs/home/rlancione/myprojects/Heart_FNIH/Objects/rv_sobj.rds"
# Assay to bulk with
ASSAY="ATAC"
# Dataset column name in object metadata. E.G. "orig.ident"
DSET_COL="donor_demux"
COND_COL="disease_status_sub"
CT_COL="celltypes" 
#Output path
ODIR="/tscc/nfs/home/rlancione/myprojects/Heart_FNIH/DASeq/RMatrices/"
# Project: To name files
PROJ="HeartFNIH"
# Extra covariate columns for meta data. 
CV1="diabetes_status"
CV2="sex"
CV3="age"

# Path to script
echo "run Rscript"
Script='/tscc/nfs/home/rlancione/scToolbox/DESeq2-GO/bmat_create.R'
# args 1. input object path 2. output path 3. project name  
Rscript $Script -sobj $IDIR -a $ASSAY -dcol $DSET_COL -ccol $COND_COL -ctcol $CT_COL --project_name $PROJ -o $ODIR \
	-cov1 $CV1 -cov2 $CV2 -cov3 $CV3


