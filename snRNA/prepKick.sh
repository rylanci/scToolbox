#!/bin/bash

#flagged required arguments:
#  the script will fail if these args are not included
#
#  -h5 H5_FILE, --h5_file H5_FILE
#                        path to dataset h5 file for object creation
#  -l LIBRARY_ID, --library_id LIBRARY_ID
#                        library id of the sample (JB_232)
#  -s SPECIES, --species SPECIES
#                        species of origin for the dataset (mm10 or hg38)
#  -o OUTPUT_PATH, --output_path OUTPUT_PATH
#                        output path for plots and objects

H5="/projects/ps-epigen/users/rlan/Karin_Liver/Objects/cellranger_out/JB_922_1_2/outs/filtered_feature_bc_matrix.h5"
OUT="/projects/ps-epigen/users/rlan/sandbox/seurat_pipe/prep_out/"
 
#Rscript scRNA_prep.R H5 $H5 SP "mm10" LIB "JB_922_1_2" #OUT $OUT
#Rscript scRNA_prep.R H5 $H5 SPECIES "mm10" OUT $OUT
Rscript scRNA_prep.R -h5 $H5 -l "JB_922_1_2" -s "mm10" -o $OUT
