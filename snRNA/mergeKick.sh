#!/bin/bash

#Create merged seurat object and create QC report
#optional arguments:
#  -h, --help            show this help message and exit
#  -ip2 INPUT_PATH2, --input_path2 INPUT_PATH2
#                        optional second path to dir of seurat objects for
#                        merging

#flagged required arguments:
#  the script will fail if these args are not included
#
#  -ip INPUT_PATH, --input_path INPUT_PATH
#                        path to dir of seurat objects for merging
#  -pr PROJECT_NAME, --project_name PROJECT_NAME
#                        Project data belongs too. e.g. Lung_BPD
#  -s SPECIES, --species SPECIES
#                        species of origin for the dataset (mm10 or hg38)
#  -o OUTPUT_PATH, --output_path OUTPUT_PATH
#                        output path for plots and objects

INPUT="/projects/ps-epigen/users/rlan/Karin_Liver/Preprocess/doubfinder_out/1kc_10m_500f/df_objs/"
OUT="/projects/ps-epigen/users/rlan/sandbox/seurat_pipe/merg_out/"
 
#Rscript scRNA_prep.R H5 $H5 SP "mm10" LIB "JB_922_1_2" #OUT $OUT
#Rscript scRNA_prep.R H5 $H5 SPECIES "mm10" OUT $OUT
Rscript scRNA_mergeV2.R -ip $INPUT -pr "pipe_test" -s "mm10" -o $OUT
