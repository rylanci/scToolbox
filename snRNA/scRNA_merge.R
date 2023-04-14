#library(Seurat)
#library(dplyr)
#library(Nebulosa)
#library(patchwork)
#library(viridis)
#library(future)
#library(harmony)
#library(preprocessCore)
#library(gridExtra)
library(argparse)
source("~/scripts/utils.R")

process_args <- function(){
# create parser object
parser <- ArgumentParser()

parser$add_argument("-ip", "--input_path",
    help="path to dir of seurat objects for merging")

parser$add_argument("-op", "--output_path",
    help="path to desired output location")


#args <- parser$parse_args(c('--input_path','--output_path'))
args <- parser$parse_args()


return(args)
}

args <- process_args()

print(args$input_path)
print(args$output_path)
