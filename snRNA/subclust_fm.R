library(Seurat)
library(dplyr)
library(harmony)
library(patchwork)
library(ggplot2)
library(reshape2)
library(gridExtra)
library(argparse)
library(MAST)
source("/tscc/nfs/home/rlancione/scToolbox/snRNA/fm_report.R")


process_args <- function(){
    # create parser object
    parser <- ArgumentParser(description = "Perform DE and visualize output")

    # Create parser arg group for our required args
    required_arg_group = parser$add_argument_group('flagged required arguments',
        'the script will fail if these args are not included')
    
    ### start args 
    required_arg_group$add_argument("-o","--output_path", required = TRUE,
        help="path to output location")

    required_arg_group$add_argument("-n","--name", required = TRUE,
        help="name of analysis e.g. batch_harmony")


    # optional args
    parser$add_argument("-a","--assay",
        help="assay to run fm on", 
        default = "RNA")

    parser$add_argument("-t","--test",
        help="seurat de testing method to use", 
        default = "wilcox")

    parser$add_argument("-w","--workers",
        help="number of cores to run FindMarkers with", 
        default = 12)

    parser$add_argument("-u","--universe", type = "logical",
        help="whether to use clusterProfiler universe (background)", 
        default = "FALSE")

    parser$add_argument("-rfm","--run_fm", type = "logical",
        help="whether to run find markers", 
        default = "TRUE")

    parser$add_argument("-rpr","--run_presto", type = "logical",
        help="whether to run presto", 
        default = "TRUE")

    parser$add_argument("-rcp","--run_cp", type = "logical",
        help="whether to run clusterProfiler", 
        default = "TRUE")


    args <- parser$parse_args()
    return(args)
}

### Record argparse input 
args <- process_args()
print("Reporting FM inputs") 
print(paste0("Output path: ", args$output_path))
print(paste0("Project Name: ", args$name))
print(paste0("Target Assay: ", args$assay))
print(paste0("DE test: ", args$test))
print(paste0("n workers: ", args$workers))


########## Begin running FM report 
print("Initialize FM Report")
# Identify processd object 
file <- list.files(paste0(args$output_path, "Objects/"))
#file <- grep(pattern = paste0(args$name,".RDS"), x = files, value = TRUE)  
print(paste0("Object to read: ", file))

print(paste0("Processing: ", args$output_path, "Objects/", file))
sobj <- readRDS(paste0(args$output_path, "Objects/", file))

dir.create(paste0(args$output_path, "FM_output/"), showWarnings = FALSE)
dir.create(paste0(args$output_path, "FM_output/", args$name), showWarnings = FALSE)
odir <- paste0(args$output_path, "FM_output/", args$name, "/")



if (args$run_fm){
	print("Run FindMarkers")
	if (args$assay == "RNA"){
		DefaultAssay(sobj) <- "RNA"
			sobj <- NormalizeData(sobj)
			res <- suppressMessages(run_fm(sobj = sobj, test = args$test, 
						n.workers = 10))
	} else if (args$assay == "SCT"){
		DefaultAssay(sobj) <- "SCT"
			res <- suppressMessages(run_fm(sobj = sobj, prep.sct = TRUE,
						test = args$test, n.workers = args$workers))
	} else {print("Not supported assay")}

res.de <- fm_report(sobj = sobj, res = res, outdir = odir)
write.table(res, paste0(odir, "fm_res.txt"))

}

if (args$run_presto){
	print("Run Presto")
	run_presto(sobj = sobj, outdir = odir)
}



if (args$universe & args$run_cp){
	print("Run ClusterProfiler")
	print("Utilizing CP universe")
	res.de <- read.table(paste0(odir, "fm_res.txt")) 
	genes <- extract_detected_features(sobj)
	fm_cProfiler(res.de = res.de, outdir = odir, universe = genes)
} else if (args$run_cp){
	print("Run ClusterProfiler")
	fm_cProfiler(res.de = res.de, outdir = odir, universe = NULL)
}



