suppressPackageStartupMessages(library(argparse))
library(Seurat)
library(Signac)
library(stringr)
library(GenomicRanges)
library(parallel)


#### Functions 
### Functions
process_args <- function(){
    # create parser object
    parser <- ArgumentParser(description = "Create seurat object from 10x h5 and collect QC metrics")

    required_arg_group = parser$add_argument_group('flagged required arguments',
        'the script will fail if these args are not included')

    # These args belong to my group of required flagged arguments
    required_arg_group$add_argument("-sobj", "--seurat_object", required = TRUE,
        help="path to dataset h5 file for object creation")

    required_arg_group$add_argument("-pp","--peak_path", required = TRUE,
        help="path to called peaks")

   required_arg_group$add_argument("-dcol","--donor_column", required = TRUE,
        help="column containing id of the patient")

    required_arg_group$add_argument("-ccol","--condition_column", required = TRUE,
        help="column containing condition group")

    required_arg_group$add_argument("-ctcol","--celltype_column", required = TRUE,
        help="column containing celltype label")

    required_arg_group$add_argument("-proj","--project_name", required = TRUE,
        help="name of the DE project. e.g. Lung_BPD_subD352")

    required_arg_group$add_argument("-o","--output_directory", required = TRUE,
        help="path to output directory")

    # optional args
    parser$add_argument("-nc","--ncores", default = 1, 
        help="column containing covariate group label 1")

	parser$add_argument("-cov1","--covariate_column_1",
        help="column containing covariate group label 1")

    parser$add_argument("-cov2","--covariate_column_2",
        help="column containing covariate group label 2")

    parser$add_argument("-cov3","--covariate_column_3",
        help="column containing covariate group label 3")

    args <- parser$parse_args()
    return(args)
}


## Creates bulk matrix from celltype specific peaks
create_ctpeak_bmat <- function(ct, sobj, peak.path, fobjs, ct.col, dset.col, project, outdir){

    peak.path <- peak.path
    peak.files <- list.files(peak.path)
    ct.file <- grep(pattern = paste0("^",ct, ".filterNfixed.peakset"), x = peak.files, value = TRUE)
    ct.path <- paste0(peak.path, ct.file)

    ### Read peaks and create count matrix 
    peaks_df <- read.table(ct.path, header = TRUE)
    gr <- makeGRangesFromDataFrame(peaks_df)

#	print(sobj)
#	print(ct.col)

    Idents(sobj) <- ct.col
    sobj.c <- subset(sobj, idents = ct)

    counts <- FeatureMatrix(
        fragments = fobjs,
        features = gr,
        cells = WhichCells(sobj.c)       
    )

    
    sobj.c[[paste0(ct,"_peaks")]] <-  CreateChromatinAssay(counts = counts, genome = Seqinfo(genome = "hg38"))

    ### Build and write bulk matrix 
    dsets <- unique(sobj.c[[dset.col]])
    dset.frame <- data.frame(row.names = row.names(sobj.c[[paste0(ct,"_peaks")]]$counts))
    Idents(sobj.c) <- dset.col
    for (d in dsets[,1]){
        sobj.t <- subset(sobj.c, idents = d)
        rsums <- rowSums(sobj.t[[paste0(ct,"_peaks")]]$counts)
        dset.frame[d] <- rsums
    }

    ct <- str_replace_all(ct, "/", ".")
    print(paste0("Writting: ", outdir, project, "_", ct, "_bulk_matrix.csv"))
    write.csv(x = dset.frame, file = paste0(outdir, project, "_", ct, "_bulk_matrix.csv"))
    
}


create_meta <- function(sobj = sobj, bulk.by = "orig.ident", col.by = "condition", outdir = NULL, project = "",
                           add_cov1 = NULL, add_cov2 = NULL , add_cov3 = NULL, add_cov4 = NULL, add_cov5 = NULL){
    rows <- unique(sobj[[bulk.by]][,1])
    Idents(sobj) <- bulk.by
    condit <- c()
    ncell_list <- c()
    cov1_labels <- c()
    cov2_labels <- c()
    cov3_labels <- c()
    cov4_labels <- c()
    cov5_labels <- c()
    # Iterate over datasets "rows" to extract conditions ect...
    for (d in rows){
        sobj.t <- subset(sobj, idents = d)
        condit <- append(condit, sobj.t[[col.by]][1,])
        ncells <- length(WhichCells(sobj.t))
        ncell_list <- append(ncell_list, ncells)

        if (!is.null(add_cov1)){
            cov1_labels <- append(cov1_labels, sobj.t[[add_cov1]][1,])
        }
        if (!is.null(add_cov2)){
            cov2_labels <- append(cov2_labels, sobj.t[[add_cov2]][1,])
        }
        if (!is.null(add_cov3)){
            cov3_labels <- append(cov3_labels, sobj.t[[add_cov3]][1,])
        }
        if (!is.null(add_cov4)){
            cov4_labels <- append(cov4_labels, sobj.t[[add_cov4]][1,])
        }
        if (!is.null(add_cov5)){
            cov5_labels <- append(cov5_labels, sobj.t[[add_cov5]][1,])
        }
    }

    meta.frame <- data.frame(row.names = rows, "condition" = condit, "cellquant" = ncell_list)
    if (!is.null(add_cov1)){meta.frame[[add_cov1]] <- cov1_labels}
    if (!is.null(add_cov2)){meta.frame[[add_cov2]] <- cov2_labels}
    if (!is.null(add_cov3)){meta.frame[[add_cov3]] <- cov3_labels}
    if (!is.null(add_cov4)){meta.frame[[add_cov4]] <- cov4_labels}
    if (!is.null(add_cov5)){meta.frame[[add_cov5]] <- cov5_labels}

	meta.frame$type <- "paired-end"	

	if (!is.null(outdir)){
    	write.csv(meta.frame, file = paste0(outdir, project, "_metadata.csv"))
	}

	return(meta.frame)

}




### If name == main... Run...
if (!interactive()) {

	# Call argparse
    args <- process_args()

    # Read Object
    sobj <- readRDS(args$seurat_object)

	# Record input 
	print("****** bulk matrix creation inputs ******")
    print(paste0("sobj.path: ", args$seurat_object))
	print(paste0("peak.path: ", args$peak_path))
	print(paste0("ncores: ", args$ncores))
    print(paste0("dataset column: ", args$donor_column))
    print(unique(sobj[[args$donor_column]][,1]))
    print(paste0("condition column: ", args$condition_column))
    print(unique(sobj[[args$condition_column]][,1]))
    print(paste0("celltype column: ", args$celltype_column))
    print(unique(sobj[[args$celltype_column]][,1]))
    print(paste0(" ouput dir: ", args$output_directory))
    print(paste0("project: ", args$project_name))


    # Create Matrices
	celltypes <- unique(sobj[[args$celltype_column]][,1])
	fobjs <- Fragments(sobj)

#	create_ctpeak_bmat(ct = "Endothelial", sobj = sobj, fobjs = fobjs, project = "LFNIH", outdir = args$output_directory, peak.path = args$peak_path, 
#		ct.col = args$celltype_column, dset.col = args$donor_column)


	mclapply(X = celltypes, FUN = create_ctpeak_bmat, mc.cores = args$ncores, sobj = sobj, fobjs = fobjs , peak.path = args$peak_path,
		 project = args$project_name, outdir = args$output_directory, ct.col = args$celltype_column, dset.col = args$donor_column)

    create_meta(sobj = sobj, bulk.by = args$donor_column, col.by = args$condition_column, outdir = args$output_directory, project = args$project_name,
        add_cov1 = args$covariate_column_1, add_cov2 = args$covariate_column_2, add_cov3 = args$covariate_column_3)

}




