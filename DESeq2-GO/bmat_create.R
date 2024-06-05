library(Seurat)
library(Signac)
library(dplyr)
library(stringr)
library(argparse)


### Functions 
process_args <- function(){
    # create parser object
    parser <- ArgumentParser(description = "Create seurat object from 10x h5 and collect QC metrics")

    required_arg_group = parser$add_argument_group('flagged required arguments', 
        'the script will fail if these args are not included')

    # These args belong to my group of required flagged arguments
    required_arg_group$add_argument("-sobj", "--seurat_object", required = TRUE,
        help="path to dataset h5 file for object creation")
    
    required_arg_group$add_argument("-dcol","--donor_column", required = TRUE,
        help="column containing id of the patient") 

    required_arg_group$add_argument("-ccol","--condition_column", required = TRUE,
        help="column containing condition group")

    required_arg_group$add_argument("-ctcol","--celltype_column", required = TRUE,
        help="column containing celltype label")

    required_arg_group$add_argument("-p","--project_name", required = TRUE,
        help="name of the DE project. e.g. Lung_BPD_subD352")

    required_arg_group$add_argument("-o","--output_directory", required = TRUE,
        help="path to output directory")

    # optional args
    parser$add_argument("-a", "--assay", default = "RNA", 
        help="the assay containing the counts to bulk")

    parser$add_argument("-cov1","--covariate_column_1",
        help="column containing covariate group label 1")

    parser$add_argument("-cov2","--covariate_column_2",
        help="column containing covariate group label 2")

    parser$add_argument("-cov3","--covariate_column_3",
        help="column containing covariate group label 3")

    args <- parser$parse_args()
    return(args)
}


bulk_matrix <- function(sobj = NULL, assay = "RNA", dset.col = "experiment", condi.col = "conditions", ct.col = "seurat_clusters", outdir = NULL, project = NULL){
    dir.create(outdir, showWarnings = FALSE)

    # Begin iterating over target clusters for matrix creation
    celltypes <- unique(sobj[[ct.col]])[,1]

    Idents(sobj) <- ct.col
    for (ct in celltypes){
        print(paste0("Create matrix for: ", ct))
        # Subset the object for target cluster
        sobj.c <- subset(sobj, idents = c(ct))

        # Extract datasets
        dsets <- unique(sobj.c[[dset.col]])

        # Initialize bulk dataframe & set ident to dataset column
        dset.frame <- data.frame(row.names = row.names(sobj.c[[assay]]$counts))
        Idents(sobj.c) <- dset.col

        # Iterate over datasets in object. Sum rows of raw un-normalize matrix
        # Add result as column to bulk dataframe
        for (d in dsets[,1]){
            sobj.t <- subset(sobj.c, idents = d)
            
            if (!is.null(ncol(sobj.t[[assay]]$counts))){
                rsums <- rowSums(sobj.t[[assay]]$counts)
                dset.frame[d] <- rsums
            }
        }

        # Save bulk dataframe
    	ct <- str_replace_all(ct, "_", ".")
    	ct <- str_replace_all(ct, "/", ".")
    	ct <- str_replace_all(ct, " ", "")

        print(paste0("Dimensions of bulk matrix: ", dim(dset.frame)))
        print(paste0("writing cluster ", ct,".csv"))
        write.csv(x = dset.frame, file = paste0(outdir, project, "_", ct, "_bulk_matrix.csv"))

    }

}


create_DE_meta <- function(sobj = sobj, bulk.by = "orig.ident", col.by = "condition", 
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
    
    return(meta.frame)   
}


### If name == main... Run...
if (!interactive()) {

    # Call argparse
    args <- process_args()
    # Read Object
    sobj <- readRDS(args$seurat_object)

    print("****** bulk matrix creation inputs ******")
    print(paste0("sobj.path: ", args$seurat_object))
    print(paste0("assay: ", args$assay))
    print(paste0("dataset column: ", args$donor_column))
    print(unique(sobj[[args$donor_column]][,1]))
    print(paste0("condition column: ", args$condition_column))
    print(unique(sobj[[args$condition_column]][,1]))
    print(paste0("celltype column: ", args$celltype_column))
    print(unique(sobj[[args$celltype_column]][,1]))
    print(paste0(" ouput dir: ", args$output_directory))
    print(paste0("project: ", args$project_name))

    # Create Matrices
    bulk_matrix(sobj = sobj, assay = args$assay, dset.col = args$donor_column, ct.col = args$celltype_column, condi.col = args$condition_column, 
        outdir = args$output_directory, project = args$project_name)

    DO_meta <- create_DE_meta(sobj = sobj, bulk.by = args$donor_column, col.by = args$condition_column, 
        add_cov1 = args$covariate_column_1, add_cov2 = args$covariate_column_2, add_cov3 = args$covariate_column_3)
    DO_meta$type <-  "paired-end" 

    write.csv(DO_meta, file = paste0(args$output_directory, args$project_name, "_metadata.csv"))
}


