suppressMessages(library(Seurat))
suppressMessages(library(dplyr))
suppressMessages(library(ggplot2))
suppressMessages(library(DESeq2))
suppressMessages(library(stringr))
suppressMessages(library(tableHTML))
suppressMessages(library(argparse))
suppressMessages(library(gridExtra))

### Script purpose
# Taking bulk matrices and metadata created from RNA_MatrixCreation.R script
# Running DESeq2 and EnrichR across all clusters and every pairwise condition combination


### Functions
process_args <- function(){
    # create parser object
    parser <- ArgumentParser(description = "Create seurat object from 10x h5 and collect QC metrics")

    # Create parser arg group for our required args
    required_arg_group = parser$add_argument_group('flagged required arguments',
        'the script will fail if these args are not included')

    # These args belong to my group of required flagged arguments
    required_arg_group$add_argument("-sobj", "--seurat_object", required = TRUE,
        help="path to dataset h5 file for object creation")

    required_arg_group$add_argument("-ctcol","--celltype_column", required = TRUE,
        help="column containing id of the patient")

    required_arg_group$add_argument("-mpath","--matrix_path", required = TRUE,
        help="column containing condition group")

    required_arg_group$add_argument("-ap","--padj_threshold", required = TRUE,
        help="column containing celltype label")

    required_arg_group$add_argument("-lfc","--log2FC_threshold", required = TRUE,
        help="column containing celltype label")

    required_arg_group$add_argument("-proj","--project_name", required = TRUE,
        help="name of the DE project. e.g. Lung_BPD_subD352")

    required_arg_group$add_argument("-o","--output_path", required = TRUE,
        help="path to output directory")

    required_arg_group$add_argument("-nf","--n_features", required = TRUE,
        help="n features must be observed in every donor for at least one condition group")

    required_arg_group$add_argument("-comp","--comparisons",
        help="path to tsv containing comparisons and formulas. Use header: Target  Reference  Formula")

    
    args <- parser$parse_args()
    return(args)
}


filt.features.avg <- function(blk.matrix, meta, avg.thresh, comp_pair){
    # will store our list of features that passed theshold
    pf <- c()
    g1.donors <- rownames(meta[meta$condition == comp_pair[1],])
    g2.donors <- rownames(meta[meta$condition == comp_pair[2],])
    
    # Iterate over rows (features) of bulk matrix
    for (r in seq(1, nrow(blk.matrix))){
      
	g1.values <- as.integer(blk.matrix[r, g1.donors])
	g2.values <- as.integer(blk.matrix[r, g2.donors])
	g1.sum <- sum(g1.values)
        g2.sum <- sum(g2.values)

	# Average of 10 required in either group, and at least 1 count detected in both groups
	# Add that feature to the list
    # Needs the extra parentheses to work correctly 
	if ((mean(g1.values) >= avg.thresh | mean(g2.values) >= avg.thresh) & g1.sum >= 1 & g2.sum >= 1){
	    pf <- append(pf, rownames(blk.matrix[r,]))
	}
    }
    return(pf)
}


# Avg feature filter using avg(all donors)
filt.features.avg_v2 <- function(blk.matrix, meta, avg.thresh, comp_pair){
    # will store our list of features that passed theshold
    pf <- c()
    g1.donors <- rownames(meta[meta$condition == comp_pair[1],])
    g2.donors <- rownames(meta[meta$condition == comp_pair[2],])
    
    # Iterate over rows (features) of bulk matrix
    for (r in seq(1, nrow(blk.matrix))){
                                
        values <- as.integer(blk.matrix[r, c(g1.donors, g2.donors)])
        
        if (mean(values) >= avg.thresh){
            pf <- append(pf, rownames(blk.matrix[r,]))
        }
    }

    return(pf)
}


# Avg feature filter using avg(all donors) + require at least one count from each donor group
filt.features.avg_v3 <- function(blk.matrix, meta, avg.thresh, comp_pair){
    # will store our list of features that passed theshold
    passed.features <- c()
    # vectors of donor ids pertaining to the comparison
    g1.donors <- rownames(meta[meta$condition == comp_pair[1],])
    g2.donors <- rownames(meta[meta$condition == comp_pair[2],])
    
    # Iterate over rows (features) of bulk matrix
    for (r in seq(1, nrow(blk.matrix))){
                                
        g1.values <- as.integer(blk.matrix[r, g1.donors])
        g2.values <- as.integer(blk.matrix[r, g2.donors])
        values <- c(g1.values, g2.values)
        
        g1.sum <- sum(g1.values)
        g2.sum <- sum(g2.values)
        
        if (mean(values) >= avg.thresh & g1.sum >= 1 & g2.sum >= 1){
            passed.features <- append(passed.features, rownames(blk.matrix[r,]))
        }
    }

    return(passed.features)
}



filt.features <- function(blk.matrix, meta, n, comp_pair){
    # will store our list of features that passed theshold
    pf <- c()

    ### unique condition for healed vs healed_control in Lung_BPD
    if ("BPD_healed" %in% comp_pair){
	print("using BPD_healed n1 filter")
	healed.donors <- rownames(meta[meta$condition == "BPD_healed",])
	control.donors <- rownames(meta[meta$condition == "BPD_healed_control",])
    }

    g1.donors <- rownames(meta[meta$condition == comp_pair[1],])
    g2.donors <- rownames(meta[meta$condition == comp_pair[2],])
    
    # Iterate over rows (features) of bulk matrix
    for (r in seq(1, nrow(blk.matrix))){
      
	### Special condition says that at least 1 count must be detected in all donors of healed group
	### hopefully prefents DE of feautres with no healed counts
        if ("BPD_healed" %in% comp_pair){
	    healed.values <- blk.matrix[r, healed.donors]
	    control.values <- blk.matrix[r, control.donors]

	    n2 <- 1

	    healed.bool.list <- healed.values >= n
	    control.bool.list <- control.values >= n2
	    healed.bool.n2 <- healed.values >= n2
	    control.bool.n2 <- control.values >= n
	
	    healed.ntrue <- sum(healed.bool.list)
	    control.ntrue <- sum(control.bool.list)
	    healed.ntrue.n2 <- sum(healed.bool.n2)
	    control.ntrue.n2 <- sum(control.bool.n2)


	    if (healed.ntrue >= length(healed.donors) & control.ntrue >= length(control.donors)){
		    pf <- append(pf, rownames(blk.matrix[r,]))
	    } else if (healed.ntrue.n2 >= length(healed.donors) & control.ntrue.n2 >= length(control.donors)){
		    pf <- append(pf, rownames(blk.matrix[r,]))
	    }
	} else {
	    g1.values <- blk.matrix[r, g1.donors]
	    g2.values <- blk.matrix[r, g2.donors]

	    # Vectorize the query across the columns 
	    g1.bool.list <- g1.values >= n
	    g2.bool.list <- g2.values >= n

	    # Sum the trues from our query
	    g1.ntrue <- sum(g1.bool.list)
	    g2.ntrue <- sum(g2.bool.list)

	    # If there were at least n counts in all donors in one of the condition groups
	    # Add that feature to the list
	    if (g1.ntrue >= length(g1.donors) | g2.ntrue >= length(g2.donors)){
		    pf <- append(pf, rownames(blk.matrix[r,]))
	    }
	}
    }
    
    return(pf)
}


pseudoBulk_DE <- function(sobj.path, matrices.path, ct.col = "seurat_clusters", padj_threshold, log2FC_threshold,
				output_dir, project_name, comps,  nfeatures = 0){
    print("Initiating DESeq Pipeline")
    seuObj <- readRDS(sobj.path)
    comp.table <- read.table(comps, header = T, sep = "\t")

    # Format numerical values
    padj_threshold <- as.double(padj_threshold)
    log2FC_threshold <- as.double(log2FC_threshold)

    ### Create output outdirs
    dir.create(output_dir, showWarnings = FALSE)
    dir.create(paste0(output_dir,"DESeq_out"), showWarnings = FALSE)

    if (ct.col == "seurat_clusters"){
        cell_type <- sort.int(unique(seuObj[[ct.col]])[,1])
    } else {
        cell_type <- unique(seuObj[[ct.col]])[,1]
        cell_type <- str_replace_all(cell_type, "/", ".")
    }

    ### Iterate through celltypes. Run DESeq. Iterate through comparisons. Extract results. Call Go.
    for (c in cell_type){
        print(paste("STARTING: Celltype", c))

        ### Reading files
      	blk.matrix <- read.csv(paste(matrices.path, project_name, "_", c, "_bulk_matrix.csv", sep=""), sep=',',
            header=TRUE, row.names=1, stringsAsFactors=FALSE)

      	metadata <- read.csv(paste(matrices.path, project_name,"_metadata.csv", sep=""), sep=',', header=TRUE,
            row.names = 1, stringsAsFactors=FALSE)

        ### subset metadata for donors in celltype blk matrix
        sub.meta <- subset(metadata, row.names(metadata) %in% colnames(blk.matrix))

	    ### No longer has pairwise comparison option
	    ### Manual comparisons in our table + formulas now mandatory 
    	comparisons <- comp.table[,c("Target", "Reference")]
      	num_pairs <- nrow(comp.table)
      	pair_list <- t(comp.table)
    	formulas <- comp.table[,"Formula"]

        # initialize results summary with header 
#        header <- paste0("Celltype", "Direction", "Comparison", "nDEG", "Formula", "Comp1 Donors", "Comp2 Donors", sep = "\t")
#        write.table(header, paste(output_dir, "DEG_countSummary.txt", sep=""), quote=FALSE, row.names=FALSE, col.names=FALSE)

	
      	### Begin iterating through comparisons to extract results for each group
        for (i in 1:num_pairs){
       	    	cpair <- pair_list[,i]
	        cond1 <- cpair[1]
	        cond2 <- cpair[2]

            # subset metadata for desired conditions. ### Still need to check whether enough donors for comparison
            sub.meta.comp <- sub.meta[which(sub.meta$condition==cond1 | sub.meta$condition==cond2), ]
	        blk.mat.comp <- subset(blk.matrix, select = row.names(sub.meta.comp))
            # If only one condition is represented by donors existing in a celltype. Skp DESeq

            ###### Here is where features we test will be filtered by a specified N.
            ### These filtered matrices will be saved for convinience
            # if filter.features == TRUE do this. if not. pass blk.mat.comp to blk.mat.comp.filt
            print(paste0("identify valid features: ", cond1, "v", cond2))
            nfeatures <- as.integer(nfeatures)
            # use average feature filter. for normal filter use n. for avg use avg.thresh param
            pass.features <- filt.features.avg(blk.matrix = blk.mat.comp, meta = sub.meta.comp, avg.thresh = nfeatures, comp_pair = cpair)

            blk.mat.comp.filt <- blk.mat.comp[pass.features,]
            print(paste0("nrow filt matrix:", nrow(blk.mat.comp.filt)))

      	    # save filtered matrix 
            condPair <- paste0(cond1, "-", cond2)
            dir.create(paste0(output_dir, "Filtered_Matrices/"), showWarnings = FALSE) 	  
            filt_mat_file <- paste0(output_dir, "Filtered_Matrices/", c, "_" , condPair, "_filt_bulk_mat.csv")
            write.csv(blk.mat.comp.filt, file = filt_mat_file, quote=FALSE)

            ### instead of making sure both conditions are present. count the number of donors present in each condition
            ### Then only run de if greater than 1 donor is present for each group
            ### and append the donors for each group to the summary 
            sub.meta.comp1 <- sub.meta[which(sub.meta$condition==cond1), ]
            sub.meta.comp2 <- sub.meta[which(sub.meta$condition==cond2), ]
            blk.mat.comp1 <- subset(blk.mat.comp, select = rownames(sub.meta.comp1)) 
            blk.mat.comp2 <- subset(blk.mat.comp, select = rownames(sub.meta.comp2)) 
            comp1.donors <- colnames(blk.mat.comp1)
            comp2.donors <- colnames(blk.mat.comp2)
            print(paste0("Comp1 donors: ", comp1.donors))
            print(paste0("Comp2 donors: ", comp2.donors))

            ### ************ else if condi = active | chronic batch add TOD as covariate 
            if (length(comp1.donors) > 1 & length(comp2.donors) > 1 & length(pass.features) >= 100){
                # try to use formula supplied in table 
                print("try to use formula supplied in table")
                form_used <- formulas[i]
                formula <- formulas[i]	
                an.error.occured <- FALSE
		        tryCatch({deseq_results <- DESeqDataSetFromMatrix(countData = blk.mat.comp.filt, 
			            colData = sub.meta.comp, design = as.formula(formula)) }
		                , error = function(e) {an.error.occured <<- TRUE})
		    # if error run default
		    if (an.error.occured == TRUE){
                # Change form used back to default since input didn't work 
                print("error occured, use defualt form")	
                form_used <- "~ condition"
                deseq_results <- DESeqDataSetFromMatrix(countData = blk.mat.comp.filt, 
                        colData = sub.meta.comp, design = ~ condition)
		    }

            deseq_results <- DESeq(deseq_results)
            ############ Declare comparison so direction of results are known ##############
            final <- results(deseq_results, contrast = c("condition", cond1, cond2))
            # keep only if both padj and log2FC is NOT NA
            final <- final[is.na(final$log2FoldChange) == FALSE, ]
            final <- final[is.na(final$padj) == FALSE, ]

            # write deseq matrix to file
            deseq_matrix <- paste(output_dir, "DESeq_out/", "DESeq2_", c, "_" , condPair, ".csv", sep="")
            write.csv(final, file=deseq_matrix, quote=FALSE)

            ### Create DE Summary
            summarylist <- c()
            # set adjusted pval cutoff.
            pfinal <- final[final$padj < padj_threshold, ]
            # use deseq results and separate by UP and DOWN-regulated for enrichr
            dwn <- pfinal[pfinal$log2FoldChange < log2FC_threshold, ]
            summarylist <- c(summarylist, paste(c, "DOWN", condPair, nrow(dwn), form_used, sep = "\t"))

            up <- pfinal[pfinal$log2FoldChange > log2FC_threshold, ]
            summarylist <- c(summarylist, paste(c, "UP", condPair, nrow(up), form_used, sep = "\t"))

            combined <- pfinal[pfinal$log2FoldChange < -log2FC_threshold | pfinal$log2FoldChange > log2FC_threshold, ]
            summarylist <- c(summarylist, paste(c, "COMBINED", condPair, nrow(combined), form_used, sep = "\t"))

            write.table(summarylist, paste(output_dir, "DEG_countSummary.txt", sep=""), quote=FALSE, row.names=FALSE, col.names=FALSE, append=TRUE)

            } else {
                # Deseq won't run if these conditions aren't met
                print("not enough replicates to run DESeq for comparison")
                print("or less than 100 features passed filter")
                # Write not enough donors to make comparison
                summarylist <- c()
                summarylist <- c(summarylist, paste(c, "COMBINED", condPair, "NA", "NA", sep = "\t"))
                write.table(summarylist, paste(output_dir, "DEG_countSummary.txt", sep=""), quote=FALSE, row.names=FALSE, col.names=FALSE, append=TRUE)
            }
	  
            # end condition loop 
            print(paste("DONE: Celltype ", c, " | Conditions ", cond1, " vs. ", cond2, sep=""))
        }    
        # end celltype loop
    }
    # end pseudobulk de func
}



########### MAIN ###########
## Aquire inputs
args <- process_args()

### Record input
print("****** DE Pipeline Inputs ******")
print(paste0("sobj.path: ", args$seurat_object))
print(paste0("matrices.path: ", args$matrix_path))
print(paste0("adjusted.pval: ", args$padj_threshold))
print(paste0("log2FC.thres: ", args$log2FC_threshold))
print(paste0("output.dir: ", args$output_path))
print(paste0("project.name: ", args$project_name))
print(paste0("celltype column: ", args$celltype_column))
print(paste0("comparison file: ", args$comparisons))
print(paste0("n features: ", args$n_features)) 
print("comparisons: ")
print(read.table(args$comparisons, header = T, sep = "\t"))

pseudoBulk_DE(sobj.path = args$seurat_object, ct.col = args$celltype_column, matrices.path = args$matrix_path,
                padj_threshold = args$padj_threshold, log2FC_threshold = args$log2FC_threshold, output_dir = args$output_path,
                project_name = args$project_name, comps = args$comparisons, nfeatures = args$n_features)


# Features to add:
# Enrichr summary
# Fix DESeq sig/no logfc count summary
