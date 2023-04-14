suppressMessages(library(Seurat))
suppressMessages(library(dplyr))
suppressMessages(library(ggplot2))
suppressMessages(library(DESeq2))
suppressMessages(library(enrichR))
suppressMessages(library(stringr))
suppressMessages(library(tableHTML))
suppressMessages(library(argparse))
suppressMessages(library(org.Hs.eg.db))
suppressMessages(library(clusterProfiler))
suppressMessages(library(enrichplot))

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

# optional args
    parser$add_argument("-comp","--comparisons",
        help="path to tsv containing comparisons. Use header: Target  Reference")

    parser$add_argument("-form","--formula",
        help="a string representing the DESeq formula. Default is: ~ condition")

    
    args <- parser$parse_args()
    return(args)
}


filt.features <- function(blk.matrix, meta, n, comp_pair){
    # will store our list of features that passed theshold
    pf <- c()

    g1.donors <- rownames(meta[meta$condition == comp_pair[1],])
    g2.donors <- rownames(meta[meta$condition == comp_pair[2],])
    
    # Iterate over rows (features) of bulk matrix
    for (r in seq(1, nrow(blk.matrix))){
      
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
    
    return(pf)
}


clusterProfiler_GO <- function(DE.table, cell_label, universe = NULL, log2FC_threshold = .58,
                            padj_threshold = 0.1, outdir = "", organism = "human", comp){
    print("Intiating cluster profiler")
    # Create out dir 
    outdir <- paste0(outdir, "CP_OUT/")
    dir.create(outdir, showWarnings = FALSE)
    dir.create(paste0(outdir,"csv_tables/"), showWarnings = FALSE)
    dir.create(paste0(outdir,"html_tables/"), showWarnings = FALSE)
    dir.create(paste0(outdir,"plots/"), showWarnings = FALSE)

    # Extract significant DE's from DE table
    final <- DE.table[DE.table$padj < padj_threshold,]

    up.df <- final[final$log2FoldChange > log2FC_threshold,]
    dwn.df <- final[final$log2FoldChange < -log2FC_threshold,]
    comb.df <- final[final$log2FoldChange < -log2FC_threshold | 
                            final$log2FoldChange > log2FC_threshold,]

    results_list <- list(up.df, dwn.df, comb.df)
    results_names <- c("up", "dwn", "cmb")
    
    # handle organism here
    if(organism == "human"){
        org.db <- org.Hs.eg.db
    } else if(organism == "mouse"){
        org.db <- org.Mm.eg.db
    }
    
    # Iterate through results list
    # Run each go list for each results list 
    for (i in seq(1,3)){
        print(paste0("Starting ", cell_label,"_", results_names[i], " GO Enrichment"))

        # results_list[i]$X extracts gene names from each DE table
        ego_cc <- enrichGO(gene = rownames(results_list[[i]]),
                    OrgDb = org.db,
                    keyType = 'SYMBOL',
                    universe = universe, 
                    ont = "CC",
                    pAdjustMethod = "BH",
                    pvalueCutoff = 1,
                    qvalueCutoff = 1)

        ego_bp <- enrichGO(gene = rownames(results_list[[i]]),
                    OrgDb = org.db,
                    keyType = 'SYMBOL',
                    universe = universe, 
                    ont = "BP",
                    pAdjustMethod = "BH",
                    pvalueCutoff = 1,
                    qvalueCutoff = 1)

        ego_mf <- enrichGO(gene = rownames(results_list[[i]]),
                    OrgDb = org.db,
                    keyType = 'SYMBOL',
                    universe = universe, 
                    ont = "MF",
                    pAdjustMethod = "BH",
                    pvalueCutoff = 1,
                    qvalueCutoff = 1)
	
	# null is returned by enrichGO if no ontologies are detected. These conditionals catch that so we can write stuff 
	if (!is.null(ego_cc)){
		write.table(ego_cc@result, file = paste0(outdir,"csv_tables/", cell_label, "_", comp, "_", results_names[i], "_CC_go.csv"))
		write_tableHTML(tableHTML(ego_cc@result), file = 
			paste0(outdir,"html_tables/", cell_label, "_", comp, "_", results_names[i], "_CC_go.html"))

		p1 <- dotplot(ego_cc, showCategory=20) + ggtitle("cellular component")
		fn1 <- paste0(cell_label, "_", comp, "_", results_names[i], "_CC_DTPLT.pdf")
                ggsave(filename = fn1, plot = p1, device = "pdf", path = paste0(outdir, "plots/"), height = 10)

        }
	if (!is.null(ego_bp)){
		 write.table(ego_bp@result, file = paste0(outdir,"csv_tables/", cell_label, "_", comp, "_", results_names[i], "_BP_go.csv"))
		 write_tableHTML(tableHTML(ego_bp@result), file = 
			paste0(outdir,"html_tables/", cell_label, "_", comp, "_", results_names[i], "_BP_go.html"))

		 p2 <- dotplot(ego_bp, showCategory=20) + ggtitle("biological process")
		 fn2<- paste0(cell_label, "_", comp, "_", results_names[i], "_BP_DTPLT.pdf")
		 ggsave(filename = fn2, plot = p2, device = "pdf", path = paste0(outdir, "plots/"),  height = 10)

        }
	if (!is.null(ego_mf)) {
		write.table(ego_mf@result, file = paste0(outdir,"csv_tables/", cell_label, "_", comp, "_", results_names[i], "_MF_go.csv"))
		write_tableHTML(tableHTML(ego_mf@result), file = 
			paste0(outdir,"html_tables/", cell_label, "_", comp, "_", results_names[i], "_MF_go.html"))

		p3 <- dotplot(ego_mf, showCategory=20) + ggtitle("molecular function")
                fn3 <- paste0(cell_label, "_", comp, "_", results_names[i], "_MF_DTPLT.pdf")
		ggsave(filename = fn3, plot = p3, device = "pdf", path = paste0(outdir, "plots/"), height = 10)
        }
        
                        
    }
}


runEnrichR <- function(deseqCSV_out, clusterNum, condPair, adj_pval, log2FC_threshold, output_dir){
  # Create enichR output dir
  dir.create(paste0(output_dir,"Enrichr_out"), showWarnings = FALSE)

  ### create UP- and DOWN-regulated csv files for EnrichR
  deseq_results <- read.csv(deseqCSV_out)
  # set adjusted pval cutoff
  deseq_results <- deseq_results[deseq_results$padj < adj_pval, ]

  # use deseq results and separate by UP and DOWN-regulated for enrichr
  dwn <- deseq_results[deseq_results$log2FoldChange < -log2FC_threshold, ]
  up <- deseq_results[deseq_results$log2FoldChange > log2FC_threshold, ]
  combined <- deseq_results[ which( deseq_results$log2FoldChange < -log2FC_threshold | deseq_results$log2FoldChange > log2FC_threshold ), ]

  # Save enrichr input in separate dir
  dir.create(paste0(output_dir,"DESeq_SIG/"), showWarnings = FALSE)
  up_fp <- paste(output_dir, "DESeq_SIG/", "DESeq_sigUP_c", clusterNum, "_",condPair,".csv", sep = "")
  dwn_fp <- paste(output_dir, "DESeq_SIG/", "DESeq_sigDWN_c", clusterNum, "_",condPair, ".csv", sep = "")
  write.csv(up, up_fp, row.names=FALSE, quote=FALSE)
  write.csv(dwn, dwn_fp, row.names=FALSE, quote=FALSE)

  ############ RUNNING ENRICHR ############
  E.libs <- c("MSigDB_Hallmark_2020", "GO_Biological_Process_2018", "GO_Biological_Process_2021", "KEGG_2021_Human", "GO_Molecular_Function_2021",
		"Human_Gene_Atlas", "Genes_Associated_with_NIH_Grants")

  # iterate through enrichr gene lists. Run Enrichr on DESeq out and write files
  for (i in seq(1,6)){
	    dir.create(paste0(output_dir, "Enrichr_out/", E.libs[i]))

      # DOWNREG
	    dwn_df <- suppressMessages(data.frame(enrichr(dwn[,1], E.libs[i])))
    	if (nrow(dwn_df) > 0 ){
         # Correct for commas in GO terms
    		 dc1.R <- str_replace_all(string = dwn_df[,1], pattern = ",", replacement = "/")
         dwn_df[,1] <- dc1.R
    	}
     	dwnfile <- paste(output_dir, "Enrichr_out/", E.libs[i], "/DWNreg_", clusterNum, "_", condPair, ".csv", sep = "")
    	write.csv(dwn_df, dwnfile, row.names=FALSE, quote=FALSE)
    	if (nrow(dwn_df > 1)){
    	   dwnfile_html <- paste(output_dir, "Enrichr_out/", E.libs[i], "/DWNreg_", clusterNum, "_", condPair, ".html", sep = "")
    	   write_tableHTML(tableHTML(dwn_df), file = dwnfile_html)
    	}

    	# UPREG
    	up_df <- suppressMessages(data.frame(enrichr(up[,1], E.libs[i])))
    	if (nrow(up_df) > 0 ){
    	   uc1.R <- str_replace_all(string = up_df[,1], pattern = ",", replacement = "/")
    	   up_df[,1] <- uc1.R
    	}
    	upfile <- paste(output_dir, "Enrichr_out/", E.libs[i], "/UPreg_", clusterNum, "_", condPair, ".csv", sep = "")
    	write.csv(up_df, upfile, row.names=FALSE, quote=FALSE)
        if (nrow(up_df > 1)){
    	   upfile_html <- paste(output_dir, "Enrichr_out/", E.libs[i], "/UPreg_", clusterNum, "_", condPair, ".html", sep = "")
    	   write_tableHTML(tableHTML(up_df), file = upfile_html)
    	}

    	# COMBINED
     	comb_df <- suppressMessages(data.frame(enrichr(combined[,1], E.libs[i])))
    	if (nrow(comb_df) > 0 ){
    	   cc1.R <- str_replace_all(string = comb_df[,1], pattern = ",", replacement = "/")
    	   comb_df[,1] <- cc1.R
    	}
    	cfile <- paste(output_dir, "Enrichr_out/", E.libs[i], "/COMB_", clusterNum, "_", condPair, ".csv", sep = "")
    	write.csv(comb_df, cfile, row.names=FALSE, quote=FALSE)
     	if (nrow(comb_df > 1)){
        cfile_html <- paste(output_dir, "Enrichr_out/", E.libs[i], "/COMB_", clusterNum, "_", condPair, ".html", sep = "")
        write_tableHTML(tableHTML(comb_df), file = cfile_html)
      }
   }
}


pseudoBulk_DE <- function(sobj.path, matrices.path, ct.col = "seurat_clusters", padj_threshold, log2FC_threshold,
				output_dir, project_name, comps = "pairwise", formula = "~ condition", nfeatures = 0){
  print("Initiating DESeq Pipeline")
  seuObj <- readRDS(sobj.path)

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

      ### Extract our comparisons it iterate through. Either pairwise combo of all or user defined.
      # if comps is null, run pairwise comparisons. null argparse as input to pseudobulk_DE overwrites default "pairwise"
      if (is.null(comps)){
      	        # DEFAULT. Run pairwise
        	pair_list <- combn(unique(metadata$condition), 2)
         	num_pairs <- ncol(pair_list)
      } else {
        	comp.table <- read.table(comps, header = T, sep = "\t")
        	num_pairs <- nrow(comp.table)
        	pair_list <- t(comp.table)
      }
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
	  pass.features <- filt.features(blk.matrix = blk.mat.comp, meta = sub.meta.comp, n = nfeatures, comp_pair = cpair)
	  
	  blk.mat.comp.filt <- blk.mat.comp[pass.features,]
	  print(paste0("nrow filt matrix:", nrow(blk.mat.comp.filt))) 

      	  # save filtered matrix 
	  condPair <- paste0(cond1, "v", cond2)
	  dir.create(paste0(output_dir, "Filtered_Matrices/"), showWarnings = FALSE) 	  
          filt_mat_file <- paste0(output_dir, "Filtered_Matrices/", c, "_" , condPair, "_filt_bulk_mat.csv")
          write.csv(blk.mat.comp.filt, file = filt_mat_file, quote=FALSE)


      	  if (length(unique(sub.meta.comp$condition)) > 1 ){
 	           ### If no formula is input. Only test condition. If formula is present, try to use it. If it doesn't work, only test condition
      		   form_used <- "~ condition"
		   if(is.null(formula)){
          	        # DEFAULT. Simply test each condition
          		deseq_results <- DESeqDataSetFromMatrix(countData = blk.mat.comp.filt, colData = sub.meta.comp, design = ~ condition)
      	 	   } else {
          		# try to use formula  
			print("try to use form")
			form_used <- formula
			an.error.occured <- FALSE
			tryCatch({deseq_results <- DESeqDataSetFromMatrix(countData = blk.mat.comp.filt, 
					colData = sub.meta.comp, design = as.formula(formula)) }
			, error = function(e) {an.error.occured <<- TRUE})
		
			# if error run default
			if (an.error.occured == TRUE){
				# Change form used back to default since input didn't work 
				print("error occured, use defualt form")	
				form_used <- "~ condition"
				deseq_results <- DESeqDataSetFromMatrix(countData = blk.mat.comp.filt, colData = sub.meta.comp, design = ~ condition)
			}
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

  	           ### Running EnrichR
  	           # saving all DEG counts and adding to list - text file for each cluster
#  	           runEnrichR(deseq_matrix, c, condPair, padj_threshold, log2FC_threshold, output_dir)
		   ### Running Cluster profiler
	   	   
		   clusterProfiler_GO(DE.table = final, cell_label = c, universe = pass.features, log2FC_threshold = log2FC_threshold,
                            padj_threshold = padj_threshold, outdir = output_dir, organism = "human", comp = condPair)

		   ### Create DE Summary
               	   summarylist <- c()
                   # set adjusted pval cutoff.
               	   final <- final[final$padj < padj_threshold, ]
               	   # use deseq results and separate by UP and DOWN-regulated for enrichr
               	   dwn <- final[final$log2FoldChange < log2FC_threshold, ]
               	   summarylist <- c(summarylist, paste(c, "DOWN", condPair, nrow(dwn), form_used, sep = "\t"))

               	   up <- final[final$log2FoldChange > log2FC_threshold, ]
               	   summarylist <- c(summarylist, paste(c, "UP", condPair, nrow(up), form_used, sep = "\t"))

               	   combined <- final[which(final$log2FoldChange < -log2FC_threshold | final$log2FoldChange > log2FC_threshold ), ]
               	   summarylist <- c(summarylist, paste(c, "COMBINED", condPair, nrow(combined), form_used, sep = "\t"))

  		   write.table(summarylist, paste(output_dir, "DEG_countSummary.txt", sep=""), quote=FALSE, row.names=FALSE, col.names=FALSE, append=TRUE)

	       } else{
		   print("not enough replicates to run DESeq for comparison")
              	   # Write not enough donors to make comparison
               	   summarylist <- c(summarylist, paste(c, "COMBINED", condPair, "Not enough donors per comparison", sep = "\t"))
  		   write.table(summarylist, paste(output_dir, "DEG_countSummary.txt", sep=""), quote=FALSE, row.names=FALSE, col.names=FALSE, append=TRUE)

	       }

          print(paste("DONE: Celltype ", c, " | Conditions ", cond1, " vs. ", cond2, sep=""))
        }
    }
}



########### MAIN ###########
## Aquire inputs
args <- process_args()

### Record input
print("****** bulk matrix creation inputs ******")
print(paste0("sobj.path: ", args$seurat_object))
print(paste0("matrices.path: ", args$matrix_path))
print(paste0("adjusted.pval: ", args$padj_threshold))
print(paste0("log2FC.thres: ", args$log2FC_threshold))
print(paste0("output.dir: ", args$output_path))
print(paste0("project.name: ", args$project_name))
print(paste0("celltype column: ", args$celltype_column))
print(paste0("comparison file: ", args$comparisons))
print(paste0("n features: ", args$n_features)) 
if(is.null(args$formula) == TRUE){
	print("formula: ~ condition")
} else{
	print(paste0("formula: ", args$formula))
}
if(is.null(args$comparisons) == TRUE){
	print("comparisons: pairwise")
} else{
	print("comparisons: ")
  print(read.table(args$comparisons, header = T, sep = "\t"))
}

pseudoBulk_DE(sobj.path = args$seurat_object, ct.col = args$celltype_column, matrices.path = args$matrix_path,
                padj_threshold = args$padj_threshold, log2FC_threshold = args$log2FC_threshold, output_dir = args$output_path,
                project_name = args$project_name, comps = args$comparisons, formula = args$formula, nfeatures = args$n_features)


# Features to add:
# Enrichr summary
# Fix DESeq sig/no logfc count summary
