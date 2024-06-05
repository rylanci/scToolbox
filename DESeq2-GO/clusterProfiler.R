suppressPackageStartupMessages(library(clusterProfiler))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(org.Hs.eg.db))
suppressPackageStartupMessages(library(org.Mm.eg.db))

clusterProfiler_GO <- function(DE.table, fname, universe = NULL, log2FC_threshold = .58,
                            padj_threshold = 0.1, outdir = "", organism = "human"){
    print("Intiating cluster profiler")
    # Create out dir
    outdir <- paste0(outdir, "CP_OUT/")
    dir.create(outdir, showWarnings = FALSE)
    # create dirs for each catagory of output
    dir.create(paste0(outdir,"csv_tables/"), showWarnings = FALSE)
    dir.create(paste0(outdir,"html_tables/"), showWarnings = FALSE)
    dir.create(paste0(outdir,"plots/"), showWarnings = FALSE)
    dir.create(paste0(outdir,"cp_objects/"), showWarnings = FALSE)
    # create separate directories for the various ontology platforms (GO, KEGG, Pathway)
    dir.create(paste0(outdir, "csv_tables/GO/"), showWarnings = FALSE)
    dir.create(paste0(outdir, "csv_tables/KEGG/"), showWarnings = FALSE)
    dir.create(paste0(outdir, "csv_tables/REACT/"), showWarnings = FALSE)
    dir.create(paste0(outdir, "html_tables/GO/"), showWarnings = FALSE)
    dir.create(paste0(outdir, "html_tables/KEGG/"), showWarnings = FALSE)
    dir.create(paste0(outdir, "html_tables/REACT/"), showWarnings = FALSE)
    dir.create(paste0(outdir, "cp_objects/GO/"), showWarnings = FALSE)
    dir.create(paste0(outdir, "cp_objects/KEGG/"), showWarnings = FALSE)
    dir.create(paste0(outdir, "cp_objects/REACT/"), showWarnings = FALSE)

    ### exclude mito features from DE table
#    mito.features <- c("MT-ND1", "MT-ND2", "MT-CO1", "MT-CO2", "MT-ATP6", "MT-CO3", "MT-ND3", 
#			"MT-ND4L", "MT-ND4", "MT-ND5", "MT-CYB") 
#    DE.table <- DE.table[!(DE.table$gene %in% mito.features),]
     
    # Extract significant DE's from DE table
    final <- DE.table[DE.table$padj < padj_threshold,]

    log2FC_threshold <- as.double(log2FC_threshold)
    up.df <- final[final$log2FoldChange > log2FC_threshold,]
    dwn.df <- final[final$log2FoldChange < -log2FC_threshold,]
    comb.df <- final[final$log2FoldChange < -log2FC_threshold |
                            final$log2FoldChange > log2FC_threshold,]

    results_list <- list(up.df, dwn.df, comb.df)
    results_names <- c("up", "dwn", "cmb")

    # handle organism here
    if(organism == "human"){
        org.db <- org.Hs.eg.db
    	kegg.org.id <- "hsa"
    } else if(organism == "mouse"){
        org.db <- org.Mm.eg.db
    	kegg.org.id <- "mmu"
    }
	

    # Iterate through results list
    # Run each go list for each results list
    for (i in seq(1,3)){
    	print(paste0("Starting ", fname, "_", results_names[[i]],  " GO Enrichment"))
        print(paste0("Num DEG input: ", nrow(results_list[[i]])))
        print(head(rownames(results_list[[i]])))
        print("Run GO")
        ### This trycatch will prevent incompatible symbols from throwing errors
        an.error.occured <- FALSE
        tryCatch({
            ### First convert gene symbols for input and universe
            # CP GO doesn't need it but KEGG and REACT do
            gene.df <- suppressMessages(suppressWarnings(
                bitr(rownames(results_list[[i]]),
                    fromType = "SYMBOL",
                    toType = c("UNIPROT", "ENTREZID"),
                    OrgDb = org.db)))
            print("dims uniprot df")
            print(dim(gene.df))
            print(head(gene.df))

            ### **** add kegg bitr for kegg pathway analysis
            kegg.gene.df <- suppressMessages(suppressWarnings(
                bitr_kegg(gene.df$UNIPROT,
                    fromType='uniprot',
                    toType='kegg',
                    organism=kegg.org.id)))
            print("dims kegg df")
            print(dim(kegg.gene.df))
            print(head(kegg.gene.df))

            gene.df.unv <- suppressMessages(suppressWarnings(
                bitr(universe,
                    fromType = "SYMBOL",
                    toType = c("UNIPROT", "ENTREZID"),
                    OrgDb = org.db)))
           print("dims uniprot universe df")
           print(dim(gene.df.unv))

           ### **** add kegg bitr for kegg pathway analysis
           kegg.gene.df.unv <- suppressMessages(suppressWarnings(
                bitr_kegg(gene.df.unv$UNIPROT,
                    fromType='uniprot',
                    toType='kegg',
                    organism=kegg.org.id)))
           print("dims kegg universe df")
           print(dim(kegg.gene.df.unv))

           #############  Run CP GO  ##############
           # Cellular Component
           ego_cc <- enrichGO(gene = gene.df$SYMBOL,
                    OrgDb = org.db,
                    keyType = 'SYMBOL',
                    universe = universe,
                    ont = "CC",
                    pAdjustMethod = "BH",
                    pvalueCutoff = 1,
                    qvalueCutoff = 1)

           # Biological Process
           ego_bp <- enrichGO(gene = gene.df$SYMBOL,
                    OrgDb = org.db,
                    keyType = 'SYMBOL',
                    universe = universe,
                    ont = "BP",
                    pAdjustMethod = "BH",
                    pvalueCutoff = 1,
                    qvalueCutoff = 1)

           # Molecular Function
           ego_mf <- enrichGO(gene = gene.df$SYMBOL,
                    OrgDb = org.db,
                    keyType = 'SYMBOL',
                    universe = universe,
                    ont = "MF",
                    pAdjustMethod = "BH",
                    pvalueCutoff = 1,
                    qvalueCutoff = 1)

            ##### Write GO tables
            ##### Cellular component
            ### apply simplify to reduce redundant terms
            #               ego_cc <- pairwise_termsim(ego_cc)
            #               ego_cc <- simplify(ego_cc, cutoff=0.7, by="p.adjust", select_fun=min)
            if(!is.null(ego_cc)){
                write.table(ego_cc@result, file = paste0(outdir,"csv_tables/GO/", fname, "_", results_names[[i]], "_CC_go.csv"))
                write_tableHTML(tableHTML(ego_cc@result), file =
                                    paste0(outdir, "html_tables/GO/", fname, "_", results_names[[i]], "_CC_go.html"))
                p1 <- dotplot(ego_cc, showCategory=20) + ggtitle("cellular component")
                fn1 <- paste0(fname, "_", results_names[[i]], "_CC_DTPLT.pdf")
                ggsave(filename = fn1, plot = p1, device = "pdf", path = paste0(outdir, "plots/"), height =10, width = 8)
                saveRDS(ego_cc, paste0(outdir,"cp_objects/GO/", fname, "_", results_names[[i]], "_CC_go.RDS"))
				
				cpobj_conv <- setReadable(ego_cc, org.db, 'ENTREZID')
				p2 <- cnetplot(cpobj_conv, showCategory = 5,
#                   max.overlaps = 1000,
                   cex_label_category = 1,
                   cex_label_gene = 1,
                   colorEdge = FALSE) + labs(title = str_wrap(paste0(i, " net plot of top GO Terms"), 60))
			
                fn2 <- paste0(fname, "_", results_names[[i]], "_CC_NTPLT.pdf")
                ggsave(filename = fn2, plot = p2, device = "pdf", path = paste0(outdir, "plots/"), height =10, width = 14)

				p3 <- heatplot(cpobj_conv, showCategory = 30) + labs(title = str_wrap(paste0(i, " heatmap of top GO Terms"), 60))
                fn3 <- paste0(fname, "_", results_names[[i]], "_CC_HTPLT.pdf")
                ggsave(filename = fn3, plot = p3, device = "pdf", path = paste0(outdir, "plots/"), height = 10, width = 14)

#				p4 <- treeplot(cpobj_conv, hclust_method = "average")
#                fn4 <- paste0(fname, "_", results_names[[i]], "_CC_TRPLT.pdf")
#                ggsave(filename = fn4, plot = p4, device = "pdf", path = paste0(outdir, "plots/"), height = 10, width = 8)
#

            }
            #### Biological Process
            if(!is.null(ego_bp)){
                #               ego_bp <- pairwise_termsim(ego_bp)
                #               ego_bp <- simplify(ego_bp, cutoff=0.7, by="p.adjust", select_fun=min)
                write.table(ego_bp@result, file = paste0(outdir,"csv_tables/GO/", fname, "_", results_names[[i]], "_BP_go.csv"))
                write_tableHTML(tableHTML(ego_bp@result), file =
                            paste0(outdir, "html_tables/GO/", fname, "_", results_names[[i]], "_BP_go.html"))
                p2 <- dotplot(ego_bp, showCategory=20) + ggtitle("biological process")
                fn2<- paste0(fname, "_", results_names[[i]], "_BP_DTPLT.pdf")
                ggsave(filename = fn2, plot = p2, device = "pdf", path = paste0(outdir, "plots/"),  height =10, width = 8)
                saveRDS(ego_bp, paste0(outdir,"cp_objects/GO/", fname, "_", results_names[[i]], "_BP_go.RDS"))

				cpobj_conv <- setReadable(ego_bp, org.db, 'ENTREZID')
				p2 <- cnetplot(cpobj_conv, showCategory = 5,
#                   max.overlaps = 1000,
                   cex_label_category = 1,
                   cex_label_gene = 1,
                   colorEdge = FALSE) + labs(title = str_wrap(paste0(i, " net plot of top GO Terms"), 60))
			
                fn2 <- paste0(fname, "_", results_names[[i]], "_BP_NTPLT.pdf")
                ggsave(filename = fn2, plot = p2, device = "pdf", path = paste0(outdir, "plots/"), height =10, width = 14)
    			
				p3 <- heatplot(cpobj_conv, showCategory = 30) + labs(title = str_wrap(paste0(i, " heatmap of top GO Terms"), 60))
                fn3 <- paste0(fname, "_", results_names[[i]], "_BP_HTPLT.pdf")
                ggsave(filename = fn3, plot = p3, device = "pdf", path = paste0(outdir, "plots/"), height = 10, width = 14)
    			
#				p4 <- treeplot(cpobj_conv, hclust_method = "average")
#                fn4 <- paste0(fname, "_", results_names[[i]], "_BP_TRPLT.pdf")
#                ggsave(filename = fn4, plot = p4, device = "pdf", path = paste0(outdir, "plots/"), height = 10, width = 8)
#

            }
            #### Molecular Function
            if(!is.null(ego_mf)){
                #               ego_mf <- simplify(ego_mf, cutoff=0.7, by="p.adjust", select_fun=min)
                #               ego_mf <- pairwise_termsim(ego_mf)
                write.table(ego_mf@result, file = paste0(outdir,"csv_tables/GO/",fname, "_", results_names[[i]], "_MF_go.csv"))
                write_tableHTML(tableHTML(ego_mf@result), file =
                                    paste0(outdir, "html_tables/GO/", fname, "_", results_names[[i]], "_MF_go.html"))
                p3 <- dotplot(ego_mf, showCategory=20) + ggtitle("molecular function")
                fn3 <- paste0(fname, "_", results_names[[i]], "_MF_DTPLT.pdf")
                ggsave(filename = fn3, plot = p3, device = "pdf", path = paste0(outdir, "plots/"), height = 10, width = 8)
                saveRDS(ego_mf, paste0(outdir,"cp_objects/GO/", fname, "_", results_names[[i]], "_MF_go.RDS"))

				cpobj_conv <- setReadable(ego_mf, org.db, 'ENTREZID')
				p2 <- cnetplot(cpobj_conv, showCategory = 5,
#                   max.overlaps = 1000,
                   cex_label_category = 1,
                   cex_label_gene = 1,
                   colorEdge = FALSE) + labs(title = str_wrap(paste0(i, " net plot of top GO Terms"), 60))
			
                fn2 <- paste0(fname, "_", results_names[[i]], "_MF_NTPLT.pdf")
                ggsave(filename = fn2, plot = p2, device = "pdf", path = paste0(outdir, "plots/"), height = 10, width = 14)
    			
				p3 <- heatplot(cpobj_conv, showCategory = 30) + labs(title = str_wrap(paste0(i, " heatmap of top GO Terms"), 60))
                fn3 <- paste0(fname, "_", results_names[[i]], "_MF_HTPLT.pdf")
                ggsave(filename = fn3, plot = p3, device = "pdf", path = paste0(outdir, "plots/"), height = 10, width = 14)

#				p4 <- treeplot(cpobj_conv, hclust_method = "average")
#                fn4 <- paste0(fname, "_", results_names[[i]], "_MF_TRPLT.pdf")
#                ggsave(filename = fn4, plot = p4, device = "pdf", path = paste0(outdir, "plots/"), height = 10, width = 8)
#
            }
            ############# Run KEGG and PATHWAY ##############
            print("run pathway")
            rc <- enrichPathway(gene = gene.df$ENTREZID,
                        pvalueCutoff = 0,
                        universe = gene.df.unv$ENTREZID,
                        readable = TRUE,
                        organism = organism)


            ###### Run Kegg. with kegg ids
            print("run kegg")
            kg <- enrichKEGG(gene = kegg.gene.df$kegg,
                        keyType = 'kegg',
                        organism = kegg.org.id,
                        universe = kegg.gene.df.unv$kegg,
                        pvalueCutoff = 0)

            #### Write KEGG table
            if(!is.null(kg)){
                write.table(kg@result, file = paste0(outdir,"csv_tables/KEGG/", fname, "_", results_names[[i]], "_KG_go.csv"))
                write_tableHTML(tableHTML(kg@result), file =
                            paste0(outdir, "html_tables/KEGG/", fname, "_", results_names[[i]], "_KG_go.html"))
                saveRDS(kg, paste0(outdir,"cp_objects/KEGG/", fname, "_", results_names[[i]], "_KG_go.RDS"))

            }
            #### Write PATHWAY table
            if(!is.null(rc)){
                write.table(rc@result, file = paste0(outdir,"csv_tables/REACT/", fname, "_", results_names[[i]], "_RCT_go.csv"))
                write_tableHTML(tableHTML(rc@result), file =
                            paste0(outdir, "html_tables/REACT/", fname, "_", results_names[[i]], "_RCT_go.html"))
                saveRDS(rc, paste0(outdir,"cp_objects/KEGG/", fname, "_", results_names[[i]], "_RCT_go.RDS"))

            }
        }, error = function(e) {an.error.occured <<- TRUE})
        print(paste0("Error during GO: ", an.error.occured))
    }
}	


#### Function that iterates over DE output matrices and runs cluster profiler 
run_CP <- function(path, organism = organism, log2FC_threshold = .58, padj_threshold = 0.1, run_cp){

	if (run_cp == FALSE){return(print("skipping cluster profiler"))}

	## Used to read deseq output 
	de_path <- paste0(path, "/DESeq_out/")
	## Will draw features for CP universe
	mat_path <- paste0(path, "Filtered_Matrices/")

	res_files <- list.files(de_path)
	mat_files <- list.files(mat_path)

	for (f in res_files){
	   #print(f)
	   fpre <- str_remove(string = f, pattern = (".csv")) 
	   fpre <- str_remove(string = fpre, pattern = ("DESeq2_")) 
	   #print(fpre)

	   fmat <- grep(pattern = paste0("^", fpre), x = mat_files, value = TRUE)
	   #print(fmat)
	   # read f and fmat into dataframes 
		
	   ### Reading files
      	   res <- read.csv(paste0(de_path, f), sep = ',', header=TRUE, row.names=1, 
				stringsAsFactors=FALSE)
           filt.mat <- read.csv(paste0(mat_path, fmat), sep=',', header=TRUE, row.names=1, 
				stringsAsFactors=FALSE)
	   pass.features <- rownames(filt.mat)

	   print(head(res))
	   print(length(pass.features))
	   ### apply input to cluster profiler 
 	   clusterProfiler_GO(DE.table = res, fname = fpre, universe = pass.features, log2FC_threshold = log2FC_threshold,
                             padj_threshold = padj_threshold, outdir = path, organism = organism) 

	}
}




