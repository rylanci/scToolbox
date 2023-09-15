#library(Seurat)
#library(ggplot2)
library(clusterProfiler)
library(ReactomePA)
library(argparse)
library(stringr)
library(tableHTML)
library(ggplot2)
library(org.Hs.eg.db)
library(org.Mm.eg.db)
library(gridExtra)

#Note: This module contains functions and performs actions for  post DE reporting

### Functions 
### Functions
process_args <- function(){
    # create parser object
    parser <- ArgumentParser(description = "Create seurat object from 10x h5 and collect QC metrics")

    # Create parser arg group for our required args
    required_arg_group = parser$add_argument_group('flagged required arguments',
        'the script will fail if these args are not included')

    required_arg_group$add_argument("-rp","--results_path", required = TRUE,
        help="path to root dir of Results")

    required_arg_group$add_argument("-ap","--padj_threshold", required = TRUE,
        help="column containing celltype label")

    required_arg_group$add_argument("-lfc","--log2FC_threshold", required = TRUE,
        help="column containing celltype label")


    # optional args
    parser$add_argument("-comp","--comparisons",
        help="path to tsv containing comparisons. Use header: Target  Reference")


    args <- parser$parse_args()
    return(args)
}


createReport <- function(res, fname, p.thresh = .05, lfc.thresh = .58){

    ## Set up and determine passed features in our results for plotting
    res.df <- as.data.frame(res)
    ## Sort results by adjusted p-values
    ord <- order(res.df$padj, decreasing = FALSE)
    res.df <- res.df[ord, ]
    res.df <- cbind(data.frame(Feature = rownames(res.df)), res.df)
    rownames(res.df) <- NULL

    # Create vec of features that pass padj and lfc thresh. also rm NA
    pass.padj.lfc <- (res.df$log2FoldChange >= lfc.thresh | res.df$log2FoldChange <= -lfc.thresh) & res.df$padj <= p.thresh
    res.df$threshold <- pass.padj.lfc
    res.df <- na.omit(res.df)

    ### Create MA Plot
    maplot <- ggplot(data=res.df, aes(x=log10(baseMean), y=log2FoldChange, colour=threshold)) +
          geom_point(alpha=0.4, size=1.8) +
          geom_hline(aes(yintercept = 0), colour = "black", size = 1.2) +
    #      ylim(c(min(res.df$log2FoldChange), max(res.df$log2FoldChange))) +
          xlab("Log10 Mean expression") +
          ylab("Log2 Fold Change") +
          theme(axis.title.x = element_text(face = "bold", size = 15),
                axis.text.x = element_text(face = "bold", size = 12)) +
          theme(axis.title.y = element_text(face = "bold", size = 15),
                axis.text.y = element_text(face = "bold", size = 12)) +
          scale_colour_discrete(name = "pass thresh") +
          theme(legend.title = element_text(face = "bold", size = 15)) +
          theme(legend.text = element_text(size = 14)) +
          ggtitle(fname)

    ### Create volcano Plot
    vplot <- ggplot(data=res.df, aes(x=log2FoldChange, y=-log10(padj), col=threshold)) +
          geom_point(alpha=0.4, size=1.8) +
          geom_vline(xintercept=c(-lfc.thresh, lfc.thresh), col="black") +
          geom_hline(yintercept=-log10(p.thresh), col="black") +
#          ylim(c(0, max(res.df$log2FoldChange))) +
          xlab("Log2FoldChange") +
          ylab("-Log10(p.adj)") +
          theme(axis.title.x = element_text(face = "bold", size = 15),
                axis.text.x = element_text(face = "bold", size = 12)) +
          theme(axis.title.y = element_text(face = "bold", size = 15),
                axis.text.y = element_text(face = "bold", size = 12)) +
          scale_colour_discrete(name = "pass thresh") +
          theme(legend.title = element_text(face = "bold", size = 15)) +
          theme(legend.text = element_text(size = 14)) +
          ggtitle(fname)

    ## Split features by different p-value cutoffs
    pval_table <- lapply(c(0.001, 0.01, 0.025, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5), function(x) {
        res.df.t1 <- res.df[res.df[["padj"]] <= x,]
        res.df.t2 <- res.df.t1[res.df.t1[["log2FoldChange"]] >= lfc.thresh | res.df.t1[["log2FoldChange"]] <= -lfc.thresh,]
        data.frame('Cut' = x, 'Count' = nrow(res.df.t2))
    })
    pval_table <- do.call(rbind, pval_table)
    tgrob <- tableGrob(pval_table)

    ### add plots to list and return list
    plot.list <- list(maplot, vplot, tgrob)
    return(plot.list)

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
    # create separate directories for the various ontology platforms (GO, KEGG, Pathway)
    dir.create(paste0(outdir, "csv_tables/GO/"), showWarnings = FALSE)
    dir.create(paste0(outdir, "csv_tables/KEGG/"), showWarnings = FALSE)
    dir.create(paste0(outdir, "csv_tables/REACT/"), showWarnings = FALSE)
    dir.create(paste0(outdir, "html_tables/GO/"), showWarnings = FALSE)
    dir.create(paste0(outdir, "html_tables/KEGG/"), showWarnings = FALSE)
    dir.create(paste0(outdir, "html_tables/REACT/"), showWarnings = FALSE)
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
	if(nrow(results_list[[i]]) > 3){

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
		fn1 <- paste0(fname, "_CC_DTPLT.pdf")
		ggsave(filename = fn1, plot = p1, device = "pdf", path = paste0(outdir, "plots/"), height = 10)
	    }
	    #### Biological Process
	    if(!is.null(ego_bp)){
	        #               ego_bp <- pairwise_termsim(ego_bp)
	        #               ego_bp <- simplify(ego_bp, cutoff=0.7, by="p.adjust", select_fun=min)
	        write.table(ego_bp@result, file = paste0(outdir,"csv_tables/GO/", fname, "_", results_names[[i]], "_BP_go.csv"))
	        write_tableHTML(tableHTML(ego_bp@result), file =
			    paste0(outdir, "html_tables/GO/", fname, "_", results_names[[i]], "_BP_go.html"))
	        p2 <- dotplot(ego_bp, showCategory=20) + ggtitle("biological process")
	        fn2<- paste0(fname, "_BP_DTPLT.pdf")
	        ggsave(filename = fn2, plot = p2, device = "pdf", path = paste0(outdir, "plots/"),  height = 10)
	    }
	    #### Molecular Function
	    if(!is.null(ego_cc)){
		#               ego_mf <- simplify(ego_mf, cutoff=0.7, by="p.adjust", select_fun=min)
		#               ego_mf <- pairwise_termsim(ego_mf)
		write.table(ego_mf@result, file = paste0(outdir,"csv_tables/GO/",fname, "_", results_names[[i]], "_MF_go.csv"))
		write_tableHTML(tableHTML(ego_mf@result), file =
				    paste0(outdir, "html_tables/GO/", fname, "_", results_names[[i]], "_MF_go.html"))
		p3 <- dotplot(ego_mf, showCategory=20) + ggtitle("molecular function")
		fn3 <- paste0(fname, "_MF_DTPLT.pdf")
		ggsave(filename = fn3, plot = p3, device = "pdf", path = paste0(outdir, "plots/"), height = 10)
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
	    }
	    #### Write PATHWAY table 
	    if(!is.null(rc)){
	        write.table(rc@result, file = paste0(outdir,"csv_tables/REACT/", fname, "_", results_names[[i]], "_RCT_go.csv"))
	        write_tableHTML(tableHTML(rc@result), file =
			    paste0(outdir, "html_tables/REACT/", fname, "_", results_names[[i]], "_RCT_go.html"))
	    }
        }
    }
}



#### Function that iterates over DE output matrices and runs cluster profiler 
run_CP <- function(path, organism, log2FC_threshold = .58,
                            padj_threshold = 0.1){

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

	   fmat <- grep(pattern = fpre, x = mat_files, value = TRUE)
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



#### Function that iterates over celltypes and runs DE report 
run_Report <- function(path, p.thresh = .05, lfc.thresh = .58){

    de_path <- paste0(path, "/DESeq_out/")
    res_files <- list.files(de_path)
    deg_sum <- read.table(paste0(path, "DEG_countSummary.txt"), sep = '\t', header=TRUE)
#    print(head(deg_sum))
    celltypes <- unique(deg_sum[,1]) 
    #read DEG table, extract unqiue celltypes from table 
    # iterate over celltypes, collect outfiles for each celltype 
    # iterate over outfiles. Create report for each outfile 
    # if report list has len > 1. Save report 
    for (c in celltypes){
#	print(c)
	ct.files <- grep(pattern = paste0("^DESeq2_", c, "_"), x = res_files, value = TRUE)
#    	print(ct.files)
        plot.list <- list()
        for (f in ct.files){
	    res <- read.csv(paste0(de_path, f), sep = ',', header=TRUE, row.names=1, 
			stringsAsFactors=FALSE)
 
 	    fpre <- str_remove(string = f, pattern = (".csv")) 
	    fpre <- str_remove(string = fpre, pattern = ("DESeq2_")) 
	   
	    plots <- createReport(res = res, fname = fpre, p.thresh = as.double(p.thresh), 
			lfc.thresh = as.double(lfc.thresh)) 
            plot.list <- append(plot.list, plots)	   
	}

     rdir <- paste0(path, "DEreports/")
     dir.create(rdir, showWarnings = FALSE)
     pdf(file = paste0(rdir,  c, "_report.pdf"), width = 16, height = length(ct.files) * 4)
                 do.call("grid.arrange", c(plot.list, ncol=3))
     dev.off()
         
    }
}


### Run post DESeq processing 
## Aquire inputs
args <- process_args()


run_CP(path = args$results_path, organism = "mouse", log2FC_threshold = args$log2FC_threshold, 
       padj_threshold = args$padj_threshold)


run_Report(path = args$results_path, lfc.thresh = args$log2FC_threshold,
  	   p.thresh = args$padj_threshold)


 




