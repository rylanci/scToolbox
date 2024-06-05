suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(tableHTML))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(gridExtra))
suppressPackageStartupMessages(library(viridis))
suppressPackageStartupMessages(library(reshape2))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(ComplexHeatmap))
suppressPackageStartupMessages(library(ggrepel))

#Note: This module contains functions and performs actions for  post DE reporting
# As well as launch point for GO 

### Functions 
### Functions
process_args <- function(){
    # create parser object
    parser <- ArgumentParser(description = "Create seurat object from 10x h5 and collect QC metrics")

    # Create parser arg group for our required args
    required_arg_group = parser$add_argument_group('flagged required arguments',
        'the script will fail if these args are not included')

    required_arg_group$add_argument("-rp","--results_path", required = TRUE,
        help="path to root dir of Results", type="character")
    
    required_arg_group$add_argument("-mp","--matrix_path", required = TRUE,
        help="path to bulk matrices", type="character")

    required_arg_group$add_argument("-ap","--padj_threshold", required = TRUE,
        help="column containing celltype label", type="double")

    required_arg_group$add_argument("-lfc","--log2FC_threshold", required = TRUE,
        help="column containing celltype label", type="double")

    required_arg_group$add_argument("-org","--organism", required = TRUE,
        help="human or mouse", type="character")

    required_arg_group$add_argument("-pn","--project_name", required = TRUE,
        help="ex. LungBPD", type="character")

    required_arg_group$add_argument("-sp","--sobj_path", required = TRUE,
        help="human or mouse", type="character")

    required_arg_group$add_argument("-ct","--celltype_col", required = TRUE,
        help="human or mouse", type="character")

    required_arg_group$add_argument("-a","--assay", required = TRUE,
        help="RNA or ATAC", type="character")


    # optional args
    parser$add_argument("-comp","--comparisons",
        help="path to tsv containing comparisons. Use header: Target  Reference")

    parser$add_argument("-run_cp","--run_ClusterProfiler", default = TRUE,
        help="if FALSE skips clusterProfiler")


    args <- parser$parse_args()
    return(args)
}



create_vplot <- function(res.df, comparison, lfc.thresh = 1, p.thresh = 0.05){
    
#    res.df <- read.csv(res.file, sep = ',', header=TRUE, row.names=1, stringsAsFactors=FALSE)

    pass.padj.lfc <- (res.df$log2FoldChange >= lfc.thresh | res.df$log2FoldChange <= -lfc.thresh) & res.df$padj <= p.thresh
    res.df$threshold <- pass.padj.lfc

    res.df.sig <- res.df[res.df$threshold == "TRUE",]
	# Format Labels. Top up and down LFC and padj 
	top_lfc_up <-  res.df.sig[order(res.df.sig$log2FoldChange, decreasing = T),]
	top_lfc_dwn <- res.df.sig[order(res.df.sig$log2FoldChange),]
	
	top_padj_up <- top_lfc_up[order(top_lfc_up$padj),]
	top_padj_dwn <- top_lfc_dwn[order(top_lfc_dwn$padj),]

	top_up <- rbind(head(top_lfc_up, n=10), head(top_padj_up, n=10))
#	ugenes <- unique(row.names(top_up))
#	top_up <- top_up[ugenes,]

	top_dwn <- rbind(head(top_lfc_dwn,n=10), head(top_padj_dwn, n=10))
#	ugenes <- unique(row.names(top_dwn))
#	top_dwn <- top_dwn[ugenes,]


#    res.df.sig.o <- res.df.sig[order(res.df.sig$log2FoldChange),]
#    label_mtx <- rbind(head(res.df.sig.o, n=10), # grabbing top 10 up and down log2FC
#                       head(res.df.sig.o[order(-log10(res.df.sig.o$padj), decreasing = T),], n=10))

	label_mtx <- rbind(top_up, top_dwn)
    label_mtx <- unique(label_mtx)
    label_mtx$labels <- rownames(label_mtx)

    
    
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
      ggtitle(comparison) + theme_classic() + theme(aspect.ratio=1) +
      scale_colour_manual(values = c("grey", "blue")) + 
	  geom_text_repel(data=label_mtx, aes(x=log2FoldChange, y=-log10(padj), 
			label=labels),colour="black", vjust="inward",hjust="inward")
    

    
#    ggsave(filename = paste0("test_out/", comparison, ".pdf"), plot = vplot)
    return(vplot)
}


createReport <- function(res, fname, p.thresh = .05, lfc.thresh = .58, outdir){
	
	print(fname)
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
		  theme_classic() + 
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
          ggtitle(fname) + theme_classic() + theme(aspect.ratio=1) 

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
          ggtitle(fname) + theme_classic() + theme(aspect.ratio=1) 



	vplot2 <- create_vplot(res = res, comparison = fname)

	ggsave(filename = paste0(fname, "_vplot.pdf"), plot = vplot2, device = "pdf", path = paste0(outdir, "Volcano_Plots/"), 
		width = 8, height = 8)

    ## Split features by different p-value cutoffs
    pval_table <- lapply(c(0.001, 0.01, 0.025, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5), function(x) {
        res.df.t1 <- res.df[res.df[["padj"]] <= x,]
        res.df.t2 <- res.df.t1[res.df.t1[["log2FoldChange"]] >= lfc.thresh | res.df.t1[["log2FoldChange"]] <= -lfc.thresh,]
        data.frame('Cut' = x, 'Count' = nrow(res.df.t2))
    })
    pval_table <- do.call(rbind, pval_table)
    tgrob <- tableGrob(pval_table)

    ### add plots to list and return list
    plot.list <- list(maplot, vplot2, tgrob)
    return(plot.list)

}



#### Function that iterates over celltypes and runs DE report 
run_Report <- function(path, p.thresh = .05, lfc.thresh = .58){

	print("Generating Reports")

    de_path <- paste0(path, "/DESeq_out/")
    rdir <- paste0(path, "DEreports/")
    dir.create(rdir, showWarnings = FALSE)
	dir.create(paste0(rdir, "Volcano_Plots/"), showWarnings = FALSE)
    res_files <- list.files(de_path)

    deg_sum <- read.table(paste0(path, "DEG_countSummary.txt"), sep = '\t', header=TRUE)
    celltypes <- unique(deg_sum[,1]) 
    #read DEG table, extract unqiue celltypes from table 
    # iterate over celltypes, collect outfiles for each celltype 
    # iterate over outfiles. Create report for each outfile 
    # if report list has len > 1. Save report 
    for (c in celltypes){
    	print(c)
    	ct.files <- grep(pattern = paste0("^DESeq2_", c, "_"), x = res_files, value = TRUE)
#    	print(ct.files)
        plot.list <- list()
        for (f in ct.files){
            res <- read.table(paste0(de_path, f), sep = ',', header=TRUE, row.names=1)

            fpre <- str_remove(string = f, pattern = (".csv")) 
            fpre <- str_remove(string = fpre, pattern = ("DESeq2_")) 

            plots <- createReport(res = res, fname = fpre, p.thresh = as.double(p.thresh), 
                                  lfc.thresh = as.double(lfc.thresh), outdir = rdir) 

            plot.list <- append(plot.list, plots)	   
    	}

    	if(length(plot.list) > 0){
        	pdf(file = paste0(rdir,  c, "_report.pdf"), width = 16, height = length(plot.list) * 2)
                 do.call("grid.arrange", c(plot.list, ncol=3))
	     	dev.off()
    	}         
    }
}


createHeatmap <- function(res, bulk.mat, meta, conditions, q.thresh = 0.05, lfc.thresh = 1){
    print("create heatmap")

    # Extract donors per condition in metadata
    donor.meta <- meta[meta$condition %in% conditions,]
    donor.meta <- donor.meta[order(donor.meta$condition),]
    donor.mat <- bulk.mat[,colnames(bulk.mat) %in% rownames(donor.meta)]
    donor.mat.cpm <- data.frame(matrix(data = 0, nrow = nrow(donor.mat),ncol = ncol(donor.mat)))

    for (col in 1:ncol(donor.mat)) {
        csum <- sum(donor.mat[col])
        #m.cpm[,col] <- m[,col] / (csum / 10^6)
        donor.mat.cpm[,col] <- donor.mat[,col] / csum * 10^6
    }
    colnames(donor.mat.cpm) <- colnames(donor.mat)
    rownames(donor.mat.cpm) <- rownames(donor.mat)

    # Filter by .1
    dtable.f <- res[res$padj < q.thresh,]

    ### Gather our Combined features, UP features, and DWN features
    dtable.comb <- dtable.f[dtable.f$log2FoldChange > lfc.thresh | dtable.f$log2FoldChange < -lfc.thresh,]
    dtable.comb <- dtable.comb[order(dtable.comb$log2FoldChange, decreasing = TRUE), ]
	print(paste0("# deg: ", nrow(dtable.comb)))
	
	if (nrow(dtable.comb) > 1){

    	donor.mat.cpm$features <- rownames(donor.mat.cpm)

		m.cpm.p.comb <- donor.mat.cpm %>%
			slice(match(rownames(dtable.comb), features))

	    m.cpm.p.comb <- m.cpm.p.comb[,-ncol(m.cpm.p.comb)]


		### Z-score. z-score after feature selection
		m.cpm.p.z.comb <- t(scale(t(m.cpm.p.comb)))

		### Legend 
		#col_fun = colorRamp2(c(0, 0.5, 1), c("blue", "white", "red"))
		if("batch" %in% colnames(donor.meta)){
			ha.b = HeatmapAnnotation(Condition = donor.meta$condition, Batch = donor.meta$batch,
				annotation_legend_param = list(
					Batch = list(
						title = "Batch",
						at = unique(donor.meta$batch),
						labels = unique(donor.meta$batch)
						),
					Condition = list(
						title = "Condition",
						at = unique(donor.meta$condition),
						labels = unique(donor.meta$condition)
						)
			))
		} else {
			ha.b = HeatmapAnnotation(Condition = donor.meta$condition,
				annotation_legend_param = list(
					Condition = list(
						title = "Condition",
						at = unique(donor.meta$condition),
						labels = unique(donor.meta$condition)
						)
			))

		}	

		hm <- Heatmap(m.cpm.p.z.comb, cluster_rows = FALSE, cluster_columns = FALSE, col = plasma(256),
				show_row_names = FALSE, column_names_rot = 60, bottom_annotation = ha.b)

	    return(hm)


	} else {
 		print("Cannot create heatmap with 0 DEG")
 		return(NULL)
 	}		

}


runHeatmaps <- function(path, mat.path, q.thresh = 0.05, lfc.thresh = 1, proj_name){
	print("Begin Creating Heatmaps")

    rdir <- paste0(path, "Heatmaps/")
    dir.create(rdir, showWarnings = FALSE)
    de_path <- paste0(path, "/DESeq_out/")
    res_files <- list.files(de_path)
    mat_files <- list.files(mat.path)
    meta_file <- grep(pattern = "meta", x = mat_files, value = TRUE)
    meta <- read.table(paste0(mat.path, meta_file), sep = ",", row.names = 1, header = TRUE)

    deg_sum <- read.table(paste0(path, "DEG_countSummary.txt"), sep = '\t', header=TRUE)
    celltypes <- unique(deg_sum[,1])
    for (c in celltypes){
		print("------")
        print(paste0("Creating Heatmaps for: ", c))
        ct.files <- grep(pattern = paste0("^DESeq2_", c, "_"), x = res_files, value = TRUE)
        bmat_file <- grep(pattern = paste0(proj_name, "_", c, "_"), x = mat_files, value = TRUE)
        bmat <- read.table(paste0(mat.path, bmat_file), sep = ",", row.names = 1, header = TRUE,
	    		,check.names=FALSE)

         for (f in ct.files){
            res <- read.csv(paste0(de_path, f), sep = ',', header=TRUE, row.names=1,
                             stringsAsFactors=FALSE)
            # Removes everything before the first two instances of _
            # Then again in case there's more _'s 
            fpre1 <- str_remove(string = f, pattern = "[^_]*_[^_]*_")
#           fpre1 <- str_remove(string = fpre1, pattern = "[^_]*_[^_]*_")
#           fpre1 <- str_remove(string = fpre1, pattern = "[^_]*_")
            fpre2 <- str_remove(string = fpre1, pattern = (".csv"))
            conditions <- str_split(string = fpre2, pattern = "-")[[1]]
 			print("Using conditions:")
            print(conditions)

	    	#### Order heatmap by condition
		    # meta.nf <- meta.nf[order(meta.nf$condition),]
			hm <- createHeatmap(res = res, bulk.mat = bmat, meta = meta, conditions = conditions,
                                      q.thresh = as.double(q.thresh), lfc.thresh = as.double(lfc.thresh))

			if (!is.null(hm)){ 
            	pdf(file = paste0(rdir,  c, "_", fpre2, "_heatmaps.pdf"), width = 10, height = 12)
             	print(hm)
	            dev.off()
			}
        
        }
    }
}


mirror.DEG.barplot <- function(DEG_sum, seuratObj, celltype.col, comparison){
    colnames(DEG_sum) <- c("celltype", "regulation", "pairwise", "counts", "formula")

    # creating up- and down-regulated DEG dataframe
    up_down_df <- data.frame()
    ct_counts <- as.data.frame(table(seuratObj[[celltype.col]]))
    up_counts <- c()
    down_counts <- c()

    # go through each celltype, count how number of Up-/Down-DEGs, and add to up_counts/down_counts list
    # doesn't consider different comparisons, just adding all of them together
    celltypes <- sort(unique(seuratObj[[celltype.col]][,1]))
    for(i in celltypes){
        ct_down <- DEG_sum[which(DEG_sum$celltype == i & DEG_sum$regulation == "DOWN"),]
        down_counts <- c(down_counts, -sum(ct_down$counts))

        ct_up <- DEG_sum[which(DEG_sum$celltype == i & DEG_sum$regulation == "UP"),]
        up_counts <- c(up_counts, sum(ct_up$counts))
    }

    up_down_df <- data.frame(ct_counts, up_counts, down_counts)
    colnames(up_down_df) <- c("celltypes", "counts", "up_deg", "down_deg")
    #print(head(up_down_df))

    # creating basic ggplot obj and returning to user - can customize outside of function call
    up.down.Plot <- ggplot(data=up_down_df, aes(x = reorder(celltypes, -counts), y = up_deg, fill = celltypes)) +
        geom_bar(data=up_down_df, aes(x = reorder(celltypes, -counts), y = down_deg, fill = celltypes),
                 stat="identity", alpha = 0.5) +
        geom_col(position = "stack") + coord_cartesian(ylim = c(min(up_down_df$down_deg), max(up_down_df$up_deg))) +
        theme_classic() +
        geom_text(aes(label=up_deg), hjust=0.5, vjust=-0.1) +
        geom_text(data=up_down_df, aes(x = reorder(celltypes, -counts), y = down_deg, fill = celltypes,
                                       label=-down_deg), hjust=0.5, vjust=1.2) +
        xlab("Cell Types") + ylab("Up- & Down-regulated DEGs")  +
        scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
        geom_hline(yintercept = 0,colour = "black") + 
        theme(axis.text.x = element_text(size=14, angle=60, vjust = 1, hjust=1)) +
        ggtitle(paste0("n = ", sum(up_down_df$up_deg, -up_down_df$down_deg), " DEGs:    ", comparison))

    
    return(up.down.Plot)
}


runBarplot <- function(path, sobj.path, celltype.col){
    ### Create barplots for each comparison

	print("Generating barplots")
    sobj <- readRDS(sobj.path)
    dir.create(paste0(path, "SummaryBarplots/"), showWarning = FALSE)
    rdir <- paste0(path, "SummaryBarplots/")
    deg_sum <- read.table(paste0(path, "DEG_countSummary.txt"), sep = '\t')
   
    comparisons <- unique(deg_sum$V3)
    celltypes <- unique(deg_sum$V1)

    for (c in comparisons){
        deg_sum.c <- deg_sum[deg_sum$V3 == c,]
        bp <- mirror.DEG.barplot(DEG_sum= deg_sum.c, seuratObj = sobj, celltype.col = celltype.col , comparison = c)

        print(paste0(rdir, c, "_barplot.pdf"))
        pdf(file = paste0(rdir, c, "_barplot.pdf"), width = 8 + (.3 * length(celltypes)) , height = 8)
        #do.call("grid.arrange", c(plot.list, ncol = 1))
        print(bp)
        dev.off()
    }
}



### Run post DESeq processing 
## Aquire inputs
args <- process_args()


if (args$assay == "RNA"){

	source("/tscc/nfs/home/rlancione/scToolbox/DESeq2-GO/clusterProfiler.R")
	run_CP(path = args$results_path, organism = args$organism, log2FC_threshold = args$log2FC_threshold, 
    	padj_threshold = args$padj_threshold, run_cp = args$run_ClusterProfiler)
} else {
	
	source("/tscc/nfs/home/rlancione/scToolbox/DESeq2-GO/great.R")
	iterate_rGreat(path_to_res = paste0(args$results_path, "Significant_Results/"), outdir = args$results_path)
}

run_Report(path = args$results_path, lfc.thresh = args$log2FC_threshold,
           p.thresh = args$padj_threshold)


runBarplot(path = args$results_path, sobj.path = args$sobj_path, celltype.col = args$celltype_col) 


#runHeatmaps(path = args$results_path, mat.path  = args$matrix_path, q.thresh = args$padj_threshold, 
#             lfc.thresh = args$log2FC_threshold, proj_name = args$project_name)







