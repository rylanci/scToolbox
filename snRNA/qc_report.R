### This rscript contains fuctions to create qc reports for scRNA
### Future changes:
	# 1. Will need to add a doublet report
	# 2. Might want to restructure so there is a report for each metric.. mito, feature, doublet, ect
	# 	vs now there is a report per plot type. but ehh.c

norm_cellquant_bplot2.0 <- function (sobj, dset.col = "orig.ident", xlab = "seurat_clusters", stack.by = "condition", rand.cols = FALSE){
    ### This function needs to do 2 things
        # 1. Normalize the cell quanities by dataset.
        # 2. Plot the normized values using specified variables
    
    # Create table of cells objserved for our dataset and xaxis variable
    table <- table(sobj@meta.data[[dset.col]], sobj@meta.data[[xlab]])
    
    # Normaize this table by dataset cell quanity, so all datasets are scaled by largest. 
    # Coult use any value but largest of the group works 
    nmax <- max(rowSums(table))
    norm.df <- data.frame(row.names = rownames(table))
    for (r in seq(1, nrow(table))) {
        rsum <- sum(table[r, ])
        cfactor <- nmax/rsum
        trow <- table[r, ] * cfactor
        # the normaized quanitites are then expresed as a percent of the total accross the x-axis groups
        trow <- trow/sum(trow) * 100
        norm.df <- rbind(norm.df, trow)
    }
    
    colnames(norm.df) <- colnames(table)
    rownames(norm.df) <- rownames(table)
    
    ### Normalization is done.
        # Next step is to add our stack.by variable
        # If it is not dataset... in which case skip
    vars <- c()
    for (r in rownames(norm.df)){
      # draw the index associated with the dataset in the stack.by column
      var <-  sobj@meta.data[sobj@meta.data[[dset.col]] == r, stack.by][1]
      vars <- append(vars, var)
    }
    norm.df[[stack.by]] <- vars
    
    # melt dataframe for ggplot
    norm.df.m <- melt(norm.df, id = stack.by)
    
  #  print(head(norm.df.m))
    ### plotting 
    cols <- c("cadetblue4", "lightgoldenrod", "salmon", 
        "paleturquoise3","palegreen3", "mediumpurple1", "salmon", 
        "lightblue4", "navajowhite1", "magenta", "coral2", 
        "mediumorchid1", "midnightblue", "lightgoldenrodyellow", 
        "black", "lightgrey", "mistyrose4","darkcyan", "steelblue2", 
        "darkolivegreen3", "mediumpurple1", "lightskyblue", "firebrick2",
        "burlywood", "chartreuse1", "deeppink2", "khaki", "powderblue",
        "slategrey", "springgreen", "yellow3", "orange2", "lightsteelblue3", 
    	"tomato3", "palegreen4", "grey27", "darkseagreen", "blue", "darkorchid",
        "snow2", "peachpuff2", "magenta2", "yellowgreen", "cornflowerblue",
        "chocolate", "blueviolet", "lighblue1", "plum2")

   
    # randomize colors for fun
    if (rand.cols == TRUE){
        cols <- sample(x = cols, size = ncol(norm.df), replace = F)
    }
    
    # Create ggplot ****** note the bizzare ass method for using variables in ggplot... nice
    bp <- ggplot(norm.df.m, aes(fill = .data[[stack.by]] , y = value, x = variable)) + 
        geom_bar(position="fill", stat="identity") + scale_fill_manual(values = cols) +
        theme(axis.text.x = element_text(angle = 40, vjust = 0.5, hjust=1)) + 
        ggtitle("Normalized Stacked Barplot") + xlab(xlab) + ylab("pct.cell.quantity")
    
    
    return(bp)
 #   return(norm.df.m)

}


dset_barplot1.0 <- function(sobj, dset.col = "orig.ident", stack.by = "seurat_clusters", rand.cols = FALSE){
    
    table <- table(sobj@meta.data[[dset.col]], sobj@meta.data[[stack.by]])
    
    ### express columns as percent of a total for each dataset
    for (r in seq(1, nrow(table))) {
        trow <- table[r, ]
        table[r,] <- trow/sum(trow) * 100
    }

    
  #  print(table)
    m.table <- melt(table, id = rownames(table))
    colnames(m.table) <- c("dataset", "cluster", "pct.quantity")
    m.table$cluster <- as.factor(m.table$cluster)
    
    ### plotting 
    cols <- c("cadetblue4", "lightgoldenrod", "salmon", 
        "paleturquoise3","palegreen3", "mediumpurple1", "salmon", 
        "lightblue4", "navajowhite1", "magenta", "coral2", 
        "mediumorchid1", "midnightblue", "lightgoldenrodyellow", 
        "black", "lightgrey", "mistyrose4","darkcyan", "steelblue2", 
        "darkolivegreen3", "mediumpurple1", "lightskyblue", "firebrick2",
        "burlywood", "chartreuse1", "deeppink2", "khaki", "powderblue",
        "slategrey", "springgreen", "yellow3", "orange2", "lightsteelblue3", 
    	"tomato3", "palegreen4", "grey27", "darkseagreen", "blue", "darkorchid",
        "snow2", "peachpuff2", "magenta2", "yellowgreen", "cornflowerblue",
        "chocolate", "blueviolet", "lighblue1", "plum2")

   
    # randomize colors for fun
    if (rand.cols == TRUE){
        cols <- sample(x = cols, size = ncol(norm.df), replace = F)
    }
    
    
    bp <- ggplot(m.table, aes(fill = cluster , y = pct.quantity, x = dataset)) + 
        geom_bar(position="fill", stat="identity") + scale_fill_manual(values = cols) +
        theme(axis.text.x = element_text(angle = 40, vjust = 0.5, hjust=1)) + 
        ggtitle("Stacked Barplot") + xlab("dataset") + ylab("pct.cell.quantity")
    
    return(bp)
}


qc_report <- function(sobj, cluster.col = "seurat_clusters", condition.col = "condition", dataset.col = "orig.ident",
                     outdir){
    
    ### UMAPS 
    clust_umap <- DimPlot(sobj, group.by = cluster.col, label = TRUE, label.size = 6) +
        ggtitle("color by celltype/cluster") + theme(text = element_text(size = 20, family="ArialMT"))  
    
    cond_umap <- DimPlot(sobj, group.by = condition.col) + 
        ggtitle("color by condition") + theme(text = element_text(size = 20, family="ArialMT")) 
    
    dset_umap <- DimPlot(sobj, group.by = dataset.col) + 
        ggtitle("color by dataset") + theme(text = element_text(size = 20, family="ArialMT")) +
        theme(legend.position = "none")
    
    # save split by as independent plot. or subset
    cond_split_umap <- DimPlot(sobj, split.by = condition.col) + 
        ggtitle("") + theme(text = element_text(size = 20, family="ArialMT")) 
    
    
    ### Feature Plots
    count_fp <- FeaturePlot(sobj, features = "nCount_RNA") + 
        ggtitle("nCount_RNA") + theme(text = element_text(size = 20, family="ArialMT"))
    
    feat_fp <- FeaturePlot(sobj, features = "nFeature_RNA") + 
       ggtitle("nFeature_RNA") + theme(text = element_text(size = 20, family="ArialMT"))
   
    mito_fp <- FeaturePlot(sobj, features = "percent.mt") + 
        ggtitle("percent.mt") + theme(text = element_text(size = 20, family="ArialMT"))

    
    ### Box Plots by cluster
    count_bp_c <- ggplot(sobj@meta.data, aes(x = .data[[cluster.col]], y = nCount_RNA)) + 
      geom_boxplot() + ggtitle("nCounts per dataset") +
        theme(axis.text.x = element_text(angle = 40, vjust = 0.5, hjust=1))
    
    feat_bp_c <- ggplot(sobj@meta.data, aes(x = .data[[cluster.col]], y = nFeature_RNA)) + 
      geom_boxplot() + ggtitle("nFeatures per cluster") +
        theme(axis.text.x = element_text(angle = 40, vjust = 0.5, hjust=1))
    
    mt_bp_c <- ggplot(sobj@meta.data, aes(x = .data[[cluster.col]], y = log2(percent.mt))) + 
      geom_boxplot() + ggtitle("percent.mt per cluster") +
        theme(axis.text.x = element_text(angle = 40, vjust = 0.5, hjust=1))
    
    
    ### Box Plots by dataset
    count_bp_d <- ggplot(sobj@meta.data, aes(x = .data[[dataset.col]], y = nCount_RNA)) + 
      geom_boxplot() + ggtitle("nCounts per dataset") +
        theme(axis.text.x = element_text(angle = 20, vjust = 1, hjust=.5))
    
    feat_bp_d <- ggplot(sobj@meta.data, aes(x = .data[[dataset.col]], y = nFeature_RNA)) + 
      geom_boxplot() + ggtitle("nFeatures per dataset") +
        theme(axis.text.x = element_text(angle = 20, vjust = 0.5, hjust=.5))
    
    mt_bp_d <- ggplot(sobj@meta.data, aes(x = .data[[dataset.col]], y = log2(percent.mt))) + 
      geom_boxplot() + ggtitle("log2 percent.mt per dataset") +
        theme(axis.text.x = element_text(angle = 20, vjust = 0.5, hjust=.5))
    
    
    ## Normalized barplots (stack by condition & normalize by condition)
    ct_barplot <- norm_cellquant_bplot2.0(sobj, dset.col = condition.col, xlab = cluster.col, stack.by = condition.col)
    ## Normalized barplots (stack by dataset & normalize by dataset)
    dset_barplot <- norm_cellquant_bplot2.0(sobj, dset.col = dataset.col, xlab = cluster.col, stack.by = dataset.col)
    ## Non normalized barplot (dataset by cluster) 
    d_barplot <- dset_barplot1.0(sobj, dset.col = dataset.col, stack.by = cluster.col)
    
    ### Create a simplified report for each catagory
    dir.create(path = outdir, showWarnings = FALSE)
    ### UMAP based report 
    pdf(file = paste0(outdir, "umap_summary.pdf"), width = 22, height = 16)
        grid.arrange(clust_umap, cond_umap, dset_umap, 
                    count_fp, feat_fp, mito_fp, ncol = 3)
    dev.off()
    ### Cluster based boxplot summary 
    pdf(file = paste0(outdir, "celltype_boxplot_summary.pdf"), width = 18, height = 16)
        grid.arrange(count_bp_c, feat_bp_c, mt_bp_c, ncol = 1)
    dev.off()
    ### Dataset based 
    pdf(file = paste0(outdir, "dataset_boxplot_summary.pdf"), width = 18, height = 16)
        grid.arrange(count_bp_d, feat_bp_d, mt_bp_d, ncol = 1)
    dev.off()
    ### barplot summary 
    pdf(file = paste0(outdir, "composition_barplot_summary.pdf"), width = 18, height = 24)
        grid.arrange(ct_barplot, dset_barplot, d_barplot, ncol = 1)
    dev.off()
    
    
    ### Save individual plots at some point
    
    
}

### Changes to make 
# mandatory rasterization for umap summary and feature plots
