library(Seurat)
library(stats)
library(gprofiler2)

### This is the binned permutation version


### version 2. Remove abs value 
sen_sig_score_rna_V2 <- function(seu_obj, gene_list) {

    DefaultAssay(seu_obj) <- 'RNA'
    # if RNA assay is unprocessed, run for sensig
    if (sum(seu_obj@assays$RNA@counts) == sum(seu_obj@assays$RNA@data)) {
        seu_obj <- NormalizeData(seu_obj)
    }
    
    # the scale data assay will change with every gene list
    seu_obj <- suppressMessages(ScaleData(seu_obj, features = gene_list))

   # print("Calculating SenSig score")
    cts_sct <- seu_obj@assays$RNA@scale.data
    deg_cts <- subset(cts_sct, rownames(cts_sct) %in% gene_list)

    scores <- colSums(deg_cts)
    scores <- trunc(scores*10^2)/10^2
   # seu_obj[[paste0(gl_col_name, "_SenSig")]] <- scores
   # return(seu_obj)
    return(scores)
    
}

sen_sig_gene <- function(seu_obj, gene_list) {

    DefaultAssay(seu_obj) <- 'RNA'
    # if RNA assay is unprocessed, run for sensig
    if (sum(seu_obj@assays$RNA@counts) == sum(seu_obj@assays$RNA@data)) {
        seu_obj <- NormalizeData(seu_obj)
    }
    # the scale data assay will change with every gene list
    seu_obj <- suppressMessages(ScaleData(seu_obj, features = gene_list))

    # print("Calculating SenSig score")
    cts_sct <- seu_obj@assays$RNA@scale.data
    deg_cts <- subset(cts_sct, rownames(cts_sct) %in% gene_list)
    
    means <- apply(deg_cts, MARGIN = 1, FUN = mean)
    medians <- apply(deg_cts, MARGIN = 1, FUN = median)
    
    m.frame <- data.frame(row.names = row.names(deg_cts), 
                          "mean" = means,
                          "median" = medians)
    return(m.frame)
}

# All mf is our data frame containing all our features as rows. And the median scaled score for each gene
permute <- function(all.mf, bin.table){
    g1.sample <- sample(x = row.names(all.mf[all.mf$bin == 1,]), size = bin.table[1], replace = FALSE)
    g2.sample <- sample(x = row.names(all.mf[all.mf$bin == 2,]), size = bin.table[2], replace = FALSE)
    g3.sample <- sample(x = row.names(all.mf[all.mf$bin == 3,]), size = bin.table[3], replace = FALSE)
    g4.sample <- sample(x = row.names(all.mf[all.mf$bin == 4,]), size = bin.table[4], replace = FALSE)
    g.sample <- c(g1.sample, g2.sample, g3.sample, g4.sample)

    p.scores <- sen_sig_score_rna_V2(sobj, g.sample)
  }


sensig_binned_permutation <- function(sobj, n_permutations, sig.thresh, score.col, name, gene_list){
    # First get the median summed score for each gene in our list	
    sig.mf <- sen_sig_gene(sobj, gene_list)
    # Extract the quantile values
    sig.min <- min(sig.mf$median)
    sig.25 <- as.double(quantile(sig.mf$median, probs = .25))
    sig.50 <- median(sig.mf$median)
    sig.75 <- as.double(quantile(sig.mf$median, probs = .75))
    sig.max <- max(sig.mf$median)
    
    # Count the sampled present in each quantile 
    # Would be able to divide by 4 in even gene samples 
    # But in odd there is an extra gene in one of the bins. So do this
    bins <- c()
    for (s in sig.mf$median){
	if(s >= sig.min & s <= sig.25){
	    bins <- append(bins, 1)
	} else if (s > sig.25 & s <= sig.50){
	    bins <- append(bins, 2)
	} else if (s > sig.50 & s <= sig.75){
	    bins <- append(bins, 3)
	} else if (s > sig.75 & s <= sig.max){
	    bins <- append(bins, 4)
	}
    }
    
    # Info will be used later when taking permuted samples 
    sig.mf$bin <- bins
    bin.table <- table(sig.mf$bin)

    # Using the quantile values. Assign all genes into a bin 
    all.mf <- sen_sig_gene(sobj, row.names(sobj@assays$RNA@counts))

    bins <- c()
    for (s in all.mf$median){
        if (s < sig.min){
            bins <- append(bins, 0)
        } else if (s >= sig.min & s <= sig.25){
            bins <- append(bins, 1)
        } else if (s > sig.25 & s <= sig.50){
            bins <- append(bins, 2)
        } else if (s > sig.50 & s <= sig.75){
            bins <- append(bins, 3)
        } else if (s > sig.75 & s <= sig.max){
            bins <- append(bins, 4)
        } else if (s > sig.max){
            bins <- append(bins, 5)
        }
    }
    all.mf$bin <- bins




    # collect permuted scores for each cell
    all.genes <- row.names(sobj@assays$RNA@counts)
    perm.df <- data.frame(row.names = row.names(sobj@meta.data))
    for (i in seq(1, n_permutations)){
	p.scores <- permute(all.mf = all.mf, bin.table= bin.table)
	perm.df[paste0("p_", i)] <- p.scores
    }

    ## Calculate pvalues
    pvalues <- c()
    scores <- sobj[[score.col]][,1]
    for (r in seq(1,nrow(perm.df))){
	# perm row pertains to a cell in perm.df
	perm.row <- perm.df[r,]

	# score pertains to the sensig score of that cell
	score <- scores[r]
	obs <- 0
	for (p in perm.row){if (p > score){obs <- obs + 1}}

	pvalues <- append(pvalues, obs/n_permutations)
    }

    psig <- c()
    for (i in pvalues){
	if (i < sig.thresh){psig <- append(psig, "senescent")
	} else {psig <- append(psig, "normal")}
    }
   
    
    ### Multiple test correction
    qvalues <- p.adjust(pvalues, method = "BH", n = length(pvalues))
	
    qsig <- c()
    for (i in qvalues){
	if (i < sig.thresh){qsig <- append(qsig, "senescent")
	} else {qsig <- append(qsig, "normal")}
    }


    sobj[[paste0(name, "_pval")]] <- pvalues
    sobj[[paste0(name, "_padj")]] <- qvalues
    sobj[[paste0(name, "_", sig.thresh, "_psig")]] <- psig
    sobj[[paste0(name, "_", sig.thresh, "_qsig")]] <- qsig

    return(sobj)

}


sensig_report <- function(sobj, name, sig.thresh, score.col, condi.col, sensig_genes, file){
    
    
    # FeaturePlot of scores 
    feat.s <- FeaturePlot(sobj, features = paste0(score.col))
    # DimPlot of conditions 
    dim.c <- DimPlot(sobj, group.by = paste0(condi.col), shuffle = TRUE)
    
    #file.create(file, showWarnings = FALSE)
    sensig.col.pval <- paste0(name, "_pval")
    sensig.col.qval <- paste0(name, "_padj")
    sensig.col.psig <- paste0(name, "_", sig.thresh, "_psig")
    sensig.col.qsig <- paste0(name, "_", sig.thresh, "_qsig")
    
    # FeaturePlots of stat values 
    ### Adding 10^-7 so no NA for 0 
    sobj[[paste0(sensig.col.pval, "_log10")]] <- -log10(sobj[[sensig.col.pval]] + 10^-4) + 10^-4 
    sobj[[paste0(sensig.col.qval, "_log10")]] <- -log10(sobj[[sensig.col.qval]] + 10^-4) + 10^-4
    
   # print(colnames(sobj@meta.data))
    
    feat.p <- FeaturePlot(sobj, features = paste0(sensig.col.pval, "_log10"), 
                          cols = c("lightgrey", "blue"), order = TRUE)
    feat.q <- FeaturePlot(sobj, features = paste0(sensig.col.qval, "_log10"), 
                          cols = c("lightgrey", "blue"), order = TRUE)
    
    
    # DimPlots of sig values 
    dim.p <- DimPlot(sobj, group.by = sensig.col.psig, cols = c("lightgrey", "cyan4"), 
                     order = c("senescent", "normal"))
    dim.q <- DimPlot(sobj, group.by = sensig.col.qsig, cols = c("lightgrey", "cyan4"), 
                     order = c("senescent", "normal"))
    
    # Volcano of score by p
    ####### Shoot columns are hardcoded. Fix. #########
    vplot.p <- ggplot(data=sobj@meta.data, aes(x=sensig_150_v2, y=,.data[[paste0(sensig.col.pval, "_log10")]], col=.data[[sensig.col.psig]])) +
          geom_point(alpha=0.4, size=1.8) +
       #   geom_vline(xintercept=c(-lfc.thresh, lfc.thresh), col="black") +
          geom_hline(yintercept=-log10(0.05), col="black") +
          ylim(c(0, max(sobj[[paste0(sensig.col.pval, "_log10")]][,1]))) +
          xlab("sensig score") +
          ylab("-Log10(pval)") +
          theme(axis.title.x = element_text(face = "bold", size = 15),
                axis.text.x = element_text(face = "bold", size = 12)) +
          theme(axis.title.y = element_text(face = "bold", size = 15),
                axis.text.y = element_text(face = "bold", size = 12)) +
          scale_colour_discrete(name = "status") +
          theme(legend.title = element_text(face = "bold", size = 15)) +
          theme(legend.text = element_text(size = 14)) +
          ggtitle("score by pval") + theme_light(aspect.ratio = 1)
    # Volcano of score by q 
    vplot.q <- ggplot(data=sobj@meta.data, aes(x=sensig_150_v2, y=.data[[paste0(sensig.col.qval, "_log10")]], col=.data[[sensig.col.qsig]])) +
          geom_point(alpha=0.4, size=1.8) +
       #   geom_vline(xintercept=c(-lfc.thresh, lfc.thresh), col="black") +
          geom_hline(yintercept=-log10(0.05), col="black") +
          ylim(c(0, max(sobj[[paste0(sensig.col.pval, "_log10")]][,1]) + .5)) +
          xlab("sensig score") +
          ylab("-Log10(qval)") +
          theme(axis.title.x = element_text(face = "bold", size = 15),
                axis.text.x = element_text(face = "bold", size = 12)) +
          theme(axis.title.y = element_text(face = "bold", size = 15),
                axis.text.y = element_text(face = "bold", size = 12)) +
          scale_colour_discrete(name = "status") +
          theme(legend.title = element_text(face = "bold", size = 15)) +
          theme(legend.text = element_text(size = 14)) +
          ggtitle("score by qval") + theme_light(aspect.ratio = 1)
    
    
    # scater plots. Score by count colored by sig 
    splot.p <- ggplot(sobj@meta.data %>% arrange(sobj[[sensig.col.psig]]), 
                      aes(x= sensig_150_v2, y=nCount_RNA, color=.data[[sensig.col.psig]])) + 
        geom_point() + theme_light(aspect.ratio = 1)


    splot.q <- ggplot(sobj@meta.data %>% arrange(sobj[[sensig.col.qsig]]), 
                      aes(x= sensig_150_v2, y=nCount_RNA, color=.data[[sensig.col.qsig]])) + 
        geom_point() + theme_light(aspect.ratio = 1)

    
    #print(max(sobj[[paste0(sensig.col.pval, "_log10")]]))
    # Rank and DotPlot of top genes by mean sensig score
    sobj <- suppressMessages(ScaleData(sobj, features = sensig_genes))
    sen.sobj <- subset(sobj, subset = sensig_150_bin_0.05_qsig == "senescent")
    
    avg_zscore_sen_cells <- sort(rowMeans(sen.sobj@assays$RNA@scale.data), decreasing = T)
    avg_zscore_sen_cells_df  <- data.frame(avg_zscore_sen_cells)
    
    colnames(avg_zscore_sen_cells_df) <- c("mean_score")
    avg_zscore_sen_cells_df$gene <- rownames(avg_zscore_sen_cells_df)
    
    adf <- avg_zscore_sen_cells_df[1:25,]
    adf <- adf %>% mutate(gene = fct_reorder(gene, mean_score))
    
    
    rankplot <- ggplot(adf, aes(x=gene, y=mean_score)) + 
      geom_bar(stat = "identity") + coord_flip() + ggtitle("Top 25 sensig genes (padj)") +
      theme_light(aspect.ratio = 1)

    
    dotplot <- suppressWarnings(DotPlot(sobj, features = rev(adf$gene), 
                    group.by = "sensig_150_bin_0.05_qsig", cols = c("lightgrey", "cyan4"))) + 
                    coord_flip() + theme_light(aspect.ratio = 1)
    
    
    pdf(file = file, width = 24, height = 16)
    grid.arrange(feat.s, dim.c, feat.p, feat.q, dim.p, dim.q, vplot.p, vplot.q, splot.p, splot.q,
                rankplot, dotplot, ncol = 4 )
    dev.off()
}





