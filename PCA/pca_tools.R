suppressMessages(library(pcaExplorer))
suppressMessages(library(DESeq2))
suppressMessages(library(graphics))
suppressMessages(library(ggplot2))
suppressMessages(library(tidyr))
library(Seurat)
library(ggfortify)
library(datasets)
library(ggrepel)
library(stringr)
library(gridExtra)
library(pheatmap)
library(viridis)
library(ComplexHeatmap)
library(factoextra)


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


PCA <- function(sobj = sobj, bulk.by = "patient", col.by = "condition", meta.frame = meta.frame,
                       nFeatures = 500, assay = "RNA", normalization = "RLG",
                       batch.col = NULL, batch2.col = NULL, batchCorr = FALSE, 
		       outdir = NULL, name = "", coolHeatmap = TRUE){
    ### Define PCA resolution... Donor, Condition, Ect...
    res <- unique(sobj[[bulk.by]][,1])
    col.groups <- meta.frame$condition
    
    ### Bulk cells based on res
    if (assay == "RNA"){
        bdf <- data.frame(row.names = rownames(sobj@assays$RNA@counts))
        print("Using RNA Assay Data")
    } else if(assay == "SCT"){
        bdf <- data.frame(row.names = rownames(sobj@assays$SCT@data))
        print("Using SCT Assay Data")
    } else {
        print("Invalid assay specified")
        print("Choose RNA (default) or SCT")
    }
    Idents(sobj) <- bulk.by
    for (c in res){
        t.sobj <- subset(sobj, idents = c)
        # Run the condition again for t.sobj
        if (assay == "RNA"){
            rsums <- rowSums(t.sobj@assays$RNA@counts)
        } else if(assay == "SCT"){
            rsums <- rowSums(t.sobj@assays$SCT@data)
        }
        bdf <- cbind(bdf, rsums)
    }
    colnames(bdf) <- res
    
    ### Subset metadata to contain only rows found in bulk count matrix
    meta.frame <- meta.frame[colnames(bdf),]
   
    ### **** Design formula only effects results reporting by adjusting fold change.i
    ### Not necessary to have full design here because it does not impact pca. 
    ### To add covariates to correct for, use limma:removeBatchEffects 

    ### Initialize DESEq dataset with design formula 
    #ff = as.formula ("~ condition + cellquant + sex + donation_age + race + exp_batch")
    ff = as.formula("~ condition")
    dds  <- suppressWarnings(DESeqDataSetFromMatrix(round(bdf), colData = meta.frame, design = ff))
    
    ### Normalize by VST/RLG. 
    ### DESeq guy said he preffers VST cause its fast but it fails on small matrixes with not enough var features. 
    if (normalization == "RLG"){
        dds <- rlogTransformation(dds)
        print("Using Rlog Normalization")
    } else if (normalization == "VST"){
        dds <- vst(dds, blind=FALSE)
        print("Using VST Normalization")
    } else {
        print("Invalid normalization method. Use RLG (default) or VST.")
    }
    
    
    ### Batch correction
    if (batchCorr == TRUE & is.null(batch2.col) == FALSE){
        print("Using batch corrected matrix")
    	print(paste0("Batch 1:", batch.col))
	print(paste0("Batch 2:", batch2.col))

        mm <- model.matrix(~ condition, colData(dds))
        btch_mat <- limma::removeBatchEffect(assay(dds), batch=meta.frame[[batch.col]], 
				batch2 = meta.frame[[batch2.col]], design=mm)
        assay(dds) <- btch_mat
    } else if (batchCorr == TRUE){
        print("Using batch corrected matrix")
    	print(paste0("Batch: ", batch.col))
        mm <- model.matrix(~ condition, colData(dds))
        btch_mat <- limma::removeBatchEffect(assay(dds), batch=meta.frame[[batch.col]], design=mm)
        assay(dds) <- btch_mat
    }

    
    ### DESeq:PlotPCA source code 
    ### Using instead to acess PC's 3/4 & increase ggrepel overlaps. Not available with PlotPCA() 
    rv <- rowVars(assay(dds))
    #select the ntop genes by variance.
    select <- order(rv, decreasing = TRUE)[seq_len(min(nFeatures, length(rv)))]
    # perform a PCA on the data in assay(x) for the selected genes
    pca <- prcomp(t(assay(dds)[select,]))
    # aquire percent variation for each PC
    percentVar <- pca$sdev^2/sum(pca$sdev^2)
    
    cols <- c("cadetblue3", "coral3", "darkolivegreen3", "black", "mediumpurple1",
          "lightgreen", "lightgoldenrod", "lightslateblue", "mistyrose", "lightblue4",
          "navajowhite1", "magenta", "lightsalmon", "mediumorchid1", "midnightblue",
          "lightskyblue", "lightgoldenrodyellow", "black", "lightgrey", "mistyrose4", "darkcyan")

    col.groups <- meta.frame[[col.by]]
 #   col.groups <- meta.frame$condition
#    col.groups <- meta.frame$TOD

    ### Select the PCs and percentVar that you like instead of 1 and 2
    d12 <- data.frame(PC1 = pca$x[, 1], PC2 = pca$x[, 2], 
                      group = col.groups, 
                      name = rownames(dds@colData))

    d34 <- data.frame(PC3 = pca$x[, 3], PC4 = pca$x[, 4], 
                    group = col.groups, 
                    name = rownames(dds@colData))
    
    myPC12 <- ggplot(data = d12, aes_string(x = "PC1", y = "PC2", color = "group", label = "name")) + 
        geom_point(size = 3) + xlab(paste0("PC1: ", round(percentVar[1] * 100), "% variance")) + 
        ylab(paste0("PC2: ", round(percentVar[2] * 100), "% variance")) +
        coord_fixed() + theme(aspect.ratio=1) + geom_text_repel(label=rownames(dds@colData),
        max.overlaps = getOption("ggrepel.max.overlaps", default = 100))

    myPC34 <- ggplot(data = d34, aes_string(x = "PC3", y = "PC4", color = "group", label = "name")) + 
        geom_point(size = 3) + xlab(paste0("PC3: ", round(percentVar[3] * 100), "% variance")) + 
        ylab(paste0("PC4: ", round(percentVar[4] * 100), "% variance")) +
        coord_fixed() + theme(aspect.ratio=1) + geom_text_repel(label=rownames(dds@colData),
        max.overlaps = getOption("ggrepel.max.overlaps", default = 100))


    ### Correlate prcomp PCs with covariates for histograms/heatmap
    cor_df <- correlatePCs(pca ,colData(dds))

    pc1.cor <- as.data.frame(cor_df[1,])
    pc1.cor$V2 <- row.names(pc1.cor)
    colnames(pc1.cor) = c( "value","factors")
    pc1.cor <- pc1.cor %>% drop_na()
    pc1.corr.plot <- ggplot(data=pc1.cor, aes(x=factors, y=value)) +
        geom_bar(stat="identity", width=0.5) + scale_fill_hue(c = 40) +
        ggtitle("PC 1 vs covariates") + ylab("-log10(pval)") + 
        theme(axis.text.x = element_text(angle = 30, vjust = 0.65, hjust=.5))
    
    pc2.cor <- as.data.frame(cor_df[2,])
    pc2.cor$V2 <- row.names(pc2.cor)
    colnames(pc2.cor) = c( "value","factors")
    pc2.cor <- pc2.cor %>% drop_na()
    pc2.corr.plot <- ggplot(data= pc2.cor, aes(x=factors, y=value)) +
        geom_bar(stat="identity", width=0.5) + scale_fill_hue(c = 40) + 
        ggtitle("PC 2 vs covariates")+ ylab("-log10(pval)") + 
        theme(axis.text.x = element_text(angle = 30, vjust = 0.65, hjust=.5))
    
    pc3.cor <- as.data.frame(cor_df[3,])
    pc3.cor$V2 <- row.names(pc3.cor)
    colnames(pc3.cor) = c( "value","factors")
    pc3.cor <- pc3.cor %>% drop_na()
    pc3.corr.plot <- ggplot(data= pc3.cor, aes(x=factors, y=value)) +
        geom_bar(stat="identity", width=0.5) + scale_fill_hue(c = 40) + 
        ggtitle("PC 3 vs covariates")+ ylab("-log10(pval)") + 
        theme(axis.text.x = element_text(angle = 30, vjust = 0.65, hjust=.5))
    
    pc4.cor <- as.data.frame(cor_df[4,])
    pc4.cor$V2 <- row.names(pc4.cor)
    colnames(pc4.cor) = c( "value","factors")
    pc4.cor <- pc4.cor %>% drop_na()
    pc4.corr.plot <- ggplot(data= pc4.cor, aes(x=factors, y=value)) +
        geom_bar(stat="identity", width=0.5) + scale_fill_hue(c = 40) + 
        ggtitle("PC 4 vs covariates")+ ylab("-log10(pval)") + 
        theme(axis.text.x = element_text(angle = 30, vjust = 0.65, hjust=.5))
    
    ### Plot corr heatmap. Extract data from DESEq rlog matrix. Create correlation between replicates 
    bmat <- assay(dds)
    cor.mat <- cor(bmat)
   
    if(coolHeatmap){
	col_fun = colorRamp2(c(0, 0.5, 1), c("blue", "white", "red"))
	ha.t = HeatmapAnnotation(Condition = meta.frame$condition, TOD = meta.frame$TOD,
 	   annotation_legend_param = list(
                	TOD = list(
                 title = "TOD",
                 at = unique(meta.frame$TOD),
                 labels = unique(meta.frame$TOD)
             ),
		Condition = list(
                title = "Condition",
                at = unique(meta.frame$condition),
                labels = unique(meta.frame$condition)
             )

	))


#       ha.l = rowAnnotation(CellQuantity = meta.frame$cellquant, Batch = meta.frame$exp_batch,
#    	    annotation_legend_param = list(
#           	CellQuantity = list(
#                title = "CellQuantity",
#                at = c(min(meta.frame$cellquant), median(meta.frame$cellquant), max(meta.frame$cellquant)),
#                labels = c("min", "median", "max")
#            ),
#        	Batch = list(
#                title = "Batch",
#                at = unique(meta.frame$exp_batch),
#                labels = unique(meta.frame$exp_batch)
#            )
#	))
   
      heatmap <- grid.grabExpr(draw(Heatmap(cor.mat, column_names_rot = 60, col = viridis(256), bottom_annotation = ha.t)))#, left_annotation = ha.l)))
    } else {
      heatmap <- grid.grabExpr(draw(Heatmap(cor.mat, column_names_rot = 60, col = viridis(256))))
    }
   # heatmap <- grid.grabExpr(grid.draw(Heatmap(cor.mat, column_names_rot = 60, col = viridis(256))))
   # heatmap <- grid.grabExpr(draw(pheatmap(cor.mat, annotation = meta.frame[, c("condition"), drop=F])))
    
    fig.scree <- fviz_eig(pca,
              col.ind = col.groups,
              palette = cols,
              legend.title = "Groups",
              repel = TRUE
             ) #+ labs(title = paste0("PCA Cluster: ", i -1))

    ### Write Report
    pdf(file = paste0(outdir, name, "_pca_report.pdf"), width = 16, height = 22)
    grid.arrange(myPC12, myPC34, fig.scree, heatmap, pc1.corr.plot, 
	pc3.corr.plot, pc2.corr.plot, pc4.corr.plot, ncol = 2)
  #  grid.arrange(heatmap, ncol = 1)
    dev.off()

}


