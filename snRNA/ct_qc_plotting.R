library(Seurat)
library(ggplot2)
library(gridExtra)


### For each celltype we want boxplots by donor/replicate and conditionc 
sobj <- readRDS("/projects/ps-epigen/users/rlan/Liver/manuscript_RNA/Objects/filt_sobj_anno.RDS")
dataset.col <- "orig.ident"
condition.col <- "conditions"
celltype.col <- "cellclass"
outdir <- "/projects/ps-epigen/users/rlan/Liver/manuscript_RNA/QC/celltype_qc/"

celltypes <- unique(sobj[[celltype.col]][,1])
Idents(sobj) <- celltype.col
for (c in celltypes){
    c.sobj <- subset(sobj, idents = c)

    
    count_bp_d <- ggplot(c.sobj@meta.data, aes(x = .data[[dataset.col]], y = nCount_RNA)) +
 	geom_boxplot() + ggtitle("nCounts per dataset") +
        theme(axis.text.x = element_text(angle = 20, vjust = 1, hjust=.5))

    feat_bp_d <- ggplot(c.sobj@meta.data, aes(x = .data[[dataset.col]], y = nFeature_RNA)) +
 	geom_boxplot() + ggtitle("nFeatures per dataset") +
 	theme(axis.text.x = element_text(angle = 20, vjust = 0.5, hjust=.5))

    lcount_bp_d <- ggplot(c.sobj@meta.data, aes(x = .data[[dataset.col]], y = log10(nCount_RNA))) +
 	geom_boxplot() + ggtitle("log10 nCounts per dataset") +
        theme(axis.text.x = element_text(angle = 20, vjust = 1, hjust=.5))

    lfeat_bp_d <- ggplot(c.sobj@meta.data, aes(x = .data[[dataset.col]], y = log10(nFeature_RNA))) +
 	geom_boxplot() + ggtitle("log10 nFeatures per dataset") +
 	theme(axis.text.x = element_text(angle = 20, vjust = 0.5, hjust=.5))

    mt_bp_d <- ggplot(c.sobj@meta.data, aes(x = .data[[dataset.col]], y = log2(percent.mt))) +
 	geom_boxplot() + ggtitle("log2 percent.mt per dataset") +
 	theme(axis.text.x = element_text(angle = 20, vjust = 0.5, hjust=.5))

    ### Plots by condition
    count_bp_c <- ggplot(c.sobj@meta.data, aes(x = .data[[condition.col]], y = nCount_RNA)) +
 	geom_boxplot() + ggtitle("nCounts per dataset") +
        theme(axis.text.x = element_text(angle = 20, vjust = 1, hjust=.5))

    feat_bp_c <- ggplot(c.sobj@meta.data, aes(x = .data[[condition.col]], y = nFeature_RNA)) +
 	geom_boxplot() + ggtitle("nFeatures per dataset") +
 	theme(axis.text.x = element_text(angle = 20, vjust = 0.5, hjust=.5))

    lcount_bp_c <- ggplot(c.sobj@meta.data, aes(x = .data[[condition.col]], y = log10(nCount_RNA))) +
 	geom_boxplot() + ggtitle("log10 nCounts per dataset") +
        theme(axis.text.x = element_text(angle = 20, vjust = 1, hjust=.5))

    lfeat_bp_c <- ggplot(c.sobj@meta.data, aes(x = .data[[condition.col]], y = log10(nFeature_RNA))) +
 	geom_boxplot() + ggtitle("log10 nFeatures per dataset") +
 	theme(axis.text.x = element_text(angle = 20, vjust = 0.5, hjust=.5))

    mt_bp_c <- ggplot(c.sobj@meta.data, aes(x = .data[[condition.col]], y = log2(percent.mt))) +
 	geom_boxplot() + ggtitle("log2 percent.mt per dataset") +
 	theme(axis.text.x = element_text(angle = 20, vjust = 0.5, hjust=.5))

    # saving plots
    pdf(file = paste0(outdir, c, "_dataset_qc.pdf"), width = 18, height = 24)
        grid.arrange(count_bp_d, feat_bp_d, lcount_bp_d, lfeat_bp_d, mt_bp_d, ncol = 1)
    dev.off()

    pdf(file = paste0(outdir, c, "_condition_qc.pdf"), width = 18, height = 24)
        grid.arrange(count_bp_c, feat_bp_c, lcount_bp_c, lfeat_bp_c, mt_bp_c, ncol = 1)
    dev.off()
}


