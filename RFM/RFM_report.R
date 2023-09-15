library(Seurat)
library(ggplot2)


### This report will display the quality of the features returned by RFM
### Each comparison has a report
# Will contain the followning plots:
#    1. FeaturePlots for each feature on the whole dataset. eh maybe later   
#    2. FeaturePlots for each feature on the celltype of interest ///
#    3. Ontologies ///
#    4. Boxplots of gene expression percell in each group ///


rfm_report <- function(sobj, ipath, opath){
	
    dir.create(path = paste0(outdir, "BoxPlots/"), showWarnings = FALSE)
    dir.create(path = paste0(outdir, "DotPlots/"), showWarnings = FALSE)
    dir.create(path = paste0(outdir, "FeaturePlots/"), showWarnings = FALSE)

    files <- list.files(ipath)
    files.top100 <- grep(pattern = "top100", x = files, value = TRUE)
    
    for (file in files.top100){
        split <- strsplit(x = file, split = "_")

        celltype <- split[[1]][1]
        cond1 <- split[[1]][2]
        cond2 <- split[[1]][3]

        top100 <- t(read.table(paste0(ipath, file), sep = ","))
        features <- c()
        for (i in top100){
           features <- append(features, i)
        }

        Idents(sobj) <- "cellclass"
        sub.sobj <- subset(sobj, idents = celltype)
        Idents(sub.sobj) <- "conditions"
        sub.sobj <- subset(sub.sobj, idents = c(cond1, cond2))
        
        #### Boxplots 
        counts <- sub.sobj@assays$RNA@counts[features,]
        counts.df <- as.data.frame(t(as.matrix(counts)))
        counts.df$conditions <- sub.sobj@meta.data$conditions
        
        bxplot <- function(feature){
            bxplt <- ggplot(counts.df, aes(x=conditions, y=log10(.data[[feature]]))) + 
                  geom_boxplot() + ggtitle(feature)
        return(bxplt)
        }

        plots <- lapply(X = features[1:25], FUN = bxplot)
        pdf(file = paste0(outdir, "BoxPlots/", celltype, "_", cond1, "_", cond2, "_bxplts.pdf"), width = 12, height = 24)
            do.call('grid.arrange',c(plots, ncol = 5))
        dev.off()
        
        #### DotPlots & FeaturePlots
        dp <- DotPlot(sub.sobj, features = features[1:30], group.by = "conditions") + theme(axis.text.x = element_text(angle = 60, hjust = 1))
        fp <- suppressMessages(FeaturePlot(sub.sobj, features = features[1:12], raster = TRUE))


        suppressMessages(ggsave(filename = paste0(celltype, "_", cond1, "_", cond2, "_dotplot.pdf"), plot = dp, device = "pdf",
                    path = paste0(outdir, "DotPlots/"), width = 9, height = 4))
        suppressMessages(ggsave(filename = paste0(celltype, "_", cond1, "_", cond2, "_featureplot.pdf"), plot = fp, device = "pdf",
                    path = paste0(outdir, "FeaturePlots/"), width = 18, height = 18))

        
    }
}

sobj <- readRDS("/projects/ps-epigen/users/rlan/Liver/manuscript_RNA/Objects/filt_sobj_anno.RDS")
path <- "/projects/ps-epigen/users/rlan/Liver/RFM_Liver/cellclass_out/model/"
outdir <- "/projects/ps-epigen/users/rlan/Liver/RFM_Liver/cellclass_out/reports/"

rfm_report(sobj = sobj, ipath = path, opath = outdir)





