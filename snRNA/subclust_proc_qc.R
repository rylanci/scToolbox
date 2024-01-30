library(Seurat)
library(harmony)

sobj <- ""
out <- ""
celltype <- ""
name <- ""
harmony.group <- ""
res <- ""
max.dim <- ""
comparison.table <- ""


process_args <- function(){
    # create parser object
    parser <- ArgumentParser(description = "Create seurat object from 10x h5 and collect QC metrics")

    # Create parser arg group for our required args
    required_arg_group = parser$add_argument_group('flagged required arguments',
        'the script will fail if these args are not included')
    
    ### start args 
    required_arg_group$add_argument("-sobj","--sobj_path", required = TRUE,
        help="path to seurat object")

    required_arg_group$add_argument("-o","--outdir", required = TRUE,
        help="path to output location")

    required_arg_group$add_argument("-n","--name", required = TRUE,
        help="name of analysis e.g. AT2_batch_harmony")

    required_arg_group$add_argument("-comp_table","--comp_table_path", required = TRUE,
        help="path to comparison table")


    # optional args
    parser$add_argument("-res","--resolution",
        help="resolution to use for clustering")

    parser$add_argument("-dims","--max_dims",
        help="maximum dimension to use for clustering")

    args <- parser$parse_args()
    return(args)
}


processing <- function (sobj, dims = 1:20, res = 0.05, n.neigh = 30L, min.dist = 0.3, 
    spread = 1, harmony.group = "orig.ident"){
    sobj <- RunPCA(sobj, verbose = FALSE)
    
    if (harmony.group) {
        sobj <- RunHarmony(object = sobj, reduction = "pca", 
            group.by.vars = harm.group, assay.use = "SCT", 
            project.dim = FALSE)
        
        sobj <- FindNeighbors(sobj, dims = dims, verbose = FALSE, 
            reduction = "harmony") %>% FindClusters(resolution = res, 
            verbose = FALSE) %>% RunUMAP(dims = dims, n.neighbors = n.neigh, 
            min.dist = min.dist, spread = spread, verbose = FALSE, 
            reduction = "harmony")
    }
    else {
        sobj <- FindNeighbors(sobj, dims = dims, verbose = FALSE, 
            reduction = "pca") %>% FindClusters(resolution = res, 
            verbose = FALSE) %>% RunUMAP(dims = dims, n.neighbors = n.neigh, 
            min.dist = min.dist, spread = spread, verbose = FALSE, 
            reduction = "pca")
    }
    return(sobj)
}

comp_umap <- function(sobj, comparison, cond_col, anno_col){
    cond1 <- comparison[[1]]
    cond2 <- comparison[[2]]
    
    # Divide total cells in cond1 by total cells in cond2 
    base.ratio <- table(sobj[[cond_col]])[cond1] / table(sobj[[cond_col]])[cond2]
    # table of anno col * cond col
    comp.table <- table(sobj@meta.data[,anno_col], sobj@meta.data[,cond_col])
    # extract cond1 and cond2 columns as a df 
    comp.df <- data.frame(cond1 = comp.table[,cond1], cond2 = comp.table[,cond2])
    
    ### Results in obseved over expected value for each cluster 
    comp.df$ratio <- comp.df[,"cond1"] /  comp.df[,"cond2"]
    comp.df$ratio.norm <- comp.df[,"ratio"] / base.ratio
    comp.df$foldChange <- log2(comp.df$ratio.norm)
    
    ### Add the ratio to each cell in the object pertaining to that cluster 
    ratio.list <- c()
    for (c in sobj@meta.data[[anno_col]]){
      #  print(c)
        t.row <- comp.df[rownames(comp.df) == c,]
        ratio.list <- append(ratio.list, t.row[,"foldChange"])

    }
    sobj[[paste0(cond1, "_", cond2, "_foldChange")]] <- ratio.list

    FeaturePlot(object = sobj, features = paste0(cond1, "_", cond2, "_foldChange")) +
        scale_color_gradient2(midpoint=0, low="blue", mid="cornsilk", high="red", space ="Lab", name = "FoldChange") +
        ggtitle(paste0(cond1))
    

    #ggsave(filename = "noharm_comp_umaps.pdf", plot = cplot1, device = "pdf", path = "composition_umap/", width = 18, height = 9)

}


####### Processing 
p.sobj <- processing(sobj, harm.group = "exp_batch", res = 1, dims = 1:50)
saveRDS(p.sobj, paste0(sobj.out, celltype, "_bharm.RDS")) 


############ QC reports 
### no harmony 
source("~/scToolbox/snRNA/qc_report.R")
qc_report(sobj = nh.sobj, cluster.col = "seurat_clusters", condition.col = "condition", dataset.col = "orig.ident", outdir = "QC_output/no_harmony_qc_out/")


######## Compositon UMAP
# Read comparison table
# create list of comparisons 
# apply comp umap to list
#active.fp.nh <- comp_umap(sobj = nh.sobj, comparison = active, cond_col = "condition", anno_col = "seurat_clusters")
#ggsave(filename = "noharm_comp_umaps.pdf", plot = cplot1, device = "pdf", path = "composition_umap/", width = 18, height = 9)










