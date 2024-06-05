library(Seurat)
library(Signac)
library(harmony)
library(argparse)
library(ggplot2)
library(reshape2)
library(gridExtra)

# Functions
process_args <- function(){
    # create parser object
    parser <- ArgumentParser(description = "Create seurat object from 10x h5 and collect QC metrics")

    # Create parser arg group for our required args
    required_arg_group = parser$add_argument_group('flagged required arguments',
        'the script will fail if these args are not included')
    
    ### start args 
    required_arg_group$add_argument("-sobj","--sobj_path", required = TRUE,
        help="path to seurat object")

    required_arg_group$add_argument("-o","--out_dir", required = TRUE,
        help="path to output location")

    required_arg_group$add_argument("-n","--name", required = TRUE,
        help="name of analysis e.g. AT2_batch_harmony")

#    required_arg_group$add_argument("-comp_table","--comp_table_path", required = TRUE,
#        help="path to comparison table")

    required_arg_group$add_argument("-cond","--condition_column", required = TRUE,
        help="column containing required conditions")

    required_arg_group$add_argument("-ds","--dataset_column", required = TRUE,
        help="column containing celltype to subset")




    # optional args
    parser$add_argument("-c","--celltype", 
        help="celltype to subcluster")

    parser$add_argument("-res","--resolution",
        help="resolution to use for clustering")

    parser$add_argument("-dims","--max_dims",
        help="maximum dimension to use for clustering")

    parser$add_argument("-p","--processing", default = TRUE, type = "logical", 
        help="TRUE or FALSE. Whether or not to reprocess or use existing clusters")

    parser$add_argument("-clust","--cluster_column", default = "seurat_clusters",
        help="If processing is FALSE one may want to specify an exisiting cluster column")

    parser$add_argument("-ct","--celltype_column",
        help="column containing celltype to subset")

    parser$add_argument("-hg","--harmony_group", default = "NONE",
        help="Set to ident to correct for e.g orig.ident or batch")

    parser$add_argument("-ur","--umap_reduction", default = "NONE",
        help="Set to ident to correct for e.g orig.ident or batch")


   args <- parser$parse_args()
    return(args)
}


processing <- function (sobj, dims = 1:20, res = 0.05, n.neigh = 30L, min.dist = 0.3, 
    spread = 1, harmony.group = NULL){
    sobj <- RunPCA(sobj, verbose = FALSE)
    
    if (harmony.group != "NONE") {
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



#### Start Subclustering and QC 
args <- process_args()
sobj <- readRDS(args$sobj_path)


### Record inputs
print(paste0("sobj path: ", args$sobj_path))
print(sobj)
print(paste0("outdir: ", args$out_dir))
print(paste0("proj name : ", args$name))
print(paste0("processing: ", args$processing))
print(paste0("harm.group: ", args$harm_group))


### Subcluster
if (is.null(args$celltype)){
	print("Skipping Subset")
	sobj.sub <- sobj
} else {
	print("Subset Object")
	sobj.sub <- subset(sobj, subset = args$celltype_column == args$celltype)
}

	

####### Processing 
if (args$processing){
	print("Processing")
	p.sobj <- processing(sobj.sub, harm.group = args$harmony_group, res = 1, dims = 1:args$max_dim)
	saveRDS(p.sobj, paste0(args$out_dir, celltype, "_", args$name, ".rds")) 
} else { 
	print("Skipping Processing")
	p.sobj <- sobj.sub 
}

############ QC reports 
source("~/scToolbox/snRNA/qc_report.R")
print("Starting QC Report")
qc_report(sobj = p.sobj, cluster.col = args$cluster_column, condition.col = args$condition_column, dataset.col = args$dataset_column,
	 outdir = args$out_dir, reduction = args$umap_reduction)


######## Compositon UMAP
# Read comparison table
# create list of comparisons 
# apply comp umap to list
#active.fp.nh <- comp_umap(sobj = nh.sobj, comparison = active, cond_col = "condition", anno_col = "seurat_clusters")
#ggsave(filename = "noharm_comp_umaps.pdf", plot = cplot1, device = "pdf", path = "composition_umap/", width = 18, height = 9)










