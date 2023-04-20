suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(Nebulosa))
suppressPackageStartupMessages(library(patchwork))
suppressPackageStartupMessages(library(viridis))
suppressPackageStartupMessages(library(future))
suppressPackageStartupMessages(library(harmony))
suppressPackageStartupMessages(library(preprocessCore))
suppressPackageStartupMessages(library(gridExtra))
suppressPackageStartupMessages(library(argparse))
source("~/scripts/utils.R")


######## Functions #########
process_args <- function(){
    parser <- ArgumentParser(description = "Create merged seurat object and create QC report")
   
    parser$add_argument("-ip2", "--input_path2",
        help="optional second path to dir of seurat objects for merging")
    
    # These args belong to my group of required flagged arguments
    required_arg_group = parser$add_argument_group('flagged required arguments', 
        'the script will fail if these args are not included')

    required_arg_group$add_argument("-ip", "--input_path", required = TRUE,
        help="path to dir of seurat objects for merging")

    required_arg_group$add_argument("-pr","--project_name", required = TRUE,
        help="Project data belongs too. e.g. Lung_BPD") 

    required_arg_group$add_argument("-s","--species", required = TRUE,
        help="species of origin for the dataset (mm10 or hg38)")

    required_arg_group$add_argument("-o","--output_path", required = TRUE,
        help="output path for plots and objects")

    args <- parser$parse_args()
    return(args)
}

read_objects <- function(ipath1, ipath2 = NULL){
    objs <- list.files(ipath1)
    obj_list <- c()
    for (o in objs) {
        temp_obj <- readRDS(paste0(ipath1, o))
        obj_list <- append(obj_list, temp_obj)
    }
    if (is.null(ipath2) == FALSE) {
        objs2 <- list.files(ipath2)
        for (o in objs2) {
            temp_obj2 <- readRDS(paste0(ipath2, o))
            obj_list <- append(obj_list, temp_obj2)
        }
    }
    return(obj_list)
}

adj_doub_colnames <- function(obj_list){
    obj_list2 <- c()
    for (sobj in obj_list){
        colnames(sobj@meta.data)[grep(pattern = "pANN", x = colnames(sobj@meta.data))] <- "DF.scores"
        colnames(sobj@meta.data)[grep(pattern = "DF.class", x = colnames(sobj@meta.data))] <- "DF.classifications"
    
	obj_list2 <- append(obj_list2, sobj)
    }
    return(obj_list2)
}

merging <- function(obj_list, project){
    print("Extracting sample names")
    samp_nms <- c()
    for (o in obj_list){
        samp_nm <- as.character(o@meta.data$orig.ident[1])
        samp_nms <- append(samp_nms, samp_nm)                        
    }
    
    print("Merging")
    m.sobj <- merge(x = obj_list[[1]], y = obj_list[2:length(obj_list)],
                add.cell.ids = samp_nms,
                project = project)
    
    return(m.sobj)
}

sct.harm.processing.m <- function (sobj, dims = 1:20, res = 0.05, n.neigh = 30L, min.dist = 0.3, 
    spread = 1, harmony = FALSE, CM = FALSE, sct.skip = FALSE) {
    if (sct.skip ==FALSE){
        sobj <- SCTransform(sobj, vst.flavor = "v2", conserve.memory = CM, 
            verbose = FALSE) %>% RunPCA(verbose = FALSE)
        }
    if (harmony) {
        sobj <- RunHarmony(object = sobj, reduction = "pca", 
            group.by.vars = "orig.ident", assay.use = "SCT", 
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

plotting <- function(sobj, outdir){
    print("Plotting")
    dir.create(outdir, showWarnings = FALSE)
    dir.create(paste0(outdir, "merge_plots"), showWarnings = FALSE)
    outdir <- paste0(outdir, "merge_plots/")
  
    Idents(sobj) <- "seurat_clusters"
    umap <- DimPlot(sobj, label = T, label.size = 6) + 
            ggtitle(paste0("Merged SCT Embedding   nCells = ", nrow(sobj@meta.data)))
    
    dset_umap <- DimPlot(sobj, label = F, label.size = 6, group.by = "orig.ident") + 
            ggtitle("Colored by Dataset")
    
    # Add estimated PC using Charlene's method
    elb_plt <- ElbowPlot(sobj, ndims = 50)
    doub_umap <- DimPlot(sobj, group.by = "DF.classifications") + 
                 ggtitle(paste0("nDoublets    n = ", nrow(sobj@meta.data[sobj@meta.data$DF.classifications == "Doublet",])))
    count.fp <- FeaturePlot(sobj, features =  "nCount_RNA", label = TRUE, label.size = 6)
    feat.fp <- FeaturePlot(sobj, features =  "nFeature_RNA", label = TRUE, label.size = 6)
    mito.fp <- FeaturePlot(sobj, features =  "percent.mt", label = TRUE, label.size = 6)
    doub.fp <- FeaturePlot(sobj, features =  "DF.scores", label = TRUE, label.size = 6)

    hx_merg <- ggplot(sobj@meta.data, aes(nFeature_RNA,  log10(nCount_RNA))) + geom_hex(bins = 100) +
                scale_fill_viridis() + theme_light() + ggtitle(paste0("Merged    n = ", nrow(sobj@meta.data))) +
                lims(y=c(2,6))
    
    pdf(file = paste0(outdir, sobj@meta.data$orig.ident[1], "_qc_report.pdf"), width = 14, height = 14)
    grid.arrange(umap, dset_umap, elb_plt, doub_umap, ncol = 2)
    grid.arrange(doub.fp, count.fp, feat.fp, mito.fp, ncol = 2)
    grid.arrange(hx_merg, ncol =1)
    dev.off()

    ### Add boxplots of qt.doubt scores, counts, features, mito per dataset
}

### Run Merging Pipeline ###
args <- process_args()

obj_list <- read_objects(ipath1 = args$input_path, ipath2 = args$input_path2)
obj_list <- adj_doub_colnames(obj_list)

m.sobj <- merging(obj_list = obj_list, project = args$project_name)

mp.sobj <- sct.harm.processing.m(sobj = m.sobj, dims = 1:25, res = 1, harmony = FALSE, CM = FALSE)

mp.sobj <- qt.norm.drm.scoring(object = mp.sobj, normalize = TRUE)

plotting(sobj = mp.sobj, outdir = args$output_path)

saveRDS(mp.sobj, paste0(args$output_path, args$project_name, "_merged_sobj.RDS"))






