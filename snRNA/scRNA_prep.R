suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(DoubletFinder))
suppressPackageStartupMessages(library(plyr))
suppressPackageStartupMessages(library(Nebulosa))
suppressPackageStartupMessages(library(patchwork))
suppressPackageStartupMessages(library(viridis))
suppressPackageStartupMessages(library(future))
suppressPackageStartupMessages(library(harmony))
suppressPackageStartupMessages(library(preprocessCore))
suppressPackageStartupMessages(library(gridExtra))
suppressPackageStartupMessages(library(argparse))
source("~/scripts/utils.R")


######### Functions ###########
process_args <- function(){
# create parser object
parser <- ArgumentParser(description = "Create seurat object from 10x h5 and collect QC metrics")

required_arg_group = parser$add_argument_group('flagged required arguments', 
        'the script will fail if these args are not included')

# These args belong to my group of required flagged arguments
required_arg_group$add_argument("-h5", "--h5_file", required = TRUE,
    help="path to dataset h5 file for object creation")

required_arg_group$add_argument("-l","--library_id", required = TRUE,
     help="library id of the sample (JB_232)") 

required_arg_group$add_argument("-s","--species", required = TRUE,
    help="species of origin for the dataset (mm10 or hg38)")

required_arg_group$add_argument("-o","--output_path", required = TRUE,
    help="output path for plots and objects")

args <- parser$parse_args()
return(args)
}


# Read cellbender h5 as seurat object 
H5_to_sobj <- function(h5_path, species, sample){
    print("Creating sobj from h5")
    H5_data <- Read10X_h5(filename = h5_path, use.names = TRUE)
    sobj_raw <- CreateSeuratObject(counts = H5_data, project = sample)

    ## Add percent mito here as well to plot pre-filter
    if(species == "mm10") {
       sobj_raw[["percent.mt"]] <- PercentageFeatureSet(sobj_raw, pattern = "^mt-")
    } else if(species == "hg38") {
       sobj_raw[["percent.mt"]] <- PercentageFeatureSet(sobj_raw, pattern = "^MT-")
    }

    return(sobj_raw)
}


sobj_filter <- function(sobj_raw, species, outdir, sample){
    print("Filtering")
    ### Filter features by 5,95 quantiles
    quants <- quantile(sobj_raw@meta.data$nFeature_RNA, c(0.05, 0.95))
    #quants
    pct_five <- quants[1]
    nFeatRNA_min <- ceiling(pct_five)
    pct_ninefive <- quants[2]
    nFeatRNA_max <- ceiling(pct_ninefive)
    filt.feat <- subset(sobj_raw, nFeature_RNA > nFeatRNA_min & nFeature_RNA < nFeatRNA_max)   
    
    ### Filter pct.mito by 95 quantile (rounded up)
    quant_mito <- quantile(sobj_raw@meta.data$percent.mt, c(0.95)) #shave top 5%
    percent_mt <- ceiling(quant_mito)
    if (percent_mt == 0) {
       percent_mt = 1
    }
    filt.mito <- subset(filt.feat, percent.mt < percent_mt)

    ### Create filter summary.txt 
    dir.create(paste0(outdir, "QC_data"), showWarnings = FALSE)
    dir.create(paste0(outdir,"QC_data/", sample), showWarnings = FALSE)
    outdir <- paste0(outdir, "QC_data/", sample, "/")
    fileConn <- file(paste0(outdir, "filter_summary.txt"))
    writeLines(c("Filtering Summary",
                 "",
                 paste0("Feature filtering threshold (q5, q95): ", nFeatRNA_min, " , ", nFeatRNA_max ),
                 paste0("Mito filtering threshold (q95): ", percent_mt),
                 "", 
                 paste0("# Nuclei Before Filtering: ", nrow(sobj_raw@meta.data)),
                 paste0("# Nuclei After Feature Filter: ", nrow(filt.feat@meta.data)),
                 paste0("# Nuclei After Mito Filter: ", nrow(filt.mito@meta.data)),
                 paste0("# Pct Remaining: ", nrow(filt.mito@meta.data) / nrow(sobj_raw@meta.data) * 100)     
    ), fileConn)
    close(fileConn)
    
    ### Create csv summary of each metric for prefilt and postfilt objects 
    metrics <- c("nCount_RNA", "nFeature_RNA", "percent.mt")
    r.qc.table <- data.frame()
    f.qc.table <- data.frame()
    qc.table <- data.frame()
    for (m in metrics){
        # r is for raw, f is for filtered
        # [[1]] is to convert list to double for summary() formating 
        r.s <- summary(sobj_raw[[m]][[1]])
        r.s <- c("raw", m, r.s)
        f.s <- summary(filt.mito[[m]][[1]])
        f.s <- c("filt", m, f.s)
        # Row is summary(min, q1, median, mean, q3, max)
        r.qc.table <- rbind(r.qc.table, r.s)
        f.qc.table <- rbind(f.qc.table, f.s)
    }
    colnames(r.qc.table) <- c("Status", "Metric", "MinScore", "Q1", "Median", "Mean", "Q3", "MaxScore")
    colnames(f.qc.table) <- c("Status", "Metric","MinScore", "Q1", "Median", "Mean", "Q3", "MaxScore")
    m.qc.table <- rbind(r.qc.table, f.qc.table)
    write.table(x = m.qc.table, file = paste0(outdir, "qc_table.txt"), quote = FALSE, sep = "    ", row.names = F)
    
    ### Return filtered object                     
    return(filt.mito)
}

find_nexp <- function(inp_cells, outdir, sample) {
    print("DoubletFinding: Calculating nexp")
    # Inp_cells is the number of cells in the object after seurat standard QC filtering for mito and nFeatures
    # Rec_cells is the number of input cells rounded down by 1000,
    # Rec_cells will be used to find the corresponding chosen percent in the multiplet_rate table
    
    # If rec cells is above 10k or below 1k choose these params 
    # else use the 10x doublet rate table
    rec_cells <- round_any(inp_cells, 1000, f=floor)
    if(inp_cells > 10000) {
        chos_pt <- 0.076
        rec_cells <- 10000
    } else if(inp_cells < 1000){
        chos_pt <- 0.004
        rec_cells <- 500
    } else {
        # Calculate percentage based on 10X multiplet rate table
        mtable <- read.delim("/home/rlancione/misc/10x_doublet_rate.txt", sep=" ", header = T, stringsAsFactors = F)
        chos_pt <- mtable[which(mtable$recov_cells == rec_cells),c("percent")]
    }
    
    calc_pt <- (inp_cells*chos_pt)/rec_cells
    calc_nexp <- round(calc_pt*inp_cells)
    
    
    # Create doub summary 
    dir.create(paste0(outdir, "QC_data"), showWarnings = FALSE)
    dir.create(paste0(outdir,"QC_data/", sample), showWarnings = FALSE)
    outdir <- paste0(outdir, "QC_data/", sample, "/")
    fileConn <- file(paste0(outdir, "exp_doublet_summary.txt"))
    writeLines(c("Expected Doublet Summary",
                 "",
                 paste0("Number of cells after filtering: ", inp_cells),
                 paste0("Chosen percent: ", chos_pt),
                 paste0("Calculated percentage of doublets: ", calc_pt),
                 paste0("Calculated expected doublets: ", calc_nexp)
                   
    ), fileConn)
    close(fileConn)
    
    # Return calc expected # of doublets 
    return(calc_nexp)
}

Doublet_Finder <- function(sobj, pN_inp, outdir, sample){
    print("DoubletFinding: Calculating pK...")
    # Long std output on paramSweep. Tried to supress. Experiment with cores. 
    sweep.sobj <- paramSweep_v3(sobj, PCs = 1:50, sct = TRUE)
    sweep.stats.sobj <- summarizeSweep(sweep.sobj, GT = FALSE)
    bcmvn.sobj <- find.pK(sweep.stats.sobj)
    #datatable(bcmvn.sobj, rownames = F)
    
    # choosing maxima of BCmetric
    pK <- as.numeric(as.character(bcmvn.sobj$pK))
    top_pK <- pK[which(bcmvn.sobj$BCmetric %in% max(bcmvn.sobj$BCmetric))]
    
    # Calculating nEXP. Apply find_nexp funtion
    cells_aft_filt <- nrow(sobj@meta.data)
    calc_nExp <- find_nexp(inp_cells = cells_aft_filt, outdir = outdir, sample = sample)
    
    ####### Run Doublet Finder ############
    sobj_seurat_DF <- doubletFinder_v3(sobj, PCs = 1:50, pN = pN_inp , 
                                   pK = top_pK, nExp = calc_nExp, reuse.pANN = FALSE, sct = TRUE)
    
    # Change DF colnames with regex and this crazy syntax  
    colnames(sobj_seurat_DF@meta.data)[grep(pattern = "pANN", x = colnames(sobj_seurat_DF@meta.data))] <- "DF.scores"
    colnames(sobj_seurat_DF@meta.data)[grep(pattern = "DF.class", x = colnames(sobj_seurat_DF@meta.data))] <- "DF.classifications"

    
    return(sobj_seurat_DF)
}

plotting <- function(raw.sobj, filt.sobj, outdir, sample){
    print("Plotting")
    dir.create(paste0(outdir, "QC_plots/"), showWarnings = FALSE)
    dir.create(paste0(outdir,"QC_plots/", sample), showWarnings = FALSE)
    outdir <- paste0(outdir, "QC_plots/", sample, "/")
    
    ### Before & after filtering 
    # Viridis hexplot features by counts 
    hx_raw <- ggplot(sobj.raw@meta.data, aes(nFeature_RNA,  log10(nCount_RNA))) + geom_hex(bins = 100) +
                scale_fill_viridis() + theme_light() + ggtitle(paste0("Prefilter    n = ", nrow(sobj.raw@meta.data))) +
                lims(y=c(2,6))
    hx_filt <- ggplot(filt.sobj@meta.data, aes(nFeature_RNA,  log10(nCount_RNA))) + geom_hex(bins = 100) +
                scale_fill_viridis() + theme_light() + ggtitle(paste0("Post Filter   n = ", nrow(filt.sobj@meta.data))) + 
                lims(y=c(2,6))

    # Pre & Post filter boxplots
    r.count.bx <- ggplot(sobj.raw@meta.data, aes(y=log10(nCount_RNA))) + 
                    geom_boxplot(fill='lightblue', color="black") +
                    ggtitle("PreFilter nCounts") + theme_light() + 
                    lims(y=c(2,6))
    r.feat.bx <- ggplot(sobj.raw@meta.data, aes(y=log10(nFeature_RNA))) + 
                    geom_boxplot(fill='lightblue', color="black") +
                    ggtitle("PreFilter nFeatures") + theme_light()+ 
                    lims(y=c(2,5))
    r.mito.bx <- ggplot(sobj.raw@meta.data, aes(y=log10(percent.mt))) + 
                    geom_boxplot(fill='lightblue', color="black") +
                    ggtitle("PreFilter pct.mito") + theme_light() #+ 
                    #lims(y=c(2,6))

    f.count.bx <- ggplot(filt.sobj@meta.data, aes(y=log10(nCount_RNA))) + 
                    geom_boxplot(fill='lightblue', color="black") +
                    ggtitle("PostFilter nCounts") + theme_light() + 
                    lims(y=c(2,6))
    f.feat.bx <- ggplot(filt.sobj@meta.data, aes(y=log10(nFeature_RNA))) + 
                    geom_boxplot(fill='lightblue', color="black") +
                    ggtitle("PostFilter nFeatures") + theme_light() + 
                    lims(y=c(2,5))
    f.mito.bx <- ggplot(filt.sobj@meta.data, aes(y=log10(percent.mt))) + 
                    geom_boxplot(fill='lightblue', color="black") +
                    ggtitle("PostFilter pct.mito") + theme_light()# + 
                   # lims(y=c(2,6))
 

    ### After Filtering 
    umap <- DimPlot(filt.sobj, label = T, label.size = 6) + ggtitle("PostFilter SCT Embedding")
    doub_umap <- DimPlot(filt.sobj, group.by = "DF.classifications") + 
                 ggtitle(paste0("nDoublets    n = ", nrow(filt.sobj@meta.data[filt.sobj@meta.data$DF.classifications == "Doublet",])))
    count.fp <- FeaturePlot(filt.sobj, features =  "nCount_RNA", label = TRUE, label.size = 6)
    feat.fp <- FeaturePlot(filt.sobj, features =  "nFeature_RNA", label = TRUE, label.size = 6)
    mito.fp <- FeaturePlot(filt.sobj, features =  "percent.mt", label = TRUE, label.size = 6)
    doub.fp <- FeaturePlot(filt.sobj, features =  "DF.scores", label = TRUE, label.size = 6)

    
    pdf(file = paste0(outdir, filt.sobj@meta.data$orig.ident[1], "_qc_report.pdf"), width = 10, height = 14)
    grid.arrange(hx_raw, hx_filt,r.count.bx, f.count.bx, r.feat.bx, 
                 f.feat.bx, r.mito.bx, f.mito.bx, ncol = 2)
    grid.arrange(umap, doub_umap, doub.fp, count.fp, feat.fp, mito.fp, 
                 ncol = 2)
    dev.off()
    
    ### Save plot independently 
    dir.create(paste0(outdir, "indi_plots"))
    ggsave(filename = "raw_hx_plot.pdf", plot = hx_raw, path = paste0(outdir, "/indi_plots"), device = "pdf")
    ggsave(filename = "filt_hx_plot.pdf", plot = hx_filt, path = paste0(outdir, "/indi_plots"), device = "pdf")
    ggsave(filename = "umap.pdf", plot = umap, path = paste0(outdir, "/indi_plots"), device = "pdf")
    ggsave(filename = "doub_umap.pdf", plot = doub_umap, path = paste0(outdir, "/indi_plots"), device = "pdf")
    ggsave(filename = "ncount_fp.pdf", plot = count.fp, path = paste0(outdir, "/indi_plots"), device = "pdf")
    ggsave(filename = "nfeat_fp.pdf", plot = feat.fp, path = paste0(outdir, "/indi_plots"), device = "pdf")
    ggsave(filename = "pct.mito_fp.pdf", plot = mito.fp, path = paste0(outdir, "/indi_plots"), device = "pdf")
    ggsave(filename = "doub_fp.pdf", plot = doub.fp, path = paste0(outdir, "/indi_plots"), device = "pdf")
}


save_rds <- function(sobj, outdir, sample){
    dir.create(paste0(outdir,"Objects"), showWarnings = FALSE)
    saveRDS(sobj, paste0(outdir,"Objects/", sample,".RDS")) 
}

### main ###
args <- process_args()

#print(args$H5)
#print(args$SP)

sobj.raw <- H5_to_sobj(args$h5_file, species = args$species, sample = args$library_id)
sobj.raw

sobj.filt <- sobj_filter(sobj_raw = sobj.raw, species = args$species, outdir = args$output_path, sample = args$library_id)
sobj.filt

psobj.filt <- sct.harm.processing(sobj.filt, dims=1:25, res=.05, harmony = FALSE, CM = FALSE)
psobj.filt

psobj.filt.db <- Doublet_Finder(sobj = psobj.filt, pN_inp = .25, outdir = args$output_path, sample = args$library_id)
psobj.filt.db

plotting(raw.sobj = sobj.raw, filt.sobj = psobj.filt.db, outdir = args$output_path, sample = args$library_id)
save_rds(sobj = psobj.filt.db, outdir = args$output_path, sample = args$library_id)

print("done")

# notes:
# Save raw version of object as well as processed
