suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(DoubletFinder))
#library(Nebulosa)
#library(patchwork)
#library(viridis)
#library(future)
#library(harmony)
#library(preprocessCore)
#library(gridExtra)
library(argparse)
source("~/scripts/utils.R")


######### Functions ###########
process_args <- function(){
# create parser object
parser <- ArgumentParser(description = "Create seuat object from 10x h5 and collect QC metrics")

parser$add_argument("-h5", "--h5_path",
    help="path to dataset h5 file for object creation")

parser$add_argument("-s", "--species",
    help="species of origin for the dataset (mm10 or hg38)")

parser$add_argument("-l", "--lib_id",
    help="library id of the sample")

parser$add_argument("-o", "--out_path",
     help="desired path for output to be saved")

#args <- parser$parse_args(c('--input_path','--output_path'))
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


sobj_filter <- function(sobj_raw, species, outdir){
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
    outdir <- paste0(outdir, "QC_data/")
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

find_nexp <- function(inp_cells, outdir) {
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
    outdir <- paste0(outdir, "QC_data")
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

Doublet_Finder <- function(sobj, pN_inp, outdir){
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
    calc_nExp <- find_nexp(inp_cells = cells_aft_filt, outdir = outdir)
    
    ####### Run Doublet Finder ############
    sobj_seurat_DF <- doubletFinder_v3(sobj, PCs = 1:50, pN = pN_inp , 
                                   pK = top_pK, nExp = calc_nExp, reuse.pANN = FALSE, sct = TRUE)
    
    # Change DF colnames with regex and this crazy syntax  
    colnames(sobj@meta.data)[grep(pattern = "pANN", x = colnames(sobj@meta.data))] <- "DF.scores"
    colnames(sobj@meta.data)[grep(pattern = "DF.class", x = colnames(sobj@meta.data))] <- "DF.classifications"
    
    return(sobj_seurat_DF)
}


### main ###
args <- process_args()

print(args$h5_path)
print(args$species)

sobj.raw <- H5_to_sobj(args$h5_path, species = args$species, sample = args$lib_id)
sobj.raw

sobj.filt <- sobj_filter(sobj_raw = sobj.raw, species = args$species, outdir = args$out_path)
sobj.filt

psobj.filt <- sct.harm.processing(sobj.filt, dims=1:25, res=.05, harmony = FALSE, CM = FALSE)
psobj.filt

psobj.filt.db <- Doublet_Finder(sobj = psobj.filt, pN_inp = .25, outdir = args$out_path)
psobj.filt.db


#print(args$input_path)
#print(args$output_path)
