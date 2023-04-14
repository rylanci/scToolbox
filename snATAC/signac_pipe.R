#This script will process, remove doublets, cluster and produce plots for scATAC datasets.
### rlibs
library(Seurat)
library(Signac)
library(tidyverse)
library(reticulate)
library(gridExtra)
library(Matrix)
library(EnsDb.Mmusculus.v79)
library(EnsDb.Hsapiens.v86)

### pylibs
sci <- import("scipy")
scrub <- import("scrublet")

# Kai's scrublet script (with some convinient edits)
source_python("/home/rlancione/scripts/ATAC/doub_detect.py")

################### Functions ###################

processing <- function(sobj, assay){
  # Define assay for processing
  DefaultAssay(sobj) <- assay
  # process object
  sobj <- RunTFIDF(sobj) %>%
              FindTopFeatures(min.cutoff = 'q0') %>%
              RunSVD() %>%
              RunUMAP(reduction = 'lsi', dims = 2:30) %>%
              FindNeighbors(reduction = 'lsi', dims = 2:30) %>%
              FindClusters(algorithm = 3, resolution = 1.2, verbose = FALSE)

  return(sobj)
}


scrublet <- function(sobj, assay, meta_suffix){
  # Convert seurat assay into scipy sparse matrix
  sobj_matrix <- eval(parse(text = paste0("sobj@assays$", assay, "@counts")))
  sci_counts <- sci$sparse$csr_matrix(t(sobj_matrix))
  # Initialize scrublet object
  scrub_mat <- scrub$Scrublet(sci_counts, expected_doublet_rate=0.06)

  # Calculate doublet scores/ predicted doublets
  scores_pred <- scrub_mat$scrub_doublets(min_counts=2,
                          min_cells=3,
                          min_gene_variability_pctl=85,
                          n_prin_comps=as.integer(30))

  # Add scores/ prediction to seurat object metadata
  scores_col <- paste0("scrub_scores_",meta_suffix)
  pred_col <- paste0("scrub_predictions_",meta_suffix)
  sobj[[scores_col]] <- scores_pred[1]
  sobj[[pred_col]] <- scores_pred[2]

  return(sobj)
}


plotting <- function(sobj, assay, outdir, scrub_meta_suffix){
  # Create outdirs
  exp <- sobj@meta.data$experiment[1]
  dir.create(path= paste0(outdir), showWarnings= FALSE)
  dir.create(path= paste0(outdir,"/",exp), showWarnings= FALSE)

  # Define assay to call peaks
  DefaultAssay(sobj) <- assay
  # Create plots
  main_dp <- DimPlot(sobj, label = T) + ggtitle(paste0(exp,"_UMAP ", "ncells=",nrow(sobj@meta.data))) #+
            #geom_text(x=29, y=6, label=paste0("n = ", nrow(sobj@meta.data))) + lims(x= c(0,30),y=c(0,6))
  #dc <- DepthCor(sobj)
  #eb <- ElbowPlot(sobj, ndims = 30, reduction = "lsi")
  uur_fp <- FeaturePlot(sobj, features = "unique_usable_reads", order = T)
  tss_fp <- FeaturePlot(sobj, features = "TSS.enrichment", order = T) # change to log10
  nfw_fp <- FeaturePlot(sobj, features = "nFeature_windows", order = T)
  # scrublet plots
  scrub_dp <- DimPlot(sobj, group.by = paste0("scrub_pred_", scrub_meta_suffix), order = c('TRUE', 'FALSE'), cols = c("grey", "blue"))
  scrub_fp <- FeaturePlot(sobj, features = paste0("scrub_score_", scrub_meta_suffix), order = T)
  bgmm_dp <- DimPlot(sobj, group.by = paste0("BGMM_pred_", scrub_meta_suffix), order = c('TRUE', 'FALSE'), cols = c("grey", "blue"))
  bgmm_fp <- FeaturePlot(sobj, features = paste0("BGMM_prob_", scrub_meta_suffix), order = T)

  pdf(file=paste0(outdir,"/",exp,"/",exp,"_",assay, "_processing_report.pdf"), title= paste0(exp,"_",assay), width= 12, height= 10, compress=FALSE)
  #print(dp)
  grid.arrange(main_dp, uur_fp, tss_fp, nfw_fp, ncol =2, padding= unit(0,"line"))
  grid.arrange(scrub_dp, scrub_fp, bgmm_dp, bgmm_fp, ncol =2, padding= unit(0,"line"))
  dev.off()
}


peak_calling <- function(sobj, input_assay, output_assay, species){
  # Define assay to call peaks from
  DefaultAssay(sobj) <- input_assay
  # Call peaks. Use peaks + frags to create counts
  peaks <- CallPeaks(object = sobj, group.by = "seurat_clusters")
  frags <- Fragments(object = sobj)
  peak_counts <- FeatureMatrix(fragments = c(frags), cells = rownames(sobj@meta.data), features = peaks)

  # Define new assay for new peak calls
  if (species == "mouse"){
      peak_assay <- CreateChromatinAssay(
          counts = peak_counts,
          sep = c(":", "-"),
          genome = "mm10",
          fragments = frags,
          min.cells = 1
          )
        }
  else if(species == "human"){
      peak_assay <- CreateChromatinAssay(
         counts = peak_counts,
         sep = c(":", "-"),
         genome = "hg38",
         fragments = frags,
         min.cells = 1
        )
      }

  # Add new peak assay
  sobj[[paste0(output_assay)]] <- peak_assay
  return(sobj)
}

write_obj <- function(sobj, outdir){
    # Create outdir
    dir.create(path= paste0(outdir), showWarnings= FALSE)
    # write object
    saveRDS(sobj, paste0(outdir ,"/", sobj@meta.data$experiment[1], "_filt_psobj.RDS"))
}

signac_pipe <- function(obj_path, out_path, species){
    # Initial dim reduction & peak calling
    sobj <- readRDS(obj_path) %>%
    	processing(assay = "windows") %>%
    	peak_calling(input_assay = "windows", output_assay="cpeaks",species=species) %>%
      processing(assay = "cpeaks") %>%
    	peak_calling(input_assay = "cpeaks", output_assay="cpeaks2",species=species)

    # Kai's Scrublet + Baysian Gaussian Mixture Model Doublet Calling
    doub_data_win <- detectDoublet(input_matrix=sobj@assays$windows@counts, BGMM_pval=.05, assay="windows")
    doub_data_cpeaks <- detectDoublet(input_matrix=sobj@assays$cpeaks@counts, BGMM_pval=.05, assay="cpeaks")
    # add doub metadata
    sobj@meta.data$scrub_score_windows <- doub_data_win[[1]][[1]]
    sobj@meta.data$scrub_pred_windows <- doub_data_win[[2]][[1]]
    sobj@meta.data$BGMM_prob_windows <- doub_data_win[[3]][[1]]
    sobj@meta.data$BGMM_pred_windows <- doub_data_win[[4]][[1]]
    sobj@meta.data$scrub_score_cpeaks <- doub_data_cpeaks[[1]][[1]]
    sobj@meta.data$scrub_pred_cpeaks <- doub_data_cpeaks[[2]][[1]]
    sobj@meta.data$BGMM_prob_cpeaks <- doub_data_cpeaks[[3]][[1]]
    sobj@meta.data$BGMM_pred_cpeaks <- doub_data_cpeaks[[4]][[1]]

    # Plotting Obvi.
    plotting(sobj, assay = "windows", outdir = paste0(out_path,"processing_reports"), scrub_meta_suffix = "windows")
    plotting(sobj, assay = "cpeaks", outdir = paste0(out_path,"processing_reports"), scrub_meta_suffix = "cpeaks")

    write_obj(sobj, outdir=paste0(out_path,"signac_objects/processed"))
}


######################## run pipe ############################

# args
args = commandArgs(trailingOnly=TRUE)
ipath <- args[1]
opath <- args[2]
species <- args[3]
# Log inputs
print("***clustering inputs***")
print(paste0("Input obj path:", ipath))
print(paste0("Output path:", opath))
print(paste0("Species:", species))

# Run signac_pipe
signac_pipe(ipath, opath, species)

# Features to add:
## Future multithreading
## Doublet plot info. # called doublets.
## Maybe histogram of scores  
