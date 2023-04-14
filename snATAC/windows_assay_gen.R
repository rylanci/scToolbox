library(Seurat)
library(Signac)
library(GenomicRanges)
library(dplyr)
library(ggplot2)
library(reticulate)
library(gridExtra)
library(Matrix)
library(BSgenome.Mmusculus.UCSC.mm10)
library(BSgenome.Hsapiens.UCSC.hg38)
library(harmony)
source("~/scripts/utils.R")
#library(EnsDb.Mmusculus.v79)
#library(EnsDb.Mmusculus.v79)
#library(repr)


# Take file args
args = commandArgs(trailingOnly=TRUE)

# Format args
experiment <- args[1] # $NAME
input_path <- args[2] # outs/ for 10x
destdir <- args[3]
species <- args[4] # human or mouse

print("*** Signac Preprocessing Inputs ***")
print(paste0("Experiment: " ,experiment))
print(paste0("Input path: " ,input_path))
print(paste0("Output dir: ", destdir))
print(paste0("Species: ",species))
print("\n")

sobj_raw <- readRDS(input_path)
print("Seurat Object: ")
print(sobj_raw)
print("Assay info: All assays, Default assay, set Assay")
print(Assays(sobj_raw))
print(DefaultAssay(sobj_raw))
DefaultAssay(sobj_raw) <- "merged_cpeaks"
print(DefaultAssay(sobj_raw))

print("Creating windows/bin matrix")
if (species == "human"){
    bin_matrix <- GenomeBinMatrix(
        cells = colnames(sobj_raw[["ATAC"]]),
        process_n = 200000,
        fragments = Fragments(sobj_raw),
        genome = seqlengths(BSgenome.Hsapiens.UCSC.hg38),
        binsize = 5000
        )
} else if (species =="mouse"){
  bin_matrix <- GenomeBinMatrix(
      cells = colnames(sobj_raw[["ATAC"]]),
      process_n = 200000,
      fragments = Fragments(sobj_raw),
      genome = seqlengths(BSgenome.Mmusculus.UCSC.mm10),
      binsize = 5000
      )
}

print("Create bin assay object")
#bin_assay <- CreateAssayObject(counts = bin_matrix)
bin_assay <- CreateChromatinAssay(bin_matrix, fragments = Fragments(sobj_raw))
#saveRDS(bin_assay, paste0(destdir, experiment, "bin_assay.RDS"))


print("windows #col")
print(ncol(bin_matrix))
print("adding assay to sobj")
sobj_raw[["windows_merg"]] <- bin_assay
DefaultAssay(sobj_raw) <- "windows_merg"

print("processing")
sobj_raw <- RunTFIDF(sobj_raw) %>%
    FindTopFeatures(min.cutoff = 'q0') %>%
    RunSVD() %>%
    RunUMAP(reduction = 'lsi', dims = 2:30) %>%
    FindNeighbors(reduction = 'lsi', dims = 2:30) %>%
    FindClusters(algorithm = 3, resolution = .3, verbose = FALSE) %>%
    peak_calling("windows_merg", "merged_win_cpeaks", species)

DefaultAssay(sobj_raw) <- "merged_win_cpeaks"
sobj_raw <- RunTFIDF(sobj_raw) %>%
    FindTopFeatures(min.cutoff = 'q0') %>%
    RunSVD() %>%
    RunUMAP(reduction = 'lsi', dims = 2:30) %>%
    FindNeighbors(reduction = 'lsi', dims = 2:30) %>%
    FindClusters(algorithm = 3, resolution = .3, verbose = FALSE)



saveRDS(sobj_raw, paste0(destdir, experiment, ".merged_bins.RDS"))

# create plotting function
