# rlibs
library(Signac)
library(Seurat)
library(tidyverse)
library(gridExtra)
library(Matrix)
library(EnsDb.Mmusculus.v79)
library(EnsDb.Hsapiens.v86)
library(viridis)

# Take file args
args = commandArgs(trailingOnly=TRUE)

# Format args
experiment <- args[1] # $NAME
input_path <- args[2]
destdir <- args[3]
min_tss <- as.double(args[4])
min_reads <- as.double(args[5])
species <- args[6]

# Assign path vars from args
meta_path <- paste0(input_path,experiment,".qc_metrics.txt")
frag_path <- paste0(input_path,experiment,"_fragments_sorted_barcode_modified.bed.gz")
mat_path <- paste0(input_path,experiment,".long_fmt_mtx.txt.gz")

print("***Signac Preprocessing Inputs***")
print(paste0("Experiment: " ,experiment))
print(paste0("Output dir: ", destdir))
print(paste0("Metadata path: ",meta_path))
print(paste0("Frag file path: ", frag_path))
print(paste0("Matrix file path: ",mat_path))
print(paste0("Min reads per cell: ",min_reads))
print(paste0("Min tss per cell: ",min_tss))
print(paste0("Species: ",species))

# Create chromatin_assay and Seurat object
data <- read.delim(mat_path, header = F)
data$V1 = as.factor(data$V1)
data$V2 = as.factor(data$V2)
data$V3 = as.numeric(data$V3)
sparse_data <- with(data, sparseMatrix(i=as.numeric(V2), j=as.numeric(V1), x=V3, dimnames=list(levels(V2), levels(V1))))
# Conditional chromatin assays based on species
if (species == "mouse"){
  chrom_assay <- CreateChromatinAssay(counts = sparse_data, sep = c(":", "-"), genome = "mm10", fragments = frag_path)
} else if (species == "human"){
  chrom_assay <- CreateChromatinAssay(counts = sparse_data, sep = c(":", "-"), genome = "hg38", fragments = frag_path)
}

metadata <- read.table(file = meta_path, header = TRUE)
sobj_raw <- CreateSeuratObject(counts = chrom_assay, assay = "windows", meta.data = metadata)

# Add gene annotations to raw objects
if (species == "mouse"){
  # extract gene annotations from EnsDb
  annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Mmusculus.v79)
  # change to UCSC style since the data was mapped to hg19
  seqlevelsStyle(annotations) <- 'UCSC'
  genome(annotations) <- "mm10"
  # add the gene information to the object
  Annotation(sobj_raw) <- annotations
} else if (species == "human"){
  annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
  seqlevelsStyle(annotations) <- 'UCSC'
  genome(annotations) <- "hg38"
  Annotation(sobj_raw) <- annotations
}

# QC
# compute nucleosome signal score per cell
sobj_raw <- NucleosomeSignal(object = sobj_raw)
sobj_raw$nucleosome_group <- ifelse(sobj_raw$nucleosome_signal > 4, 'NS > 4', 'NS < 4')

# compute TSS enrichment score per cell
sobj_raw <- TSSEnrichment(object = sobj_raw, fast = FALSE)
print("TSS Done")
sobj_raw$high.tss <- ifelse(sobj_raw$TSS.enrichment > 2, 'High', 'Low')

# add experiment and mitochondrial metadata
sobj_raw[["experiment"]]<- experiment
# different case convention for mouse vs human RNA. So use condition
if (species == "mouse"){
  sobj_raw[["percent.mt"]] <- PercentageFeatureSet(sobj_raw, pattern = "^Mt-") # IF human then "Mt"
} else if(species == "human"){
  sobj_raw[["percent.mt"]] <- PercentageFeatureSet(sobj_raw, pattern = "^MT-")
}

# Filtering. Save filtered object
sobj_filt <- subset(x = sobj_raw, subset = unique_usable_reads > min_reads & TSS.enrichment > min_tss)
dir.create(path= paste0(destdir,"signac_objects/"), showWarnings= FALSE)
dir.create(path= paste0(destdir,"signac_objects/filtered"), showWarnings= FALSE)
saveRDS(sobj_filt, paste0(destdir, "signac_objects/filtered/", experiment, ".filt.RDS"))


# QC Plots
qc1 <- TSSPlot(sobj_raw) + NoLegend()
qc2 <- ggplot(sobj_raw@meta.data, aes(TSS.enrichment,	log10(unique_usable_reads))) + geom_hex(bins = 100) +
                scale_fill_viridis() + theme_light() + geom_text(x=29, y=6, label=paste0("n = ", nrow(sobj_raw@meta.data))) + lims(x= c(0,30),y=c(0,6))
#qc3 <- FragmentHistogram(object = sobj_raw, group.by = 'nucleosome_group', region = 'chr1-1-10000000')
qc4 <- TSSPlot(sobj_filt) + NoLegend()
qc5 <- ggplot(sobj_filt@meta.data, aes(TSS.enrichment, log10(unique_usable_reads))) + geom_hex(bins = 100) +
            scale_fill_viridis() + theme_light() + geom_text(x=29, y=6, label=paste0("n = ", nrow(sobj_filt@meta.data))) + lims(x= c(0,30),y=c(0,6))

dir.create(path= paste0(destdir,"QC_plots/"), showWarnings= FALSE)
pdf(file=paste0(destdir,"QC_plots/",experiment,"_qcplot.pdf"), title= experiment, width= 12, height= 5, compress=FALSE)
grid.arrange(qc1, qc2, ncol =2, padding= unit(0,"line"))
grid.arrange(qc4, qc5, ncol =2, padding= unit(0,"line"))
dev.off()

# Funct to add
#  write qc metric .csv for each dataset containing: mean tss, mean unqiue reads, mean total reads, maybe others
