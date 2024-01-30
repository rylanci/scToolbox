library(Seurat)
library(Signac)
library(UCell)
library(dplyr)
library(ggplot2)
library(patchwork)
library(gridExtra)
library(tidyr)
library(openxlsx)
library(stringr)
library(circlize)
library(ComplexHeatmap)
library(argparse)

process_args <- function(){
    # create parser object
    parser <- ArgumentParser(description = "Create seurat object from 10x h5 and collect QC metrics")

    # Create parser arg group for our required args
    required_arg_group = parser$add_argument_group('flagged required arguments',
        'the script will fail if these args are not included')

    ### start args
    required_arg_group$add_argument("-s","--sobj_path", required = TRUE,
	help="path to output location")

    required_arg_group$add_argument("-e","--excel_path", required = TRUE,
	help="path to output location")

    required_arg_group$add_argument("-o","--output_path", required = TRUE,
        help="path to output location")

    required_arg_group$add_argument("-n","--name", required = TRUE,
        help="name of analysis for output files")


    # optional args
    parser$add_argument("-sn","--sheet_number",
        help="number of cores to run FindMarkers with",
        default = 1)

    args <- parser$parse_args()
    return(args)
}

args <- process_args()

sobj <- readRDS(args$sobj_path)
DefaultAssay(sobj) <- "RNA"
cc <- read.xlsx(xlsxFile = args$excel_path, sheet = args$sheet_number)
print(head(cc))

colnames(cc) <- c("CellType", "Marker")
celltypes <- unique(cc$CellType)

# Create lists for each celltype and markers 
gene.sets <- lapply(X = celltypes, FUN = function(X){
    cc.markers <- cc[grep(pattern = paste0("^", X ,"$"), x = cc$CellType),]$Marker
    
    return(cc.markers)
})
gene.sets <- setNames(gene.sets, celltypes)
print(head(gene.sets))

#### Run UCell 
sobj <- AddModuleScore_UCell(sobj, features = gene.sets)
saveRDS(sobj, paste0(args$output_path, args$name, "_Ucell_sobj.RDS"))



# Create matrices of mean and median scores 
Idents(sobj) <- "seurat_clusters"
i.df.avg <- data.frame(row.names = celltypes)
i.df.med <- data.frame(row.names = celltypes)
celltypes <- str_replace_all(string = celltypes, pattern = "/", replacement = ".")
celltypes <- str_replace_all(string = celltypes, pattern = " ", replacement = ".")

for(i in seq(0, length(unique(sobj@meta.data$seurat_clusters)) -1 )){
    t.sobj <- subset(sobj, idents = i)
    
    avg.score.list <- list()
    med.score.list <- list()
    for(c in celltypes){
        avg <- sum(t.sobj[[paste0(c, "_UCell")]]) / length(WhichCells(t.sobj))
        med <- median(as.matrix(t.sobj[[paste0(c, "_UCell")]]))
        avg.score.list <- append(avg.score.list, avg)
        med.score.list <- append(med.score.list, med)
    }
    
    i.df.avg[[paste0("c_", i)]] <- avg.score.list
    i.df.med[[paste0("c_", i)]] <- med.score.list
    
}


# Format matrices 
i.mat.avg <- as.matrix(i.df.avg)
i.mat.med <- as.matrix(i.df.med)

i.mat.avg <- matrix(as.numeric(i.mat.avg),    # Convert to numeric matrix
                  ncol = ncol(i.mat.avg))
i.mat.med <- matrix(as.numeric(i.mat.med),    # Convert to numeric matrix
                  ncol = ncol(i.mat.med))

rownames(i.mat.avg) = rownames(i.df.avg)
colnames(i.mat.avg) = colnames(i.df.avg)
rownames(i.mat.med) = rownames(i.df.med)
colnames(i.mat.med) = colnames(i.df.med)


# Create Heatmaps 
options(repr.plot.width = 18)
options(repr.plot.height = 12)
imm.avg.hm <- Heatmap(i.mat.avg, cluster_rows = FALSE, cluster_columns = FALSE, 
                      heatmap_legend_param = list(title = "avg"), rect_gp = gpar(col = "black", lwd = 2))
imm.med.hm <- Heatmap(i.mat.med, cluster_rows = FALSE, cluster_columns = FALSE, 
                      heatmap_legend_param = list(title = "med"), rect_gp = gpar(col = "black", lwd = 2))


pdf(file = paste0(args$output_path, args$name, "_Ucell_heatmaps.pdf"), width = 18, height = 12, title = "Annotation Heatmaps")
imm.avg.hm + imm.med.hm
dev.off()








