library(Seurat)
library(stringr)
library(scCustomize)
source("~/scToolbox/utils.R")


sobj <- readRDS("/projects/ps-epigen/users/rlan/Liver/manuscript_RNA/Objects/filt_sobj_anno.RDS")

split.by <- "cellclass"

comps <- "/projects/ps-epigen/users/rlan/Liver/manuscript_RNA/DESeq/comparisons.txt" 

outdir <- "/projects/ps-epigen/users/rlan/Liver/RFM_Liver/cellclass_out/matrices/" 



# process comparisons 
if (is.null(comps)){
    # DEFAULT. Run pairwise
    pair_list <- combn(unique(sobj$conditions), 2)
    num_pairs <- ncol(pair_list)
} else {
    comp.table <- read.table(comps, header = T, sep = "\t")
    num_pairs <- nrow(comp.table)
    pair_list <- t(comp.table)
}

Idents(sobj) <- split.by
### Iterate over our groups and comparisons to write csvs for every comparison in each group
if (!is.null(split.by)){
    for (g in unique(sobj[[split.by]][,1])){
        print(g)  	
 	# subset by group	
 	g.sobj <- subset(sobj, idents = g)

	# normalize, scale, call var features, use top 5k
	# Reprocess subset sobj via sct to normalize and obtain variable features
	# scale.data slot is 5k * nCells 
	DefaultAssay(g.sobj) <- "RNA"
	g.sobj <- NormalizeData(g.sobj) %>% 
		  FindVariableFeatures(selection.method = "vst", nfeatures = 5000) %>%
		  ScaleData()

	### create table of percent expression across our groups (split.by)
	percent_expressed <- Percent_Expressing(seurat_object = sobj, features = rownames(g.sobj@assays$RNA@scale.data),
				 threshold = 0)
	
	for (i in 1:num_pairs){
	    cpair <- pair_list[,i]
	    cond1 <- cpair[1]
	    cond2 <- cpair[2]
	    Idents(g.sobj) <- "conditions"
    	    # subset group object by conditions 
	    c.sobj <- subset(g.sobj, idents = c(cond1, cond2))

   	    nmat <- c.sobj@assays$RNA@scale.data
	    nmeta <- c.sobj@meta.data
	    sub.nmeta <- data.frame(row.names = rownames(nmeta), "conditions" = nmeta$conditions)

	    ### subset nmat for features with > 10% representation in target celltype
	    percent_expressed.sub <- percent_expressed[percent_expressed[[g]] > 10, ]
	    nmat.sub <- nmat[rownames(percent_expressed.sub),]
	    
	    cond1 <- str_replace_all(cond1, " ", ".")
            cond2 <- str_replace_all(cond2, " ", ".")

 	    write.csv(nmat.sub, file = paste0(outdir, g, "_", cond1, "_", cond2, "_rfm_mat.csv"))
  	    write.csv(sub.nmeta, file = paste0(outdir, g, "_", cond1, "_", cond2, "_rfm_meta.csv"))

	}
    }
}




