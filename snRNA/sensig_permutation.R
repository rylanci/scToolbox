library(Seurat)
library(stats)
library(gprofiler2)


### version 2. Remove abs value 
sen_sig_score_rna_V2 <- function(seu_obj, gene_list) {

    DefaultAssay(seu_obj) <- 'RNA'
    # if RNA assay is unprocessed, run for sensig
    if (sum(seu_obj@assays$RNA@counts) == sum(seu_obj@assays$RNA@data)) {
        seu_obj <- NormalizeData(seu_obj)
    }
    
    # the scale data assay will change with every gene list
    seu_obj <- suppressMessages(ScaleData(seu_obj, features = gene_list))

   # print("Calculating SenSig score")
    cts_sct <- seu_obj@assays$RNA@scale.data
    deg_cts <- subset(cts_sct, rownames(cts_sct) %in% gene_list)

    scores <- colSums(deg_cts)
    scores <- trunc(scores*10^2)/10^2
   # seu_obj[[paste0(gl_col_name, "_SenSig")]] <- scores
   # return(seu_obj)
    return(scores)
    
}

sensig_permutation <- function(sobj, n_sample, n_permutations, sig.thresh, score.col, name){
	
    # collect permuted scores for each cell
    counts <- sobj@assays$RNA@counts
    perm.df <- data.frame(row.names = row.names(sobj@meta.data))
    for (i in seq(1, n_permutations)){
        g.sample <- sample(x = row.names(counts), size = n_sample, replace = TRUE)
        p.scores <- sen_sig_score_rna_V2(sobj, g.sample)
        perm.df[paste0("p_", i)] <- p.scores
    
    }

    ## Calculate pvalues
    pvalues <- c()
    scores <- sobj[[score.col]][,1]
    for (r in seq(1,nrow(perm.df))){
	# perm row pertains to a cell in perm.df
	perm.row <- perm.df[r,]

	# score pertains to the sensig score of that cell
	score <- scores[r]
	obs <- 0
	for (p in perm.row){if (p > score){obs <- obs + 1}}

	pvalues <- append(pvalues, obs/n_permutations)
    }

    psig <- c()
    for (i in pvalues){
	if (i < sig.thresh){psig <- append(psig, "sensecent")
	} else {psig <- append(psig, "normal")}
    }
   
    
    ### Multiple test correction
    qvalues <- p.adjust(pvalues, method = "BH", n = length(pvalues))
	
    qsig <- c()
    for (i in qvalues){
	if (i < sig.thresh){qsig <- append(qsig, "sensecent")
	} else {qsig <- append(qsig, "normal")}
    }


    sobj[[paste0(name, "_pval")]] <- pvalues
    sobj[[paste0(name, "_padj")]] <- qvalues
    sobj[[paste0(name, "_", sig.thresh, "_psig")]] <- psig
    sobj[[paste0(name, "_", sig.thresh, "_qsig")]] <- qsig

    return(sobj)

}






