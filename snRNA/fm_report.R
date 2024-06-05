suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(future))
suppressPackageStartupMessages(library(tableHTML))
suppressPackageStartupMessages(library(openxlsx))
suppressPackageStartupMessages(library(clusterProfiler))
suppressPackageStartupMessages(library(org.Hs.eg.db))
suppressPackageStartupMessages(library(org.Mm.eg.db))
### This script will contain functions to run findmarkers and write a report on the output

# Functions
run_fm <- function(sobj, future = FALSE, n.workers = 1, prep.sct = FALSE, 
    test = "wilcox", idents = "seurat_clusters", logfc.threshold = 0.25){

    if (future == TRUE){
	plan(multicore, workers = n.workers)
	options(future.globals.maxSize = 2000 * 1024^2)
    }

    if (prep.sct == TRUE){
	sobj <- PrepSCTFindMarkers(sobj, assay = "SCT", verbose = FALSE)
    }

    Idents(sobj) <- idents
    res <- FindAllMarkers(sobj, test.use = test, logfc.threshold = logfc.threshold)

    return(res)
}


# Writes spreadsheet of top DE genes
write.top.n.xlsx <- function(markers, file, group.by = "cluster", n = 100){
    wb <- createWorkbook("TopMarkers")

    for (c in unique(markers[[group.by]])){
	if (c == 0){c <- "0"}
        addWorksheet(wb, c)
        tdf <- head(markers[markers[[group.by]] == c, ], n = n)
        writeData(wb, sheet = c, x = tdf)
    }

    saveWorkbook(wb = wb, file = file, overwrite = TRUE)
}

### Changes to make 
#fm_kmeans <- function(){
	# think about this later
#}


run_presto <- function(sobj, group.by = "seurat_clusters", n = 10, outdir){
    
    pres <- presto:::wilcoxauc.Seurat(X = sobj,
                                   group_by = group.by,
                                   assay = 'data',
                                   seurat_assay = 'RNA')
    
    # sort results by auc and take top 10 
    pres.sort <- pres %>% arrange(group, desc(auc))
    top10 <- pres.sort  %>% group_by(group) %>% slice_head(n = n)

    dp <- DotPlot(object = sobj, assay = "RNA", features = unique(top10$feature), 
		cols = c("lightblue", "darkblue"), cluster.idents = TRUE) + coord_flip()

    suppressMessages(ggsave(filename = paste0("presto_dotplot.pdf"), plot = dp, device = "pdf", 
		path = paste0(outdir), width = 14, height = 30))

    write.top.n.xlsx(markers = pres.sort, file = paste0(outdir,"presto_res.xlsx"), group.by = "group", n = 100)
    
}


fm_report <- function(res, sobj, outdir){
	# summary of n diff genes per cluster, up down, combined 
	### we want a barchart that is grouped (up, down) for each cluster celltype
	### to do this begin by iterating through results
	### testing whether a feature is significantly up or down regulated, then adding that status as a new column
	dir.create(outdir, showWarnings = FALSE)
	de.status <- c()
	for (r in seq(1, nrow(res))){
		trow <- res[r,]
		if (trow$avg_log2FC > .58 & trow$p_val_adj < .05){
			de.status <- append(de.status, "up")
		} else if (trow$avg_log2FC < -.58 & trow$p_val_adj < .05){
			de.status <- append(de.status, "down")
		} else {
			de.status <- append(de.status, "notDE")
		}
	}
	res$status <- de.status
	# write res table to file
#	tbl <- table(res$cluster, res$status)[,c(1,3)]
#	write.table(tbl, paste0(outdir, "fm_summary.txt"), sep = "    ",
#	    quote = FALSE, row.names = TRUE)

	### barplot 
	res.de <- res[res$status == "up" | res$status == "down",]
	res.t <- table(res.de$cluster, res.de$status)
	res.t.m <- melt(res.t)

	colnames(res.t.m) <- c("cluster", "status", "counts")
	res.t.m$cluster <- as.factor(res.t.m$cluster)

	bp <- ggplot(res.t.m, aes(fill=status, y=counts, x=cluster)) + 
		geom_bar(position="dodge", stat="identity")  +
		theme(axis.text.x = element_text(angle = 30, hjust = 1)) + 
		ggtitle("# DEG per cluster/celltype") + 
		theme(plot.title = element_text(size = 20, face = "bold") + 
		scale_x_continuous(breaks = seq(1,2))
        )
	bp <- bp + scale_fill_manual(values=c("up" = "darkseagreen3", "down" = "cadetblue2"))
	ggsave(filename = paste0("fm_barplot.pdf"), plot = bp, device = "pdf", 
	       path = outdir, width = 12)


	# Filter res to top 10 marker per cluster
	# aleady sorted by padj
	res <- na.omit(res)
	top12 <- res %>% group_by(cluster) %>% slice_head(n = 12)

	### heatmap 
	# res %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC) -> top10
	#hmap <- DoHeatmap(sobj, features = top10$gene, slot = "scale.data") + NoLegend()
	
	# iterate through clusters and create Feature plots of top 10 markers 
	dir.create(paste0(outdir, "featureplots/"), showWarnings = FALSE)
	for (c in unique(sobj$seurat_clusters)){ 
  	    top12.sub <- top12[top12$cluster == c,]
	    if (nrow(top12.sub) > 0){
   	        fp <- suppressMessages(FeaturePlot(sobj, features = unique(top12.sub$gene), cols = c("lightblue", "darkblue"), 
			raster=TRUE, order = TRUE))# +#plot_layout(ncol = 4) 
 	        suppressMessages(ggsave(filename = paste0("c", c, "_fm_featureplots.pdf"), plot = fp, device = "pdf", 
		    path = paste0(outdir, "featureplots/"), width = 24, height = 20))
	    }
	}

	# DotPlot of top 10
	dp <- DotPlot(object = sobj, assay = "RNA", features = unique(top12$gene), 
		cols = c("lightblue", "darkblue"), cluster.idents = TRUE) + coord_flip()
	suppressMessages(ggsave(filename = paste0("fm_dotplot.pdf"), plot = dp, device = "pdf", 
		path = paste0(outdir), width = 14, height = 30))
	    
	# save nebulosa plots
	#pdf(file = paste0(outdir, "plots/ct_", c, "_nebulosa.pdf"), width = 12, height = 12)
	#    print(np)
	#dev.off()
	# save dotplot
	write.top.n.xlsx(markers = res, file = paste0(outdir,"fm_res.xlsx"), group.by = "cluster", n = 100)


	return(res.de)
}


fm_cProfiler<- function(res.de, outdir, organism = "human", universe = NULL){
	# run cluster profiler on the up/down from each cluster
	dir.create(paste0(outdir, "CP_out/"), showWarnings = FALSE)
	dir.create(paste0(outdir, "CP_out/tables/"), showWarnings = FALSE)
	dir.create(paste0(outdir, "CP_out/plots"), showWarnings = FALSE)

	de.status <- c()
	for (r in seq(1, nrow(res.de))){
		trow <- res.de[r,]
		if (trow$avg_log2FC > .58 & trow$p_val_adj < .05){
			de.status <- append(de.status, "up")
		} else if (trow$avg_log2FC < -.58 & trow$p_val_adj < .05){
			de.status <- append(de.status, "down")
		} else {
			de.status <- append(de.status, "notDE")
		}
	}
	res.de$status <- de.status

	for (c in seq(0, length(unique(res.de$cluster)))){
	    res.cp <- res.de[res.de$cluster == c & res.de$status == "up",]
	    
	    print(paste0("N input genes for cluster ", c, ":  ", nrow(res.cp))) 
	    if (nrow(res.cp) > 200){res.cp <- res.cp[1:200,]}
	    
	    # handle organism here
	    if(organism == "human"){
		    org.db <- org.Hs.eg.db
	    } else if(organism == "mouse"){
		    org.db <- org.Mm.eg.db
	    }


	    ego_bp <- suppressMessages(
			enrichGO(gene = res.cp$gene,
                    OrgDb = org.db,
                    keyType = 'SYMBOL',
                    universe = universe, 
                    ont = "BP",
                    pAdjustMethod = "BH",
                    pvalueCutoff = 1,
                    qvalueCutoff = 1))

	    if (!is.null(ego_bp)){
	        write.table(ego_bp@result, file = paste0(outdir,"CP_out/tables/ct_", c, "_CP_BP_go.csv"))
			write_tableHTML(tableHTML(ego_bp@result), file = 
		    paste0(outdir,"CP_out/tables/ct_", c, "_CP_BP_go.html"))

			p2 <- dotplot(ego_bp, showCategory=20) + ggtitle("biological process")
			fn2<- paste0("ct_", c, "_CP_BP_DTPLT.pdf")
			suppressMessages(ggsave(filename = fn2, plot = p2, device = "pdf", 
					path = paste0(outdir, "CP_out/plots/"),  height = 10))
	    }
	}
}


extract_detected_features <- function(sobj){
    sums <- rowSums(sobj@assays$RNA@counts)
    detected_genes <- names(sums[sums > 0])

    return(detected_genes)
}

# Change layering of feature plots
# Change color scale of feature plots 

