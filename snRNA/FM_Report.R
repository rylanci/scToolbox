suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(future))
suppressPackageStartupMessages(library(tableHTML))
suppressPackageStartupMessages(library(openxlsx))
suppressPackageStartupMessages(library(clusterProfiler))
suppressPackageStartupMessages(library(org.Hs.eg.db))
suppressPackageStartupMessages(library(org.Mm.eg.db))
### This script will contain functions to run findmarkers and write a report on the output

# Functions
run_fm <- function(sobj, n.workers = 1, idents = "seurat_clusters", logfc.threshold = 0.25, sct.adjust = FALSE){
	plan(multicore, workers = n.workers)
	options(future.globals.maxSize = 2000 * 1024^2)

	Idents(sobj) <- idents
	res <- FindAllMarkers(sobj, logfc.threshold = logfc.threshold)

	return(res)
}

fm_report <- function(res, outdir){
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
	tbl <- table(res$cluster, res$status)[,c(1,3)]
	write.table(tbl, paste0(outdir, "fm_de_table.tsv"), sep = "    ",
	    quote = FALSE, row.names = rownames(tbl))

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
	ggsave(filename = paste0("FM_barplot.pdf"), plot = bp, device = "pdf", 
	       path = outdir, width = 12)

	### heatmap 
	# res %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC) -> top10
	#hmap <- DoHeatmap(sobj, features = top10$gene, slot = "scale.data") + NoLegend()
	
	dir.create(paste0(outdir, "plots/"), showWarnings = FALSE)
	# nebulosa/ FeaturePlot and dotplot for top n gene per cluster
	for (c in seq(0, max(res$cluster))){
		# subset res for target cluster
		res.c <- na.omit(res[res$cluster == c,])
		# order subset by p.adj
		res.o <- res.c[order(res.c$p_val_adj),]
		# take top 12
		genes <- res.o$gene

		# create nebulosa plots 
		#np <- plot_density(sobj, features = genes[1:9])
		#np + plot_layout(ncol = 3)
		# create dotplot and feature plot
#		print(genes)
		if (length(genes) > 12){
		    dp <- DotPlot(sobj, features = genes[1:12]) + ggtitle(paste0("cluster ", c, " top 12 DEG"))
		    dp + theme(axis.text.x = element_text(angle = 30, hjust = 1))
		    fp <- suppressMessages(FeaturePlot(sobj, features = genes[1:12]) + ggtitle(paste0("cluster ", c, " top 12 DEG")))
		} else {
		    dp <- DotPlot(sobj, features = genes) + ggtitle(paste0("cluster ", c, " top 12 DEG"))
		    dp + theme(axis.text.x = element_text(angle = 30, hjust = 1))
		    fp <- suppressMessages(FeaturePlot(sobj, features = genes) + ggtitle(paste0("cluster ", c, " top 12 DEG")))
		}
		# save nebulosa plots
		#pdf(file = paste0(outdir, "plots/ct_", c, "_nebulosa.pdf"), width = 12, height = 12)
		#    print(np)
		#dev.off()
		# save dotplot
		suppressMessages(ggsave(filename = paste0("ct_", c, "_FM_dotplot.pdf"), plot = dp, device = "pdf", 
			path = paste0(outdir, "plots/"), width = 12))
		suppressMessages(ggsave(filename = paste0("ct_", c, "_FM_featureplot.pdf"), plot = fp, device = "pdf", 
 			path = paste0(outdir, "plots/"), width = 18, height = 18))

	}

	return(res.de)
}


fm_cProfiler<- function(res.de, outdir, organism = "human"){
	# run cluster profiler on the up/down from each cluster
	# top 200
	dir.create(paste0(outdir, "CP_out/"), showWarnings = FALSE)
	dir.create(paste0(outdir, "CP_out/tables/"), showWarnings = FALSE)
	dir.create(paste0(outdir, "CP_out/plots"), showWarnings = FALSE)

	for (c in seq(0, length(unique(res.de$cluster)))){
	    res.cp <- res.de[res.de$cluster == c & res.de$status == "up",]
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
                    universe = NULL, 
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

# Writes spreadsheet of top DE genes
write.top.n.xlsx <- function(markers, outdir, group.by = "cluster", n = 100){
    wb <- createWorkbook("TopMarkers")

    for (c in unique(markers[[group.by]])){
	if (c == 0){c <- "z"}
        addWorksheet(wb, c)
        tdf <- head(markers[markers[[group.by]] == c, ], n = n)
        writeData(wb, sheet = c, x = tdf)
    }

    saveWorkbook(wb = wb, file = outdir, overwrite = TRUE)
}


#fm_kmeans <- function(){
	# think about this later
#}
