### Create scCustomize grouped top 10 dotplot
cDotplot <- function(res, sobj, nFeatures, k, outdir, width = 10, height = 16){

    topN_markers <- Extract_Top_Markers(marker_dataframe = res, num_genes = nFeatures, named_vector = FALSE,
         make_unique = TRUE)

    plot <- Clustered_DotPlot(seurat_object = sobj, features = topN_markers, k = k)

    pdf(file = paste0(outdir, "FM_cDotplot.pdf"), width = width, height = height)
    print(plot)
    dev.off()

#    return(plot)
}
