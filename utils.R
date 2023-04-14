#library(Seurat, quietly = T)
#library(Signac, quietly = T)
#library(GenomicRanges, quietly = T)
#library(dplyr, quietly = T)
#library(ggplot2, quietly = T)
#library(Matrix, quietly = T)
#library(preprocessCore)
#library(harmony)
#library(openxlsx)
#library(chromVAR, quietly = TRUE)
#library(motifmatchr, quietly = TRUE)
#library(SummarizedExperiment, quietly = TRUE)
#library(JASPAR2020)
#library(TFBSTools)
#library(BSgenome.Mmusculus.UCSC.mm10)
#library(BSgenome.Hsapiens.UCSC.hg38)
#library(EnsDb.Mmusculus.v79)
#library(EnsDb.Hsapiens.v86)

# Read seurat objects. Returns list of objects
# add functionality to take 2d path
read_objects <- function(ipath1, ipath2=NULL){
    objs <- list.files(ipath1)
    obj_list <- c()

    for(o in objs){
        temp_obj <- readRDS(paste0(ipath1,o))
        obj_list <- append(obj_list, temp_obj)
    }

    if (is.null(ipath2) == FALSE){
        objs2 <- list.files(ipath2)
        for(o in objs2){
        temp_obj2 <- readRDS(paste0(ipath2,o))
        obj_list <- append(obj_list, temp_obj2)
        }
    }
    return(obj_list)
}

# Writes spreadsheet of top DE genes
write.top.n.xlsx <- function(markers, outdir, group.by = "cluster", n = 100){
    wb <- createWorkbook("TopMarkers")

    for (c in unique(markers[[group.by]])){
        addWorksheet(wb, c)
        tdf <- head(markers[markers[[group.by]] == c, ], n = n)
        writeData(wb, sheet = c, x = tdf)
    }

    saveWorkbook(wb = wb, file = outdir, overwrite = TRUE)
}


peak_calling <- function(sobj, output_assay, species, group.by = "seurat_clusters", outdir = tempdir()){
  print("Calling Peaks")
  ######## Changed  to group.by experiment ###########
  peaks <- CallPeaks(object = sobj, group.by = group.by, outdir = outdir, cleanup = FALSE)

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

chromvar <- function(sobj, species){
    pfm <- getMatrixSet(
        x = JASPAR2020,
        opts = list(collection = "CORE", tax_group = 'vertebrates', all_versions = FALSE)
        )
    if (species == "mouse"){
        # add motif information
        sobj <- AddMotifs(
                object = sobj,
                genome = BSgenome.Mmusculus.UCSC.mm10,
                pfm = pfm)
        sobj <- RunChromVAR(
                object = sobj,
                genome = BSgenome.Mmusculus.UCSC.mm10)
    } else if (species == "human"){
        sobj <- AddMotifs(
                object = sobj,
                genome = BSgenome.Hsapiens.UCSC.hg38,
                pfm = pfm)
        sobj <- RunChromVAR(
                object = sobj,
                genome = BSgenome.Hsapiens.UCSC.hg38)
    } else {print("incorrect species input for chromvar: choose 'mouse' or 'human'")}

    return(sobj)
}

activity <- function(atac.sobj, species){
    if (species == "mouse"){
        annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Mmusculus.v79)
        seqlevelsStyle(annotations) <- "UCSC"
        genome(annotations) <- "mm10"
        Annotation(atac.sobj) <- annotations
    } else if (species == "human"){
        print("human data is not currently supported. I forget why. Probably cause Ryan is lazy")
    }
    gene.activities <- GeneActivity(atac.sobj)
    # add gene activities as a new assay
    atac.sobj[["ACTIVITY"]] <- CreateAssayObject(counts = gene.activities)
    # normalize gene activities
    DefaultAssay(atac.sobj) <- "ACTIVITY"
    atac.sobj <- NormalizeData(atac.sobj)
    atac.sobj <- ScaleData(atac.sobj, features = rownames(atac.sobj))

    return(atac.sobj)
}

annotation <- function(label_transfer=FALSE, atac.sobj, rna.sobj, anno_column){
     if (label_transfer == TRUE){
        gene.activities <- GeneActivity(atac.sobj, features = VariableFeatures(rna.sobj))
        transfer.anchors <- FindTransferAnchors(reference = rna.sobj, query = atac.sobj, features = VariableFeatures(object = rna.sobj),
            reference.assay = "RNA", query.assay = "ACTIVITY", reduction = "cca")

        atac_labels <- TransferData(anchorset = transfer.anchors, refdata = rna.sobj[[anno_column]], weight.reduction = mwt[["lsi"]],
            dims = 2:30)

        atac.sobj <- AddMetaData(atac.sobj, atac_labels, col.name = NULL)
        return(atac.sobj)
    } else if(label_transfer == FALSE) {
         print("skipping label transfer")
    }
}

norm_cellquant_bplot <- function(sobj, group.by = "conditions", invert = FALSE, x.scale = 2){
    if (invert == FALSE){
        Idents(sobj) <- "seurat_clusters"
        #table <- table(Idents(sobj), sobj$conditions)
	table <- table(Idents(sobj), sobj@meta.data[[group.by]])
        nmax <- max(table(Idents(sobj)))
        #norm_df <- data.frame(row.names = seq(0, length(unique(sobj$conditions))))
	norm_df <- data.frame(row.names = seq(0, length(unique(sobj@meta.data[[group.by]]))))
    } else if (invert == TRUE){ # invert conditions to x-axis. Cluster contribution as y.
        Idents(sobj) <- "conditions"
        table <- table(Idents(sobj), sobj$seurat_clusters)
        nmax <- max(table(Idents(sobj)))
        norm_df <- data.frame(row.names = seq(0, length(unique(sobj$seurat_clusters))))
    }

    # Buld normalized cell quantity table
    for (r in seq(1,nrow(table))){
        #print(r)
        rsum <- sum(table[r,])
        cfactor <- nmax / rsum
        trow <- table[r,] * cfactor
        trow <- trow / sum(trow) * 100
        norm_df <- rbind(norm_df, trow)
    }

    colnames(norm_df) <- colnames(table)
    row.names(norm_df ) <- seq(1, nrow(norm_df))
    #norm_df
    cols <- c("cadetblue3", "coral3", "darkolivegreen3", "darkcyan", "mediumpurple1",
          "lightgreen", "lightgoldenrod", "lightslateblue", "mistyrose", "lightblue4",
          "navajowhite1", "magenta", "lightsalmon", "mediumorchid1", "midnightblue",
         "lightskyblue", "lightgoldenrodyellow", "black", "lightgrey", "mistyrose4")

    ncols <- sample(x = cols, size = ncol(norm_df), replace = F)
    #ncols
   # pdf(file = "Plots/ATAC_all_barplot.pdf", width = 14, height = 10)
    #barplot(height = as.matrix(t(norm_df)), col = cols, legend = T, , xlim = c(0, 12),
    #              names.arg = seq(0,length(unique(sobj$seurat_clusters)) -1))
   # dev.off()
     if (invert == FALSE){
         return(barplot(height = as.matrix(t(norm_df)), col = cols, legend = T,
                    xlim = c(0, (length(unique(sobj$seurat_clusters)) * x.scale)),
                    names.arg = seq(0,length(unique(sobj$seurat_clusters)) -1)))
         } else if(invert == TRUE){
         return(barplot(height = as.matrix(t(norm_df)), col = cols, legend = T,
                    xlim = c(0, (length(unique(sobj$conditions)) * x.scale)),
                    names.arg = row.names(table)))
     }
}


# Function to add conditions using a csv
add.conditions <- function(sobj, dset.col = "orig.ident", cond.csv.path){
    # read 2 col csv. C1 = dset. C2 = condition
    csv <- read.csv2(cond.csv.path, sep = ",", header = F)

    # Create function to return condition for each dset by referencing csv
    cond_list <- function(dataset){
        csv.row <- csv[csv$V1 == dataset,]
        condition <- csv.row[,2]
        return(condition)
        }

    # Apply func to each row of metadata dset col.
    conditions <- apply(X = as.data.frame(sobj@meta.data[[dset.col]]), MARGIN = 1, FUN = cond_list)
    # Store new vector of conditions in object
    sobj$conditions <- conditions
    # Return new object with condition column
    return(sobj)
}


qt.normalize <- function(object_list){
        # Find max ncells out of object list. Need that to create correct size matrix max ncells will be nrows
        # matrix is max ncells by n objects in list
    ncells.l <- sapply(object_list, FUN = function(object){
        return(nrow(object@meta.data))
        })
    max.c <- max(ncells.l)

    # initialize matrix
    dmat <- matrix(data = as.numeric("NA"), nrow = max.c, ncol = length(object_list))
    # for o in obj list: get scores, append to matrix
    cntr <- 1
    for (o in object_list){
        # extract doublet scores
        d.scores <- o@meta.data$DF.scores
        # enter scores into matrix to be normalized
        dmat[1:length(d.scores),cntr] <- d.scores

        cntr <- cntr + 1
    }
    # quantile normalize matrix
    n.dmat <- normalize.quantiles(x = dmat)

    #print(head(n.dmat))
    # for o in obj list: add col of matrix to metadata, replace itself in object list
    cntr <- 1
    object_list2 <- c()
    for (o in object_list){
        # extract normalized scores
        n.scores <- n.dmat[,cntr]
        # apply them to objects with NA's stripped
        o@meta.data$qn.DF.scores <- na.omit(n.scores)
        # add object with normalized scores to new list
        object_list2[[cntr]] <- o
        cntr <- cntr + 1
    }
    return(object_list2)
}


qt.norm.drm.scoring <- function(object, normalize = TRUE, outdir = NULL, dset_col = "orig.ident"){
# If single object is submitted with no qt normalization. Run All
    if(length(object) == 1 & normalize == TRUE){
        # 1.Subset obj by dataset
        # 2.create list of objects
        #Idents(object) <- "orig.ident"
        #dsets <- unique(object@meta.data$orig.ident)
        # changed for atac flexibility
        Idents(object) <- dset_col
        dsets <- unique(Idents(object))

        obj_list <- c()
        for (d in dsets){
            t.obj <- subset(object, ident = d)
            obj_list <- append(obj_list, t.obj)
        }

        # 3.apply qt.normalize function
        n.obj.list <- qt.normalize(obj_list)

        # 4.merge object list with norm doub scores
        mobj <- merge(x= n.obj.list[[1]],
            y= n.obj.list[2:length(n.obj.list)])

    } else if(length(object) > 1 & normalize == TRUE){
        # 3.apply qt.normalize function
        n.obj.list <- qt.normalize(obj_list)

        # 4.merge object list with norm doub scores
        mobj <- merge(x= n.obj.list[[1]],
            y= n.obj.list[2:length(n.obj.list)])

    } else if(length(object) > 1 & normalize == FALSE){
        # 4.merge object list with norm doub scores
        mobj <- merge(x= n.obj.list[[1]],
            y= n.obj.list[2:length(n.obj.list)])

    }

# If object length = 1 and norm = False. Start here
    # 5.For each cluster extract the norm doublet scores
    Idents(mobj) <- "seurat_clusters"
    clusters <- seq(0, length(unique(mobj@meta.data$seurat_clusters)) - 1 )
    clust.qn.dfscores <- lapply(X = clusters, FUN = function(c){
        # For each cluster; do this
        t.obj <- subset(mobj, idents = c)
        scores <- t.obj$qn.DF.scores
        return(scores)
        })

    # 6.Take the summary() of each distribution of doub scores
    # 7.Print summary of clust distributions to a file
    cntr <- 1
    doub.table <- data.frame()
    for (c in clust.qn.dfscores){
        s <- summary(c)
        # Row is clust #, ncells, summary scores(min, q1, median, mean, q3, max)
        t.row <- c(cntr - 1, length(c), s)
        doub.table <- rbind(doub.table, t.row)
        cntr <- cntr + 1
    }
    colnames(doub.table) <- c("Cluster", "nCells", "MinScore", "Q1", "Median", "Mean", "Q3", "MaxScore")
    dt.sort <- doub.table[order(doub.table$Mean),]


    # transfer qn.scores to original object. Create feature plot of qn.scores on orig embedding
    object[["qn.DF.scores"]] <- mobj$qn.DF.scores
    # ****** Plotting ******
    # 8.Create violin plots,  kelseys line plot, write doub table to file
    if (is.null(outdir) == FALSE){
        # Save doub.table as tsv
        write.table(x = dt.sort, file = paste0(outdir, "normalized_doublet_summary.txt"), row.names = F)

        # Boxplots
        bxp1 <- ggplot(mobj@meta.data, aes(seurat_clusters, DF.scores)) +
          geom_boxplot(varwidth=T, fill="plum") +
          labs(title="Liver RNA Box plot",
          subtitle="Doublet Score by Cluster",
          x="Cluster",
          y="Doublet Score")

        bxp2 <- ggplot(mobj@meta.data, aes(seurat_clusters, qn.DF.scores)) +
          geom_boxplot(varwidth=T, fill="plum") +
          labs(title="Liver RNA Box plot",
          subtitle="Normalized Doublet Score by Cluster",
          x="Cluster",
          y="Normalized Doublet Score")
        # Kelseys line plot
        min.mean <- min(dt.sort$Mean)
        max.mean <- max(dt.sort$Mean)
        klp1 <- ggplot(dt.sort, aes(x=1:nrow(dt.sort), y=Mean, label=Cluster)) +
            geom_point() + geom_line() +
            theme(axis.text.x=element_blank(),
                  axis.ticks.x=element_blank()) +
            ylab("Mean Normalized Doublet Score") + xlab("Clusters") +
            geom_text(hjust=-0.2, vjust=-2) +
            lims(x= c(1,nrow(dt.sort)),y=c(min.mean - .01, max.mean + .02)) +
            ggtitle("Mean Normalized Doublet Score per Cluster")
            #geom_text_repel(aes(x=1:nrow(dt.sort), y=Mean, label=Cluster))

      dt.sort.2 <- doub.table[order(doub.table$Median),]
      min.median <- min(dt.sort$Median)
      max.median <- max(dt.sort$Median)
      klp2 <- ggplot(dt.sort.2, aes(x=1:nrow(dt.sort.2), y=Median, label=Cluster)) +
          geom_point() +
          geom_line() +
          theme(axis.text.x=element_blank(),
              axis.ticks.x=element_blank()) +
            ylab("Median Normalized Doublet Score") +
            xlab("Clusters") +
            geom_text(hjust=0, vjust=-2) +
            lims(x= c(1,nrow(dt.sort.2)),y=c(min.median - .01, max.median + .02)) +
            ggtitle("Median Normalized Doublet Score per Cluster")


        Idents(object) <- "seurat_clusters"
        Idents(mobj) <- "seurat_clusters"
        vp1 <- VlnPlot(object, features = c("DF.scores"), ncol = 1)
        vp2 <- VlnPlot(mobj, features = c("qn.DF.scores"), ncol = 1)
        Idents(object) <- dset_col
        Idents(mobj) <- dset_col
        vp3 <- VlnPlot(object, features = c("DF.scores"), ncol = 1)
        vp4 <- VlnPlot(mobj, features = c("qn.DF.scores"), ncol = 1)

        # Umap of high res object
        Idents(object) <- "seurat_clusters"
        ump1 <- DimPlot(object, label = T)
        # Feature plot of normalized scores on original object
        fplt1 <- FeaturePlot(object, features = "qn.DF.scores")


        # Write plots to outdir
        ggsave(filename = paste0(outdir, "mean_ndoub_lplot.pdf"), plot = klp1, device = "pdf", width = 13, height = 9)
        ggsave(filename = paste0(outdir, "median_ndoub_ndoub_lplot.pdf"), plot = klp2, device = "pdf", width = 13, height = 9)
        ggsave(filename = paste0(outdir, "doub_bxplt.pdf"), plot = bxp1, device = "pdf", width = 20, height = 9)
        ggsave(filename = paste0(outdir, "ndoub_bxplt.pdf"), plot = bxp2, device = "pdf", width = 20, height = 9)
        ggsave(filename = paste0(outdir, "doub_clust_vln.pdf"), plot = vp1, device = "pdf", width = 14, height = 9)
        ggsave(filename = paste0(outdir, "ndoub__clust_vln.pdf"), plot = vp2, device = "pdf", width = 14, height = 9)
        ggsave(filename = paste0(outdir, "doub_dset_vln.pdf"), plot = vp3, device = "pdf", width = 14, height = 9)
        ggsave(filename = paste0(outdir, "ndoub__dset_vln.pdf"), plot = vp4, device = "pdf", width = 14, height = 9)
        ggsave(filename = paste0(outdir, "ndoub__dset_vln.pdf"), plot = vp4, device = "pdf", width = 14, height = 9)
        ggsave(filename = paste0(outdir, "high_res_umap.pdf"), plot = ump1, device = "pdf", width = 11, height = 9)
        ggsave(filename = paste0(outdir, "ndoub_fplt.pdf"), plot = fplt1, device = "pdf", width = 11, height = 9)

    }

    #return(mobj)
    return(object)
}


sct.harm.processing <- function(sobj, dims=1:20, res=.05, n.neigh = 30L, min.dist = .3, spread = 1, harmony = FALSE, CM = TRUE){
    sobj <- SCTransform(sobj, vst.flavor = "v2",
		conserve.memory = CM, verbose = FALSE) %>%
            	RunPCA(verbose = FALSE)

    if(harmony){
      sobj <- RunHarmony(object = sobj,
                reduction = "pca",
                group.by.vars = "orig.ident",
                assay.use = 'SCT',
                project.dim = FALSE)

      sobj <- FindNeighbors(sobj, dims = dims, verbose = FALSE, reduction = "harmony") %>%
              FindClusters(resolution = res, verbose = FALSE) %>%
              RunUMAP(dims = dims, n.neighbors = n.neigh, min.dist = min.dist, spread = spread, verbose = FALSE, reduction = "harmony")
    } else {
      sobj <- FindNeighbors(sobj, dims = dims, verbose = FALSE, reduction = "pca") %>%
              FindClusters(resolution = res, verbose = FALSE) %>%
              RunUMAP(dims = dims, n.neighbors = n.neigh, min.dist = min.dist, spread = spread, verbose = FALSE, reduction = "pca")
    }

    return(sobj)
}


lsi.processing <- function(sobj, TF.method = 1, dims = 1:30, res = .05, assay = "ATAC", harmony = FALSE){
    DefaultAssay(sobj) <- assay
    sobj <- RunTFIDF(sobj, method = TF.method, assay = assay) %>%
                FindTopFeatures(min.cutoff = 'q0') %>%
                RunSVD()
    if(harmony){
        sobj <- RunHarmony(object = sobj,
                    reduction = "lsi",
                    group.by.vars = "experiment",
                    assay.use = assay,
                    project.dim = FALSE)

        sobj <- FindNeighbors(sobj, dims = dims, verbose = FALSE, reduction = "harmony") %>%
                FindClusters(resolution = res, verbose = FALSE, algorithm = 3) %>%
                RunUMAP(dims = dims , verbose = FALSE, reduction = "harmony")
    } else{
        sobj <- FindNeighbors(sobj, dims = dims, verbose = FALSE, reduction = "lsi") %>%
                FindClusters(resolution = res, verbose = FALSE, algorithm = 3) %>%
                RunUMAP(dims = dims , verbose = FALSE, reduction = "lsi")
    }

    return(sobj)
}
