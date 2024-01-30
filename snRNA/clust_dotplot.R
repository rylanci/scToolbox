library(Seurat)
library(scCustomize)
library(argparse)

process_args <- function(){
      # create parser object
      parser <- ArgumentParser(description = "Create seurat object from 10x h5 and collect QC metrics")

      # Create parser arg group for our required args
      required_arg_group = parser$add_argument_group('flagged required arguments',
          'the script will fail if these args are not included')

      ### start args
      required_arg_group$add_argument("-o","--output_path", required = TRUE,
          help="path to output location")

      required_arg_group$add_argument("-n","--name", required = TRUE,
          help="name of analysis e.g. batch_harmony")


      # optional args
#      parser$add_argument("-w","--workers",
#          help="number of cores to run FindMarkers with",
#          default = 12)

      args <- parser$parse_args()
      return(args)
  }


### Create scCustomize grouped top 10 dotplot
cDotplot <- function(res, sobj, nFeatures, k, outdir, width = 12, height = 20){

      topN_markers <- Extract_Top_Markers(marker_dataframe = res, num_genes = nFeatures, named_vector = FALSE,
           make_unique = TRUE)

      plot <- Clustered_DotPlot(seurat_object = sobj, features = topN_markers, k = k)

      pdf(file = paste0(outdir, "FM_cDotplot.pdf"), width = width, height = height)
      print(plot[2])
      dev.off()

  #    return(plot)
  }


#### create scCustomize clustered dotplot
args <- process_args()

print("Initialize scCustomize Dotplot Report")
# Identify processed object
files <- list.files(paste0(args$output_path, "Objects/"))
file <- grep(pattern = paste0(args$name,".RDS"), x = files, value = TRUE)

print(paste0("Plotting: ", args$output_path, "Objects/", file))
sobj <- readRDS(paste0(args$output_path, "Objects/", file))
res <- read.table(paste0(args$output_path,  "FM_output/", args$name, "/fm_res.txt"))

cDotplot(res = res, sobj = sobj, nFeatures = 10, outdir = paste0(args$output_path, "FM_output/", args$name, "/"), height = 20, width = 12)


## Scale height by length 0f topNmarkers 
