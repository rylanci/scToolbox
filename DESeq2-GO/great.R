suppressPackageStartupMessages(library(rGREAT))
suppressPackageStartupMessages(library(GenomicRanges))
library(stringr)


res_to_bed <- function(df){
	peaks <- rownames(df)
	p_split <- str_split(peaks, "-")

	chr <- c()
	start <- c()
	end <- c()
	for (p in p_split){
		chr <- append(chr, p[1])
		start <- append(start, p[2])
		end <- append(end, p[3])
	}

	peak.df <- data.frame("chr" = chr, "start" = start, "end" = end)
	return(peak.df)

}

run_rGreat <- function(res, genome = "hg38", outdir){

	sig.peaks <- read.table(res, header = TRUE, row.names = 1)

	sig.up <- sig.peaks[sig.peaks$log2FoldChange > 0,]
	sig.dwn<- sig.peaks[sig.peaks$log2FoldChange < 0,]

	print(nrow(sig.up))
	print(nrow(sig.dwn))
	
	dir.create(paste0(outdir, "rGreat_out/"), showWarnings = FALSE)	
	bn <- basename(res)
	comp <- str_remove(bn, "Sig_")
	comp <- str_remove(bn, ".csv")


	if (nrow(sig.up) > 0){
		peak.df.up <- res_to_bed(sig.up)
		gr.up <- makeGRangesFromDataFrame(peak.df.up)

		gout_bp.up <- great(gr.up, "GO:BP", genome)
		gout_cc.up <- great(gr.up, "GO:CC", genome)
		gout_mf.up <- great(gr.up, "GO:MF", genome)

		tb_bp.up <- getEnrichmentTable(gout_bp.up)
		tb_cc.up <- getEnrichmentTable(gout_cc.up)
		tb_mf.up <- getEnrichmentTable(gout_mf.up)

		write.csv(tb_bp.up ,paste0(outdir, "/rGreat_out/RG_", comp, "_bp_up.csv"))
		write.csv(tb_cc.up ,paste0(outdir, "/rGreat_out/RG_", comp, "_cc_up.csv"))
		write.csv(tb_mf.up ,paste0(outdir, "/rGreat_out/RG_", comp, "_mf_up.csv"))
	}

	if (nrow(sig.dwn) > 0){
		peak.df.dwn <- res_to_bed(sig.dwn)
		gr.dwn <- makeGRangesFromDataFrame(peak.df.dwn)

		gout_bp.dwn <- great(gr.dwn, "GO:BP", genome)
		gout_cc.dwn <- great(gr.dwn, "GO:CC", genome)
		gout_mf.dwn <- great(gr.dwn, "GO:MF", genome)

		tb_bp.dwn <- getEnrichmentTable(gout_bp.dwn)
		tb_cc.dwn <- getEnrichmentTable(gout_cc.dwn)
		tb_mf.dwn <- getEnrichmentTable(gout_mf.dwn)

		write.csv(tb_bp.up ,paste0(outdir, "/rGreat_out/RG_", comp, "_bp_dwn.csv"))
		write.csv(tb_cc.up ,paste0(outdir, "/rGreat_out/RG_", comp, "_cc_dwn.csv"))
		write.csv(tb_mf.up ,paste0(outdir, "/rGreat_out/RG_", comp, "_mf_dwn.csv"))
	}

}


iterate_rGreat <- function(path_to_res, outdir){
	
	# to do: 
		# initialize log
		# have run_rGreat append info 

	res.files <- list.files(path_to_res)
	res.paths <- paste0(path_to_res, res.files)

	print(res.paths)
	
	lapply(X = res.paths, FUN = run_rGreat, outdir =  outdir)

	
}		


#res <- "/tscc/nfs/home/rlancione/myprojects/Liver_FNIH/Callpeaks/DESeq_Results_n20_ctpeaks/Significant_Results/Sig_Hepatocytes_HF_NASH-LF_NASH.csv"
res <- "/tscc/nfs/home/rlancione/myprojects/Liver_FNIH/Callpeaks/DESeq_Results_n20_ctpeaks/Significant_Results/Sig_TCells_HF_NASH-Normal.csv"
#run_rGreat(res = res, outdir = "/tscc/nfs/home/rlancione/myprojects/sandbox/Great/test_out/")


ptr <- "/tscc/nfs/home/rlancione/myprojects/Liver_FNIH/Callpeaks/DESeq_Results_n20_ctpeaks/Significant_Results/"
out <- "/tscc/nfs/home/rlancione/myprojects/sandbox/Great/test_out/"
iterate_rGreat(path_to_res = ptr, outdir = out)


