library(Matrix)
library(ggplot2)
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(ATACseqQC))
suppressPackageStartupMessages(library(Rsamtools))
suppressPackageStartupMessages(library(ballgown))

args <- commandArgs(trailingOnly = TRUE)
bam <- args[1] 
genome <- args[2]
name<- args[3]

#genome <- "rat"
#bam <- "/projects/ps-epigen/users/rlan/ATACSeq/zemke_scripts/Rat_Samples/Derek_QY_2054-2077//mapping//QY_2055_rat_sorted_rmdup.bam"

if (genome == "mm10") {
	cat("calculated TSS enrichment with mm10 reference...\n")
	suppressPackageStartupMessages(library(TxDb.Mmusculus.UCSC.mm10.knownGene))
	seqinformation <- seqinfo(TxDb.Mmusculus.UCSC.mm10.knownGene)
#	major_chr <- read.table("/projects/ps-renlab/y2xie/projects/genome_ref/mm10.main.chrom.sizes")
	txs <- transcripts(TxDb.Mmusculus.UCSC.mm10.knownGene)
} else if (genome == "human") {
	cat("calculated TSS enrichment with hg38 reference...\n")
	suppressPackageStartupMessages(library(TxDb.Hsapiens.UCSC.hg38.knownGene))
	seqinformation <- seqinfo(TxDb.Hsapiens.UCSC.hg38.knownGene)
#	major_chr <- read.table("/projects/ps-epigen/GENOME/hg38/hg38.chrom.sizes")
	txs <- transcripts(TxDb.Hsapiens.UCSC.hg38.knownGene)
} else if (genome == "cow") {
        cat("calculated TSS enrichment with cow reference...\n")
# 	 not using these because we need a scaffold based set of transcripts... which this database is not 
#        suppressPackageStartupMessages(library(TxDb.Btaurus.UCSC.bosTau8.refGene))
#        seqinformation <- seqinfo(TxDb.Btaurus.UCSC.bosTau8.refGene)
#        major_chr <- read.table("/projects/ps-epigen/users/rlan/ATACSeq/zemke_scripts/bt9.chrom.sizes")
#        txs <- transcripts(TxDb.Btaurus.UCSC.bosTau8.refGene)
	
	# scaffold based gtf paired with the bowtie reference  
	gtf <- "/projects/ps-epigen/users/rlan/ATACSeq/zemke_scripts/ARS-UCD1.3.gtf"
	gtf <- gffReadGR(gtf)
	txs <- gtf[gtf$type == "exon",]
} else if (genome == "rat") {
        cat("calculated TSS enrichment with rat reference...\n")
        suppressPackageStartupMessages(library(TxDb.Rnorvegicus.UCSC.rn7.refGene))
        seqinformation <- seqinfo(TxDb.Rnorvegicus.UCSC.rn7.refGene)
	# major chr used to subset seq info for chromosomes ... will skip 
        major_chr <- read.table("/projects/ps-epigen/users/rlan/ATACSeq/zemke_scripts/Rat_Samples/rn7.chrom.sizes")
        txs <- transcripts(TxDb.Rnorvegicus.UCSC.rn7.refGene)
}


# don't know what this is for... only use for scanBam 
possibleTag <- list("integer"=c("AM", "AS", "CM", "CP", "FI", "H0", "H1", "H2", 
                                "HI", "IH", "MQ", "NH", "NM", "OP", "PQ", "SM",
                                "TC", "UQ"), 
                 "character"=c("BC", "BQ", "BZ", "CB", "CC", "CO", "CQ", "CR",
                               "CS", "CT", "CY", "E2", "FS", "LB", "MC", "MD",
                               "MI", "OA", "OC", "OQ", "OX", "PG", "PT", "PU",
                               "Q2", "QT", "QX", "R2", "RG", "RX", "SA", "TS",
                               "U2"))

# what does this do?
bamTop100 <- scanBam(BamFile(bam, yieldSize = 100), param = ScanBamParam(tag=unlist(possibleTag)))[[1]]$tag
tags <- names(bamTop100)[lengths(bamTop100)>0]

gal <- readBamFile(bam, asMates=TRUE, bigFile=TRUE)
gal1 <- shiftGAlignmentsList(gal, outbam = paste0(bam, ".shifted.bam"))
tsse <- TSSEscore(gal1, txs)
#write.table(max(tsse$TSSEscore), file = paste0(bam, ".tsse.score.txt"), quote = F)
#system(paste0("rm ", bam, ".shifted.bam"))

fileConn <- file(paste0(bam, ".tsse.score.txt"))
writeLines(
        c(paste0("QC summary for: ", name),
          paste0("tsse: ", max(tsse$TSSEscore))
	  ), fileConn)
close(fileConn)

