
### Old enrichr loop
runEnrichR <- function(deseqCSV_out, clusterNum, condPair, adj_pval, log2FC_threshold, output_dir){
    # Create enichR output dir
    dir.create(paste0(output_dir,"Enrichr_out"), showWarnings = FALSE)

    ### create UP- and DOWN-regulated csv files for EnrichR
    deseq_results <- read.csv(deseqCSV_out)
    # set adjusted pval cutoff
    deseq_results <- deseq_results[deseq_results$padj < adj_pval, ]

    # use deseq results and separate by UP and DOWN-regulated for enrichr
    dwn <- deseq_results[deseq_results$log2FoldChange < -log2FC_threshold, ]
    up <- deseq_results[deseq_results$log2FoldChange > log2FC_threshold, ]
    combined <- deseq_results[ which( deseq_results$log2FoldChange < -log2FC_threshold | deseq_results$log2FoldChange > log2FC_threshold ), ]

    # Save enrichr input in separate dir
    dir.create(paste0(output_dir,"DESeq_SIG/"), showWarnings = FALSE)
    up_fp <- paste(output_dir, "DESeq_SIG/", "DESeq_sigUP_c", clusterNum, "_",condPair,".csv", sep = "")
    dwn_fp <- paste(output_dir, "DESeq_SIG/", "DESeq_sigDWN_c", clusterNum, "_",condPair, ".csv", sep = "")
    write.csv(up, up_fp, row.names=FALSE, quote=FALSE)
    write.csv(dwn, dwn_fp, row.names=FALSE, quote=FALSE)

    ############ RUNNING ENRICHR ############
    E.libs <- c("MSigDB_Hallmark_2020", "GO_Biological_Process_2018", "GO_Biological_Process_2021", "KEGG_2021_Human", "GO_Molecular_Function_2021",
                "Human_Gene_Atlas", "Genes_Associated_with_NIH_Grants")

# iterate through enrichr gene lists. Run Enrichr on DESeq out and write files
    for (i in seq(1,6)){
        dir.create(paste0(output_dir, "Enrichr_out/", E.libs[i]))

        # DOWNREG
        dwn_df <- suppressMessages(data.frame(enrichr(dwn[,1], E.libs[i])))
        if (nrow(dwn_df) > 0 ){
            # Correct for commas in GO terms
            dc1.R <- str_replace_all(string = dwn_df[,1], pattern = ",", replacement = "/")
            dwn_df[,1] <- dc1.R
        }
        dwnfile <- paste(output_dir, "Enrichr_out/", E.libs[i], "/DWNreg_", clusterNum, "_", condPair, ".csv", sep = "")
        write.csv(dwn_df, dwnfile, row.names=FALSE, quote=FALSE)
        if (nrow(dwn_df > 1)){
            dwnfile_html <- paste(output_dir, "Enrichr_out/", E.libs[i], "/DWNreg_", clusterNum, "_", condPair, ".html", sep = "")
            write_tableHTML(tableHTML(dwn_df), file = dwnfile_html)
        }

        # UPREG
        up_df <- suppressMessages(data.frame(enrichr(up[,1], E.libs[i])))
        if (nrow(up_df) > 0 ){
           uc1.R <- str_replace_all(string = up_df[,1], pattern = ",", replacement = "/")
           up_df[,1] <- uc1.R
        }
        upfile <- paste(output_dir, "Enrichr_out/", E.libs[i], "/UPreg_", clusterNum, "_", condPair, ".csv", sep = "")
        write.csv(up_df, upfile, row.names=FALSE, quote=FALSE)
        if (nrow(up_df > 1)){
           upfile_html <- paste(output_dir, "Enrichr_out/", E.libs[i], "/UPreg_", clusterNum, "_", condPair, ".html", sep = "")
           write_tableHTML(tableHTML(up_df), file = upfile_html)
        }

        # COMBINED
        comb_df <- suppressMessages(data.frame(enrichr(combined[,1], E.libs[i])))
        if (nrow(comb_df) > 0 ){
           cc1.R <- str_replace_all(string = comb_df[,1], pattern = ",", replacement = "/")
           comb_df[,1] <- cc1.R
        }
        cfile <- paste(output_dir, "Enrichr_out/", E.libs[i], "/COMB_", clusterNum, "_", condPair, ".csv", sep = "")
        write.csv(comb_df, cfile, row.names=FALSE, quote=FALSE)
        if (nrow(comb_df > 1)){
        cfile_html <- paste(output_dir, "Enrichr_out/", E.libs[i], "/COMB_", clusterNum, "_", condPair, ".html", sep = "")
        write_tableHTML(tableHTML(comb_df), file = cfile_html)
      }
   }
}



