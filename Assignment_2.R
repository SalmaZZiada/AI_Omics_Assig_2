classify_gene <- function(logFC, padj) {
  if (is.na(padj)) padj <- 1           
  if (logFC > 1 & padj < 0.05) {
    return("Upregulated")
  } else if (logFC < -1 & padj < 0.05) {
    return("Downregulated")
  } else {
    return("Not_Significant")
  }
}


input_dir  <- "Raw_Data"
output_dir <- "Results"
if (!dir.exists(output_dir)) dir.create(output_dir)

files_to_process <- c("DEGs_Data_1.csv", "DEGs_Data_2.csv")



results_list <- list()
for (file in files_to_process) {
  cat("\nProcessing:", file, "\n")
  input_file <- file.path(input_dir, file)
  data <- read.csv(input_file, header = TRUE)
  
  names(data) <- tolower(names(data))      
  
  data$padj[is.na(data$padj)] <- 1
  
  data$status <- mapply(classify_gene, data$logfc, data$padj)
  
  output_file <- file.path(output_dir, paste0("Processed_", file))
  write.csv(data, output_file, row.names = FALSE)
  cat("Saved ->", output_file, "\n")
  
  cat("Summary (status counts):\n")
  print(table(data$status))
    results_list[[file]] <- data
}



results_1 <- results_list[["DEGs_Data_1.csv"]]
results_2 <- results_list[["DEGs_Data_2.csv"]]

cat("\nAll done. Check the 'Results' folder and summaries above.\n")