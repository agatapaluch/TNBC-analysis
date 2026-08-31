library(DESeq2)
library(biomaRt)

data_dir <- "../data/rna/combined_data"
mart <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")

process_cohort <- function(cohort, output_file) {
  files <- list.files(data_dir, pattern = paste0("^", cohort, ".*\\.txt$"), full.names = TRUE)
  sample_ids <- sub("_#_.*", "", basename(files))
  
  read_count_file <- function(file, sample_name) {
    df <- read.table(file, header = FALSE, stringsAsFactors = FALSE)
    colnames(df) <- c("gene_id", sample_name)
    df <- df[!grepl("^__", df$gene_id), ]
    return(df)
  }
  
  list_of_dfs <- mapply(read_count_file, file = files, sample_name = sample_ids, SIMPLIFY = FALSE)
  merged_df <- Reduce(function(x, y) merge(x, y, by = "gene_id", all = TRUE), list_of_dfs)
  
  count_matrix <- as.matrix(merged_df[, -1])
  rownames(count_matrix) <- merged_df$gene_id
  count_matrix[is.na(count_matrix)] <- 0
  
  keep <- rowSums(count_matrix) >= 10
  count_matrix <- count_matrix[keep, ]
  
  rownames(count_matrix) <- sub("\\..*", "", rownames(count_matrix))
  count_matrix <- rowsum(count_matrix, group = rownames(count_matrix), reorder = FALSE)
  
  genes <- rownames(count_matrix)
  
  mapping <- getBM(
    attributes = c("ensembl_gene_id", "hgnc_symbol"),
    filters = "ensembl_gene_id",
    values = genes,
    mart = mart
  )
  
  mapping <- unique(mapping)
  mapping <- mapping[mapping$hgnc_symbol != "", ]
  
  mapping_counts <- table(mapping$ensembl_gene_id)
  mapping <- mapping[mapping$ensembl_gene_id %in% names(mapping_counts[mapping_counts == 1]), ]
  
  count_df <- as.data.frame(count_matrix)
  count_df$ensembl_gene_id <- rownames(count_df)
  count_df <- merge(mapping, count_df, by = "ensembl_gene_id")
  count_df <- count_df[, c("hgnc_symbol", sample_ids)]
  
  count_matrix <- as.matrix(count_df[, -1])
  rownames(count_matrix) <- count_df$hgnc_symbol
  count_matrix <- rowsum(count_matrix, group = rownames(count_matrix), reorder = FALSE)
  
  coldata <- data.frame(row.names = colnames(count_matrix))
  
  dds <- DESeqDataSetFromMatrix(
    countData = count_matrix,
    colData = coldata,
    design = ~ 1
  )
  
  vsd <- vst(dds)
  mat_vsd <- assay(vsd)
  
  mat_final <- data.frame(Name = rownames(mat_vsd), mat_vsd, check.names = FALSE)
  write.table(mat_final, file = output_file, sep = "\t", quote = FALSE, row.names = FALSE)
}

process_cohort("IMMU", "R_output/immu_expression_matrix.tsv")
process_cohort("SYNG", "R_output/synergy_expression_matrix.tsv")

