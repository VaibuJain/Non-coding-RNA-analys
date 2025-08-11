# Load required libraries
library(dplyr)

# Step 1: Load expression matrices
lncrna_counts <- read.delim("/home/vaibhavi/Downloads/spinal_cord/lncrna_count..csv", header = TRUE, sep = "\t", stringsAsFactors = FALSE)
rownames(lncrna_counts) <- lncrna_counts$lncRNA_ID
lncrna_counts <- lncrna_counts[, -1]  # Remove lncRNA_ID column

mrna_counts <- read.csv("/home/vaibhavi/Downloads/spinal_cord/filtered_mrna_counts.csv", row.names = 1, check.names = FALSE)

# Step 2: Load overlap file
overlaps <- read.csv("/home/vaibhavi/Downloads/spinal_cord/lncRNA_mRNA_overlaps_full.csv", stringsAsFactors = FALSE)

# Step 3: Rename columns if needed
colnames(overlaps)[colnames(overlaps) == "lncRNA_id"] <- "lncRNA_ID"
colnames(overlaps)[colnames(overlaps) == "mRNA_gene_id"] <- "mRNA_ID"

# Step 4: Drop rows with missing IDs
overlaps <- overlaps %>% filter(!is.na(lncRNA_ID), !is.na(mRNA_ID))

# Step 5: Get unique lncRNA–mRNA pairs
lnc_mrna_pairs <- overlaps[, c("lncRNA_ID", "mRNA_ID")] %>% distinct()

# Step 6: Initialize results dataframe
results <- data.frame(
  lncRNA = character(),
  mRNA = character(),
  correlation = numeric(),
  pvalue = numeric(),
  stringsAsFactors = FALSE
)

# Step 7: Correlation analysis
for (i in 1:nrow(lnc_mrna_pairs)) {
  lnc_id <- lnc_mrna_pairs$lncRNA_ID[i]
  mrna_id <- lnc_mrna_pairs$mRNA_ID[i]
  
  if (lnc_id %in% rownames(lncrna_counts) && mrna_id %in% rownames(mrna_counts)) {
    lnc_expr <- as.numeric(lncrna_counts[lnc_id, ])
    mrna_expr <- as.numeric(mrna_counts[mrna_id, ])
    
    # Check if both have non-zero variance
    if (sd(lnc_expr) != 0 && sd(mrna_expr) != 0) {
      cor_test <- cor.test(lnc_expr, mrna_expr, method = "pearson")
      results <- rbind(results, data.frame(
        lncRNA = lnc_id,
        mRNA = mrna_id,
        correlation = cor_test$estimate,
        pvalue = cor_test$p.value
      ))
    }
  }
}

# Step 8: Save correlation results
write.csv(results, "/home/vaibhavi/Downloads/spinal_cord/lncRNA_mRNA_correlation_results.csv", row.names = FALSE)

# Step 9: Print summary
cat("\nCorrelation analysis completed.")
cat("\nNumber of correlated pairs:", nrow(results), "\n")

