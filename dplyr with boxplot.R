# Load required libraries
library(dplyr)
library(tidyr)
library(ggplot2)
library(pheatmap)

# Step 1: Load data
logfc_df <- read.csv("/home/vaibhavi/Downloads/spinal_cord/count_matrix.csv", header = TRUE, check.names = FALSE)

# Step 2: Rename columns
colnames(logfc_df) <- c("transcript_id", "LogFC_978_980", "LogFC_978_981", "LogFC_979_980", "LogFC_979_981")

# Step 3: Filter by absolute logFC ≥ 1 in any comparison
filtered_df <- logfc_df %>%
  filter(
    abs(LogFC_978_980) >= 1 |
      abs(LogFC_978_981) >= 1 |
      abs(LogFC_979_980) >= 1 |
      abs(LogFC_979_981) >= 1
  )

# Step 4: Export filtered data
write.csv(filtered_df, "/home/vaibhavi/Downloads/spinal_cord/filtered_logfc_transcripts.csv", row.names = FALSE)

# Step 5: Heatmap of LogFC
plot_matrix <- as.matrix(filtered_df[, -1])
rownames(plot_matrix) <- filtered_df$transcript_id
pheatmap(plot_matrix, cluster_rows = TRUE, cluster_cols = TRUE, main = "LogFC Heatmap")

# -----------------------------
# Step 6: Prepare for boxplot
# -----------------------------

# Convert to long format for ggplot2
long_df <- filtered_df %>%
  pivot_longer(cols = starts_with("LogFC"), names_to = "comparison", values_to = "logFC")

# Categorize by logFC magnitude
long_df <- long_df %>%
  mutate(
    expression_change = case_when(
      abs(logFC) < 1 ~ "Low",
      abs(logFC) >= 1 & abs(logFC) < 2 ~ "Moderate",
      abs(logFC) >= 2 ~ "High"
    )
  )

# -----------------------------
# Step 7: Boxplot
# -----------------------------
ggplot(long_df, aes(x = expression_change, y = logFC, fill = expression_change)) +
  geom_boxplot() +
  theme_minimal() +
  labs(title = "Boxplot of LogFC by Expression Change Category",
       x = "Expression Change Category", y = "Log Fold Change") +
  scale_fill_manual(values = c("Low" = "lightblue", "Moderate" = "orange", "High" = "red"))

# Step 6.1: Save the categorized long-format data
write.csv(long_df, "/home/vaibhavi/Downloads/spinal_cord/logfc_expression_change_categories.csv", row.names = FALSE)


