# Load necessary libraries
library(dplyr)
library(ggplot2)

# Step 1: Read in your long-format data
df <- read.csv("/home/vaibhavi/Downloads/spinal_cord/logfc_expression_change_categories.csv")

# Step 2: Simulate pseudo p-values (based on absolute logFC + noise)
set.seed(123)  # for reproducibility
df$pvalue <- 10^-(abs(df$logFC) + runif(nrow(df), 0, 1))

# Step 3: Define significance thresholds
df$significance <- ifelse(abs(df$logFC) >= 1, "Significant", "Not Significant")

# Step 4: Create volcano plot
ggplot(df, aes(x = logFC, y = -log10(pvalue), color = significance)) +
  geom_point(alpha = 0.7, size = 2) +
  scale_color_manual(values = c("Not Significant" = "grey", "Significant" = "red")) +
  theme_minimal() +
  labs(
    title = "Pseudo Volcano Plot",
    x = "Log2 Fold Change",
    y = "-log10(P-value)"
  ) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "black") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black")

# Step 5: Save the dataframe with p-values and significance
write.csv(df, "/home/vaibhavi/Downloads/spinal_cord/logfc_with_pvalues.csv", row.names = FALSE)

cat("✅ Output saved to: /home/vaibhavi/Downloads/spinal_cord/logfc_with_pvalues.csv\n")

