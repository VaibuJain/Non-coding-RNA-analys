# Install required packages if not already installed
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install(c("clusterProfiler", "org.Hs.eg.db", "AnnotationDbi", "enrichplot", "ggplot2", "dplyr"))

# Load libraries
library(clusterProfiler)
library(org.Hs.eg.db)
library(AnnotationDbi)
library(enrichplot)
library(ggplot2)
library(dplyr)

# Step 1: Load the dataset
data <- read.csv("/home/vaibhavi/Downloads/spinal_cord/logfc_with_gene_ids.csv")

# Step 2: Clean gene IDs (remove version numbers)
data$mRNA_gene_id_clean <- sub("\\..*", "", data$mRNA_gene_id)
unique_genes <- unique(data$mRNA_gene_id_clean)

# Step 3: Convert Ensembl IDs to Entrez IDs
entrez_ids <- bitr(unique_genes,
                   fromType = "ENSEMBL",
                   toType = "ENTREZID",
                   OrgDb = org.Hs.eg.db)

# Step 4: Prepare gene list for enrichment
gene_list <- na.omit(entrez_ids$ENTREZID)

# Step 5: GO Enrichment (ALL: BP, MF, CC)
go_enrich <- enrichGO(gene          = gene_list,
                      OrgDb         = org.Hs.eg.db,
                      ont           = "ALL",
                      pAdjustMethod = "BH",
                      pvalueCutoff  = 0.05,
                      readable      = TRUE)

# Step 6: KEGG Enrichment
kegg_enrich <- enrichKEGG(gene         = gene_list,
                          organism     = 'hsa',
                          pvalueCutoff = 0.05)

# Step 7: Save Results
write.csv(as.data.frame(go_enrich), "/home/vaibhavi/Downloads/spinal_cord/GO_Enrichment_Results.csv", row.names = FALSE)
write.csv(as.data.frame(kegg_enrich), "/home/vaibhavi/Downloads/spinal_cord/KEGG_Enrichment_Results.csv", row.names = FALSE)

# Step 8: Plot GO Enrichment - Barplot
barplot(go_enrich, showCategory = 10, title = "Top 10 GO Terms") +
  theme(axis.text.x = element_text(size = 12, angle = 45, hjust = 1),
        axis.text.y = element_text(size = 12),
        plot.title = element_text(size = 16, face = "bold"),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14))

# Step 9: Plot KEGG Enrichment - Barplot
barplot(kegg_enrich, showCategory = 10, title = "Top 10 KEGG Pathways") +
  theme(axis.text.x = element_text(size = 12, angle = 45, hjust = 1),
        axis.text.y = element_text(size = 12),
        plot.title = element_text(size = 16, face = "bold"),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14))

# Step 10: GO Dotplot
dotplot(go_enrich, showCategory = 10, title = "GO Dotplot") +
  theme(axis.text.x = element_text(size = 12, angle = 45, hjust = 1),
        axis.text.y = element_text(size = 12),
        plot.title = element_text(size = 16, face = "bold"),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14))

# Step 11: KEGG Dotplot
dotplot(kegg_enrich, showCategory = 10, title = "KEGG Dotplot") +
  theme(axis.text.x = element_text(size = 12, angle = 45, hjust = 1),
        axis.text.y = element_text(size = 12),
        plot.title = element_text(size = 16, face = "bold"),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14))

# Optional: Enrichment Map
emapplot(go_enrich)
