# BiocManager::install("org.Hs.eg.db")
# Generate the main-text enrichment analysis
rm(list = ls())
library(clusterProfiler)
library(dplyr)
library(ggplot2)

# Set the working directory
setwd("/Volumes/Haleo_DK01/WD_Takakia2")

# Step 1: Read differential-abundance results and select significant KOs
res_df <- read.csv("wilcox_results_0612.csv", stringsAsFactors = FALSE)
diff_kos <- res_df %>%
  filter(padj < 0.05 & abs(log2FoldChange) > 0) %>%
  pull(KO_ID)  # All differentially abundant KOs

# Separate Soil-enriched (log2FC > 0) and Plant-enriched (log2FC < 0) KOs
soil_up_kos <- res_df %>%
  filter(padj < 0.05 & log2FoldChange > 0) %>%
  pull(KO_ID)

plant_up_kos <- res_df %>%
  filter(padj < 0.05 & log2FoldChange < 0) %>%
  pull(KO_ID)

# Use all KO IDs as the background
background_kos <- res_df$KO_ID

# Prepare a ranked KO vector for GSEA, with KO IDs as names and log2FC as values
gene_list <- res_df$log2FoldChange
names(gene_list) <- res_df$KO_ID
gene_list <- sort(gene_list, decreasing = TRUE)  # Sort in descending order


# Step 3: Compare KEGG enrichment for Soil- and Plant-enriched KOs
gene_lists <- list(Soil_Up = soil_up_kos, Plant_Up = plant_up_kos)

comp_kegg <- compareCluster(
  geneCluster = gene_lists,
  fun = "enrichKEGG",
  universe = background_kos,
  organism = "ko",
  pAdjustMethod = "none",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.2,
  minGSSize = 3
)

# write.csv(comp_kegg@compareClusterResult, "kegg_compare_enrichment_0_fdr_0612.csv", row.names = FALSE)

p2 <- dotplot(comp_kegg, by = "RichFactor", color = "qvalue", showCategory = 20) + 
  ggtitle("Comparative KEGG Pathway Enrichment") +
  theme_minimal() +
  theme(panel.border = element_rect(color = "black", fill = NA, size = 1))

p2
ggsave("kegg_plant_up_wilcox_0_fdr_0612.pdf", p2, width = 4, height = 5)


# Step 5: Run KEGG and GO GSEA
# GSEA for KEGG
gse_kegg <- gseKEGG(
  geneList = gene_list,  # Ranked log2FC values
  organism = "ko",
  keyType = "kegg",
  pAdjustMethod = "none",
  minGSSize = 3,
  pvalueCutoff = 0.05
  
)

write.csv(gse_kegg@result, "gse_kegg_0612.csv", row.names = FALSE)

p5 <- gseaplot(gse_kegg, geneSetID = 1)  # Example: the top-ranked pathway
p5
# ggsave("gse_kegg_plot.pdf", p5, width = 8, height = 6)