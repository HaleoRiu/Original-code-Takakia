# BiocManager::install("org.Hs.eg.db")
# 这个跑出来可以用在正文
rm(list = ls())
library(clusterProfiler)
library(dplyr)
library(ggplot2)

# 设置工作目录
setwd("/Volumes/Haleo_DK01/WD_Takakia2")

# 步骤1: 读取DESeq2结果CSV，筛选差异KO（padj < 0.05, |log2FC| >1）
res_df <- read.csv("wilcox_results_0612.csv", stringsAsFactors = FALSE)
diff_kos <- res_df %>%
  filter(padj < 0.05 & abs(log2FoldChange) > 0) %>%
  pull(KO_ID)  # 所有差异KO

# 分开上调（Soil高，log2FC > 0）和下调（Plant高，log2FC < 0）
soil_up_kos <- res_df %>%
  filter(padj < 0.05 & log2FoldChange > 0) %>%
  pull(KO_ID)

plant_up_kos <- res_df %>%
  filter(padj < 0.05 & log2FoldChange < 0) %>%
  pull(KO_ID)

# 背景KO：所有KO_ID
background_kos <- res_df$KO_ID

# 为GSEA准备排序基因列表（named vector: KO = log2FC）
gene_list <- res_df$log2FoldChange
names(gene_list) <- res_df$KO_ID
gene_list <- sort(gene_list, decreasing = TRUE)  # 降序排序


# 步骤3: 比较KEGG富集分析（Soil上调 vs Plant上调）
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


# 步骤5: 新增GSEA（KEGG和GO）
# GSEA for KEGG
gse_kegg <- gseKEGG(
  geneList = gene_list,  # 排序的log2FC
  organism = "ko",
  keyType = "kegg",
  pAdjustMethod = "none",
  minGSSize = 3,
  pvalueCutoff = 0.05
  
)

write.csv(gse_kegg@result, "gse_kegg_0612.csv", row.names = FALSE)

p5 <- gseaplot(gse_kegg, geneSetID = 1)  # 示例top1通路
p5
# ggsave("gse_kegg_plot.pdf", p5, width = 8, height = 6)