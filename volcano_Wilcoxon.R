rm(list = ls())
setwd("/path/to/your/files")


library(tidyverse)
library(ggplot2)
library(ggrepel)


tpm_df <- read.table("tpm.txt", header = T, row.names = 1, sep = "\t", check.names = F)

sample_names <- colnames(tpm_df)
condition <- ifelse(grepl("^P", sample_names), "Plant", "Soil")
colData <- data.frame(row.names = sample_names, condition = factor(condition, levels = c("Plant", "Soil")))

calc_diff <- function(row, group1_samples, group2_samples, method = "wilcox") {
  group1 <- row[group1_samples]
  group2 <- row[group2_samples]
  log2fc <- log2((mean(group1) + 1e-6) / (mean(group2) + 1e-6))
  if (method == "t") {
    test_res <- t.test(group1, group2)
  } else if (method == "wilcox") {
    test_res <- wilcox.test(group1, group2)
  }
  pval <- test_res$p.value
  return(c(log2FC = log2fc, pvalue = pval))
}

plant_samples <- rownames(colData)[colData$condition == "Plant"]
soil_samples <- rownames(colData)[colData$condition == "Soil"]

res_list <- apply(tpm_df, 1, calc_diff, group1_samples = soil_samples, group2_samples = plant_samples, method = "wilcox")  # Use "t" to select a t-test
res_df <- as.data.frame(t(res_list))
res_df$KO_ID <- rownames(res_df)
res_df$padj <- p.adjust(res_df$pvalue, method = "BH")
res_df$baseMean <- rowMeans(tpm_df)

write.csv(res_df, "wilcox_results.csv", row.names = FALSE)


fc_cutoff <- 0
p_cutoff <- 0.05

res_df <- res_df %>%
  mutate(group = case_when(
    padj < p_cutoff & log2FC > fc_cutoff ~ "Up",    # R enriched
    padj < p_cutoff & log2FC < -fc_cutoff ~ "Down", # E enriched
    TRUE ~ "NS"
  ))

res_df$group <- factor(res_df$group, levels = c("Up", "Down", "NS"))


p <- ggplot(res_df, aes(x = log2FC, y = -log10(padj), color = color_group)) +
  geom_point(alpha = 0.6, size = 1.5) +
  
  geom_point(data = subset(res_df, category != "Other"), alpha = 0.8, size = 2) +  # Increase opacity and point size for emphasis
  
  scale_color_manual(values = c(
    "Up" = "#01B7ED",
    "Down" = "#37B484",
    "NS" = "#999999",
  ),
  labels = c("Enriched in R", "Enriched in E", "NS", names(categories)),
  drop = FALSE) +  
  
  geom_vline(xintercept = c(-fc_cutoff, fc_cutoff), linetype = "dashed", color = "black", alpha=0.5) +
  geom_hline(yintercept = -log10(p_cutoff), linetype = "dashed", color = "black", alpha=0.5) +
  
  labs(x = expression(log[2] ~ "Fold Change"), 
       y = expression(-log[10] ~ "Adjusted P-value"),
       title = "Soil vs Plant (KO Abundance)",
       subtitle = paste0("Up: ", sum(res_df$group == "Up"), 
                         " | Down: ", sum(res_df$group == "Down")),
       color = "Group/Category") +
  

  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5),
    legend.position = "right",
    legend.background = element_rect(fill = "transparent"),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10)
  )

print(p)

ggsave("Volcano_Plot_wilcox.pdf", plot = p, width = 6, height = 5)
