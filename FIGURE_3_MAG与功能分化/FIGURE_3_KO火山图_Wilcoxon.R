# 清空环境
rm(list = ls())

# 设置工作目录
setwd("/Volumes/Haleo_DK01/WD_Takakia2")

# 加载必要的包
library(tidyverse)
library(ggplot2)
library(ggrepel)

# ------------------------------------------------------------------
# 1. 数据读取与预处理
# ------------------------------------------------------------------

# 读取 TPM 数据（替换counts_df）
tpm_df <- read.table("tpm_renormalized.txt", header = T, row.names = 1, sep = "\t", check.names = F)

# 根据列名自动生成分组信息 (P开头=Plant, S开头=Soil)
sample_names <- colnames(tpm_df)
condition <- ifelse(grepl("^P", sample_names), "Plant", "Soil")
colData <- data.frame(row.names = sample_names, condition = factor(condition, levels = c("Plant", "Soil")))

# 检查分组
print(table(colData$condition))

# log变换TPM（+1避免0）
# tpm_log <- log2(tpm_df + 1)
tpm_log <- tpm_df
# ------------------------------------------------------------------
# 2. 计算差异（用t-test或Wilcoxon替换DESeq2）
# ------------------------------------------------------------------

# 函数：计算log2FC和p值
calc_diff <- function(row, group1_samples, group2_samples, method = "wilcox") {
  group1 <- row[group1_samples]
  group2 <- row[group2_samples]
  log2fc <- log2((mean(group1) + 1e-6) / (mean(group2) + 1e-6))  # +小数避免除0，Soil/Plant
  if (method == "t") {
    test_res <- t.test(group1, group2)
  } else if (method == "wilcox") {
    test_res <- wilcox.test(group1, group2)
  }
  pval <- test_res$p.value
  return(c(log2FC = log2fc, pvalue = pval))
}

# 应用到每个KO（用tpm_log）
plant_samples <- rownames(colData)[colData$condition == "Plant"]
soil_samples <- rownames(colData)[colData$condition == "Soil"]

res_list <- apply(tpm_log, 1, calc_diff, group1_samples = soil_samples, group2_samples = plant_samples, method = "wilcox")  # 改"t"为t-test
res_df <- as.data.frame(t(res_list))  # 转置
res_df$KO_ID <- rownames(res_df)
res_df$padj <- p.adjust(res_df$pvalue, method = "BH")  # Benjamini-Hochberg校正
res_df$baseMean <- rowMeans(tpm_df)  # 平均TPM（可选过滤）

# 保存结果
write.csv(res_df, "wilcox_results_0612.csv", row.names = FALSE)  # 或t_test_results.csv

# ------------------------------------------------------------------
# 3. 准备绘图数据（其余不变）
# ------------------------------------------------------------------

# 定义阈值
fc_cutoff <- 0      # log2FC 阈值
p_cutoff <- 0.05    # padj 阈值

# 添加分类标签：Up (Soil高), Down (Plant高), NS
res_df <- res_df %>%
  mutate(group = case_when(
    padj < p_cutoff & log2FC > fc_cutoff ~ "Up",    # Soil显著高
    padj < p_cutoff & log2FC < -fc_cutoff ~ "Down", # Plant显著高
    TRUE ~ "NS"
  ))

res_df$group <- factor(res_df$group, levels = c("Up", "Down", "NS"))

# 筛选 Top 10 Up/Down（按padj排序）
top_up <- res_df %>%
  filter(group == "Up") %>%
  arrange(padj) %>%
  head(10)

top_down <- res_df %>%
  filter(group == "Down") %>%
  arrange(padj) %>%
  head(10)

top_genes <- bind_rows(top_up, top_down)

# ------------------------------------------------------------------
# 3. 定义抗逆类别KO列表并标注
# ------------------------------------------------------------------

# DNA修复
dna_repair_kos <- c("K03648", "K10771", "K03660", "K10847", "K10844", "K10843", "K08735", "K04482", "K10866", "K10865", "K10867", "K10873", "K10901", "K10605", "K08775", "K10875", "K10879", "K10880", "K10843", "K10849", "K10846", "K10848", "K10747", "K02330")

# 氧化应激
oxidative_stress_kos <- c("K04565", "K04564", "K03781", "K23856", "K24158", "K00799", "K03671", "K00384", "K00131", "K00434", "K00430", "K08232", "K00128")

# 类胡萝卜素
carotenoids_kos <- c("K02291", "K02293", "K00514", "K09835", "K06443", "K22502", "K02294", "K09838", "K09839", "K09840", "K10027", "K01669", "K06976")

# 氨基酸 (结合半胱氨酸和谷氨酸)
amino_acids_kos <- c("K00640", "K01738", "K12339", "K10150", "K01697", "K00261", "K00262", "K01425", "K01919", "K11205", "K14454", "K14455", "K01915")

# 萜类及前体
terpenoids_kos <- c("K00626", "K01641", "K00021", "K00869", "K00869", "K01597", "K01823", "K01662", "K00099", "K00991", "K00919", "K01770", "K03526", "K13789", "K00787")

# 黄酮类
flavonoids_kos <- c("K10775", "K00660", "K01859", "K13406", "K23135", "K13083", "K05278", "K00091", "K13082", "K05277", "K13081", "K21102", "K05281", "K13064")

# 冷应激
cold_stress_kos <- c("K11574", "K11575", "K11576", "K11482", "K03704")

# 铁载体相关
iron_carrier_kos <- c("K02361", "K00216", "K01252", "K02364", "K08225", "K19611", "K25290")

# 海藻糖相关
trehalose_kos <- c("K02361", "K01087", "K01194", "K01194", "K00130", "K00286")

# VB12及ACC
vb12_acc_kos <- c("K01505", "K02303", "K05934", "K05936", "K00595", "K02224", "K02232")

# 定殖、生物膜等
colonization_kos <- c("K07173", "K10909", "K10911", "K10912", "K13061", "K18099", "K19421", "K19422", "K19423")

# 所有类别列表
categories <- list(
  DNA_repair = dna_repair_kos,
  Oxidative_stress = oxidative_stress_kos,
  Carotenoids = carotenoids_kos,
  Amino_acids = amino_acids_kos,
  Terpenoids = terpenoids_kos,
  Flavonoids = flavonoids_kos,
  Cold_stress = cold_stress_kos,
  Iron_carrier = iron_carrier_kos,
  Trehalose = trehalose_kos,
  VB12_ACC = vb12_acc_kos,
  Colonization = colonization_kos
)


res_df$category <- "Other"
for (cat_name in names(categories)) {
  res_df$category[res_df$KO_ID %in% categories[[cat_name]]] <- cat_name
}

res_df$color_group <- ifelse(res_df$category == "Other", as.character(res_df$group), res_df$category)
res_df$color_group <- factor(res_df$color_group, levels = c("Up", "Down", "NS", names(categories)))

# ------------------------------------------------------------------
# 4. 使用 ggplot2 绘图（其余不变）
# ------------------------------------------------------------------

p <- ggplot(res_df, aes(x = log2FC, y = -log10(padj), color = color_group)) +
  # 1. 先绘制所有点（底层）
  geom_point(alpha = 0.6, size = 1.5) +
  
  # 2. 再单独绘制抗逆相关点（顶层），使用subset数据
  geom_point(data = subset(res_df, category != "Other"), alpha = 0.8, size = 2) +  # alpha稍高，size稍大以突出
  
  # 3. 自定义颜色 (原有group + 类别颜色，未变)
  scale_color_manual(values = c(
    "Up" = "#01B7ED",          # Soil高 (蓝)
    "Down" = "#37B484",        # Plant高 (绿)
    "NS" = "#999999",          # 无差异 (灰)
    "DNA_repair" = "purple",   # 紫色
    "Oxidative_stress" = "red",# 红色
    "Carotenoids" = "orange",  # 橙色
    "Amino_acids" = "brown",   # 棕色 (结合Cys和Glu)
    "Terpenoids" = "darkgreen",# 深绿
    "Flavonoids" = "yellow",   # 黄色
    "Cold_stress" = "cyan",    # 青色
    "Iron_carrier" = "magenta",# 品红
    "Trehalose" = "lightblue", # 浅蓝
    "VB12_ACC" = "violet",     # 紫罗兰
    "Colonization" = "darkred" # 深红
  ),
  labels = c("Enriched in Soil", "Enriched in Plant", "NS", names(categories)),
  drop = FALSE) +  # drop=FALSE确保所有levels显示在图例中，即使无数据
  
  # 4. 添加阈值辅助线
  geom_vline(xintercept = c(-fc_cutoff, fc_cutoff), linetype = "dashed", color = "black", alpha=0.5) +
  geom_hline(yintercept = -log10(p_cutoff), linetype = "dashed", color = "black", alpha=0.5) +
  
  # 5. 添加 Top 基因标签 (使用 ggrepel 防止重叠)
  geom_text_repel(
    data = top_genes,
    aes(label = KO_ID),
    size = 3,
    color = "black",      # 标签统一黑色
    box.padding = 0.5,    # 文字框周围留白
    max.overlaps = 50,    # 允许重叠的最大数量
    segment.color = "grey50" # 连接线颜色
  ) +
  
  # 6. 设置坐标轴和标题
  labs(x = expression(log[2] ~ "Fold Change"), 
       y = expression(-log[10] ~ "Adjusted P-value"),
       title = "Soil vs Plant (KO Abundance)",
       subtitle = paste0("Up: ", sum(res_df$group == "Up"), 
                         " | Down: ", sum(res_df$group == "Down")),
       color = "Group/Category") +  # 图例标题调整
  
  # 7. 主题美化 (Classic 主题最干净，适合出版)
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5),
    legend.position = "right",
    legend.background = element_rect(fill = "transparent"),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10)
  )

# ------------------------------------------------------------------
# 5. 保存结果
# ------------------------------------------------------------------

# 显示图片
print(p)


# 保存 PDF
ggsave("Volcano_Plot_wilcox_12_0612.pdf", plot = p, width = 6, height = 5)
