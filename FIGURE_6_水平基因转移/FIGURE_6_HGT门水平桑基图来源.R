# 设置工作目录
setwd("/Volumes/Haleo_DK01/WD_Takakia_F/HGT")
rm(list = ls())
# 加载必要的包
library(tidyverse)
library(ggalluvial)
library(ggraph)
library(igraph)

# ==========================================
# 第一部分：数据深度解析与整合
# ==========================================

# 1. 读取所有原始文件
list_data <- read.csv("list.csv", header = FALSE, col.names = c("filename"))
group_data <- read.table("group.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)
tax_data <- read.table("taxonomy.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)
donor_data <- read.table("donor.tsv", header = TRUE, sep = "\t", quote = "", comment.char = "")

# 2. 从 list.csv 的 SVG 文件名中精准提取供体与受体 MAG 名称
# 文件名格式例如：Bin_89_01043___Bin_121_01510.SVG
hgt_events <- list_data %>%
  mutate(
    # 去除 .SVG 后缀
    clean_name = str_remove(filename, "\\.SVG$|\\.svg$"),
    # 以 "___" 分割供体和受体基因
    donor_gene = str_split_i(clean_name, "___", 1),
    recipient_gene = str_split_i(clean_name, "___", 2),
    # 去除基因编号后缀 (例如 _01043)，保留 MAG 名称
    Donor_MAG = str_replace(donor_gene, "_[0-9]+$", ""),
    Recipient_MAG = str_replace(recipient_gene, "_[0-9]+$", "")
  )

# 3. 清洗并简化分类学门类 (Phylum) 名称
tax_clean <- tax_data %>%
  mutate(Phylum = str_remove(Phylum, "^p__")) %>%
  select(Name, Phylum)

# 4. 构建包含完整信息的 HGT 边缘数据集
hgt_full <- hgt_events %>%
  # 整合供体分组与分类
  left_join(group_data, by = c("Donor_MAG" = "name")) %>%
  rename(Donor_Group = group) %>%
  left_join(tax_clean, by = c("Donor_MAG" = "Name")) %>%
  rename(Donor_Phylum = Phylum) %>%
  # 整合受体分组与分类
  left_join(group_data, by = c("Recipient_MAG" = "name")) %>%
  rename(Recipient_Group = group) %>%
  left_join(tax_clean, by = c("Recipient_MAG" = "Name")) %>%
  rename(Recipient_Phylum = Phylum) %>%
  # 整合 COG 功能注释 (提取 COG_category 的首字母)
  left_join(donor_data %>% select(X.query, COG_category), by = c("donor_gene" = "X.query")) %>%
  mutate(COG_Main = str_sub(COG_category, 1, 1)) %>%
  # 处理缺失值
  replace_na(list(COG_Main = "Unknown", Donor_Phylum = "Unclassified", Recipient_Phylum = "Unclassified"))

# ==========================================
# 第二部分：多维信息冲积图绘制 (升级版桑基图)
# ==========================================

# 汇总各类流动路径的数量
alluvial_data <- hgt_full %>%
  count(Donor_Group, Donor_Phylum, Recipient_Phylum, Recipient_Group)

p1_alluvial <- ggplot(data = alluvial_data,
                      aes(axis1 = Donor_Group, axis2 = Donor_Phylum, 
                          axis3 = Recipient_Phylum, axis4 = Recipient_Group,
                          y = n)) +
  scale_x_discrete(limits = c("Donor Group", "Donor Phylum", "Recipient Phylum", "Recipient Group"), expand = c(.1, .05)) +
  geom_alluvium(aes(fill = Donor_Phylum), alpha = 0.7, curve_type = "cubic") +
  geom_stratum(alpha = 0.9, color = "white") +
  geom_text(stat = "stratum", aes(label = after_stat(stratum)), size = 3) +
  theme_minimal() +
  labs(title = "Taxonomic and Spatial Flow of Horizontal Gene Transfer",
       y = "Number of HGT Events", fill = "Donor Phylum") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text.x = element_text(size = 12, face = "bold"))

print(p1_alluvial)
#ggsave("Enhanced_Sankey_Taxonomy.pdf", p1_alluvial, width = 9, height = 6)

# ==========================================
# 第三部分：群落 HGT 物种互作网络图绘制 (修改版)
# ==========================================

# 1. 整理节点的 Phylum 信息，聚焦核心门类，其余归为 "Others"
# 构建网络边数据与节点元数据
edges <- hgt_full %>% count(Donor_MAG, Recipient_MAG, name = "weight")
nodes <- data.frame(MAG = unique(c(edges$Donor_MAG, edges$Recipient_MAG))) %>%
  left_join(group_data, by = c("MAG" = "name")) %>%
  left_join(tax_clean, by = c("MAG" = "Name"))

# 构建 igraph 对象
net <- graph_from_data_frame(d = edges, vertices = nodes, directed = TRUE)
# 计算节点的入度与出度作为节点大小
V(net)$degree <- degree(net, mode = "all")

p2_network <- ggraph(net, layout = "fr") + 
  geom_edge_fan(aes(edge_width = weight), alpha = 0.4, color = "gray50", 
                arrow = arrow(length = unit(2, 'mm')), end_cap = circle(3, 'mm')) +
  geom_node_point(aes(color = group, size = degree), alpha = 0.8) +
  scale_color_manual(values = c("E" = "#E69F00", "R" = "#56B4E9")) +
  theme_void() +
  labs(title = "HGT Interaction Network across Microenvironments",
       color = "Compartment", size = "HGT Events (Degree)") +
  theme(legend.position = "right")

print(p2_network)
# ggsave("HGT_Network_Hubs.pdf", p2_network, width = 10, height = 8)

target_phyla <- c("Acidobacteriota", "Actinomycetota", "Bacteroidota", 
                  "Eremiobacterota", "Gemmatimonadota", "Myxococcota", 
                  "Pseudomonadota")

nodes <- nodes %>%
  mutate(Phylum_Group = ifelse(Phylum %in% target_phyla, Phylum, "Others")) %>%
  # 确保 Others 在图例中排在最后
  mutate(Phylum_Group = factor(Phylum_Group, levels = c(target_phyla, "Others")))

# 2. 构建 igraph 对象并计算节点度（Degree）
net <- graph_from_data_frame(d = edges, vertices = nodes, directed = TRUE)
V(net)$degree <- degree(net, mode = "all")

# 3. 设置自定义的颜色和形状映射字典
phylum_colors <- c(
  "Acidobacteriota" = "#f8766d",
  "Actinomycetota"  = "#c49a00",
  "Bacteroidota"    = "#53b400",
  "Eremiobacterota" = "#00c094",
  "Gemmatimonadota" = "#00b6eb",
  "Myxococcota"     = "#a58aff",
  "Pseudomonadota"  = "#fb61d7",
  "Others"          = "#d4d4d4"
)

# 16 代表实心圆，17 代表实心三角。你可以根据偏好将它们分配给 E 和 R
group_shapes <- c("E" = 16, "R" = 17) 

# 4. 绘制 Fruchterman-Reingold 布局的网络图
p2_network_modified <- ggraph(net, layout = "fr") + 
  # 绘制连线：适当降低透明度，避免遮挡节点
  geom_edge_fan(aes(edge_width = weight), alpha = 0.25, color = "gray60", 
                arrow = arrow(length = unit(2.5, 'mm')), end_cap = circle(3, 'mm')) +
  # 绘制节点：映射颜色(门类)、形状(分组)和大小(度)
  geom_node_point(aes(color = Phylum_Group, shape = group, size = degree), alpha = 0.9) +
  # 应用自定义的映射
  scale_color_manual(values = phylum_colors) +
  scale_shape_manual(values = group_shapes) +
  scale_size_continuous(range = c(3, 11)) + # 拉开节点大小的视觉差距，凸显 Hubs
  theme_void() +
  labs(color = "Phylum", shape = "Compartment", size = "HGT Degree") +
  # 优化图例排版，使其适合放入论文正文
  theme(
    legend.position = "right",
    legend.title = element_text(face = "bold", size = 12),
    legend.text = element_text(size = 10),
    legend.key = element_blank()
  ) +
  # 让点大小的图例保持统一的颜色，避免受 Phylum 颜色干扰
  guides(
    color = guide_legend(override.aes = list(size = 5)),
    shape = guide_legend(override.aes = list(size = 5)),
    size = guide_legend(override.aes = list(color = "gray30"))
  )

print(p2_network_modified)

# 建议导出为矢量图以备后续排版
ggsave("HGT_Network_Modified_FR.pdf", p2_network_modified, width = 11, height = 8)
