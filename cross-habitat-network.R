rm(list = ls())

library(Hmisc)
library(igraph)
library(dplyr)
library(tidyr)

setwd("/path/to/your/files")

asv_data <- read.delim("asv_bacteria.txt", header = TRUE, sep = "\t", row.names = 1, check.names = FALSE)
taxonomy <- read.delim("taxonomy_bacteria.txt", header = TRUE, sep = "\t", row.names = 1, check.names = FALSE)

total_reads <- rowSums(asv_data)
asv_data <- asv_data[total_reads >= 10, ]
cat("ASV count after removing total reads < 10:", nrow(asv_data), "\n")

sample_occurrence <- rowSums(asv_data > 0)
asv_data <- asv_data[sample_occurrence >= 6, ]
cat("ASV count after requiring occurrence in at least 6 samples:", nrow(asv_data), "\n")

taxonomy <- taxonomy[rownames(asv_data), ]

soil_samples <- grep("^S", colnames(asv_data), value = TRUE)
plant_samples <- grep("^P", colnames(asv_data), value = TRUE)
asv_soil <- asv_data[, soil_samples]
asv_plant <- asv_data[, plant_samples]

prev_soil <- apply(asv_soil > 0, 1, mean) > 0.1
prev_plant <- apply(asv_plant > 0, 1, mean) > 0.1
shared_asvs <- rownames(asv_data)[prev_soil & prev_plant]

combined_asv <- t(asv_data[shared_asvs, ])

combined_asv <- combined_asv / rowSums(combined_asv)
combined_asv[is.na(combined_asv)] <- 0

corr_result <- rcorr(combined_asv, type = "spearman") 
corr_matrix <- corr_result$r
p_matrix <- corr_result$P

p_matrix_adj <- p.adjust(p_matrix, method = "fdr")

adj_matrix <- corr_matrix
adj_matrix[p_matrix_adj >= 0.05 | abs(corr_matrix) <= 0.6] <- 0
diag(adj_matrix) <- 0

edges <- as.data.frame(which(adj_matrix != 0, arr.ind = TRUE))
edges$weight <- adj_matrix[as.matrix(edges[, 1:2])]
edges <- edges %>% filter(row < col)
colnames(edges) <- c("from", "to", "weight")
edges$from <- rownames(corr_matrix)[edges$from]
edges$to <- rownames(corr_matrix)[edges$to]

net_graph <- graph_from_data_frame(edges, directed = FALSE, vertices = shared_asvs)

V(net_graph)$Genus <- taxonomy[V(net_graph)$name, "Genus"]
V(net_graph)$label <- ifelse(is.na(V(net_graph)$Genus) | V(net_graph)$Genus == "", V(net_graph)$name, V(net_graph)$Genus)

write_graph(net_graph, "cross_niche_network.graphml", format = "graphml")