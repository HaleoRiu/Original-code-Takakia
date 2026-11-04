# Set the working directory
setwd("/Volumes/Haleo_DK01/WD_Takakia_F/HGT")
rm(list = ls())
# Load required packages
library(tidyverse)
library(ggalluvial)
library(ggraph)
library(igraph)

# ==========================================
# Part 1: Parse and integrate the data
# ==========================================

# 1. Read all input files
list_data <- read.csv("list.csv", header = FALSE, col.names = c("filename"))
group_data <- read.table("group.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)
tax_data <- read.table("taxonomy.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)
donor_data <- read.table("donor.tsv", header = TRUE, sep = "\t", quote = "", comment.char = "")

# 2. Extract donor and recipient MAG names from SVG filenames in list.csv
# Example filename: Bin_89_01043___Bin_121_01510.SVG
hgt_events <- list_data %>%
  mutate(
    # Remove the .SVG extension
    clean_name = str_remove(filename, "\\.SVG$|\\.svg$"),
    # Split donor and recipient genes at "___"
    donor_gene = str_split_i(clean_name, "___", 1),
    recipient_gene = str_split_i(clean_name, "___", 2),
    # Remove the gene-number suffix (e.g. _01043), retaining the MAG name
    Donor_MAG = str_replace(donor_gene, "_[0-9]+$", ""),
    Recipient_MAG = str_replace(recipient_gene, "_[0-9]+$", "")
  )

# 3. Clean and simplify phylum names
tax_clean <- tax_data %>%
  mutate(Phylum = str_remove(Phylum, "^p__")) %>%
  select(Name, Phylum)

# 4. Build the annotated HGT edge dataset
hgt_full <- hgt_events %>%
  # Join donor groups and taxonomy
  left_join(group_data, by = c("Donor_MAG" = "name")) %>%
  rename(Donor_Group = group) %>%
  left_join(tax_clean, by = c("Donor_MAG" = "Name")) %>%
  rename(Donor_Phylum = Phylum) %>%
  # Join recipient groups and taxonomy
  left_join(group_data, by = c("Recipient_MAG" = "name")) %>%
  rename(Recipient_Group = group) %>%
  left_join(tax_clean, by = c("Recipient_MAG" = "Name")) %>%
  rename(Recipient_Phylum = Phylum) %>%
  # Join COG annotations using the first letter of COG_category
  left_join(donor_data %>% select(X.query, COG_category), by = c("donor_gene" = "X.query")) %>%
  mutate(COG_Main = str_sub(COG_category, 1, 1)) %>%
  # Handle missing values
  replace_na(list(COG_Main = "Unknown", Donor_Phylum = "Unclassified", Recipient_Phylum = "Unclassified"))

# ==========================================
# Part 2: Draw the multidimensional alluvial diagram
# ==========================================

# Count events for each flow path
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
# Part 3: Draw the community HGT interaction network
# ==========================================

# 1. Assign node phyla, retaining focal phyla and grouping the rest as "Others"
# Prepare network edges and node metadata
edges <- hgt_full %>% count(Donor_MAG, Recipient_MAG, name = "weight")
nodes <- data.frame(MAG = unique(c(edges$Donor_MAG, edges$Recipient_MAG))) %>%
  left_join(group_data, by = c("MAG" = "name")) %>%
  left_join(tax_clean, by = c("MAG" = "Name"))

# Build the igraph object
net <- graph_from_data_frame(d = edges, vertices = nodes, directed = TRUE)
# Use node in-degree and out-degree to set node sizes
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
  # Place Others last in the legend
  mutate(Phylum_Group = factor(Phylum_Group, levels = c(target_phyla, "Others")))

# 2. Build the igraph object and calculate node degree
net <- graph_from_data_frame(d = edges, vertices = nodes, directed = TRUE)
V(net)$degree <- degree(net, mode = "all")

# 3. Define custom color and shape mappings
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

# Shape 16 is a filled circle; shape 17 is a filled triangle. Assign these to E and R as needed
group_shapes <- c("E" = 16, "R" = 17) 

# 4. Draw the network using a Fruchterman-Reingold layout
p2_network_modified <- ggraph(net, layout = "fr") + 
  # Draw edges with reduced opacity to keep nodes visible
  geom_edge_fan(aes(edge_width = weight), alpha = 0.25, color = "gray60", 
                arrow = arrow(length = unit(2.5, 'mm')), end_cap = circle(3, 'mm')) +
  # Draw nodes using phylum colors, group shapes, and degree-based sizes
  geom_node_point(aes(color = Phylum_Group, shape = group, size = degree), alpha = 0.9) +
  # Apply the custom mappings
  scale_color_manual(values = phylum_colors) +
  scale_shape_manual(values = group_shapes) +
  scale_size_continuous(range = c(3, 11)) + # Increase the node-size range to emphasize hubs
  theme_void() +
  labs(color = "Phylum", shape = "Compartment", size = "HGT Degree") +
  # Arrange the legend for the main-text figure
  theme(
    legend.position = "right",
    legend.title = element_text(face = "bold", size = 12),
    legend.text = element_text(size = 10),
    legend.key = element_blank()
  ) +
  # Use one color in the size legend, independent of phylum colors
  guides(
    color = guide_legend(override.aes = list(size = 5)),
    shape = guide_legend(override.aes = list(size = 5)),
    size = guide_legend(override.aes = list(color = "gray30"))
  )

print(p2_network_modified)

# Export vector graphics for subsequent layout editing
ggsave("HGT_Network_Modified_FR.pdf", p2_network_modified, width = 11, height = 8)
