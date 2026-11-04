# Clear the workspace
rm(list = ls())

# Set the working directory
setwd("/Volumes/Haleo_DK01/WD_Takakia2")

# Load required packages
library(tidyverse)
library(ggplot2)
library(ggrepel)

# ------------------------------------------------------------------
# 1. Read and preprocess data
# ------------------------------------------------------------------

# Read the TPM data
tpm_df <- read.table("tpm_renormalized.txt", header = T, row.names = 1, sep = "\t", check.names = F)

# Infer groups from sample names: P = Plant, S = Soil
sample_names <- colnames(tpm_df)
condition <- ifelse(grepl("^P", sample_names), "Plant", "Soil")
colData <- data.frame(row.names = sample_names, condition = factor(condition, levels = c("Plant", "Soil")))

# Check group assignments
print(table(colData$condition))

# Log-transform TPM, adding 1 to handle zeros
# tpm_log <- log2(tpm_df + 1)
tpm_log <- tpm_df
# ------------------------------------------------------------------
# 2. Calculate differential abundance using t-tests or Wilcoxon tests
# ------------------------------------------------------------------

# Calculate log2 fold changes and p-values
calc_diff <- function(row, group1_samples, group2_samples, method = "wilcox") {
  group1 <- row[group1_samples]
  group2 <- row[group2_samples]
  log2fc <- log2((mean(group1) + 1e-6) / (mean(group2) + 1e-6))  # Add a small offset to avoid division by zero; Soil/Plant ratio
  if (method == "t") {
    test_res <- t.test(group1, group2)
  } else if (method == "wilcox") {
    test_res <- wilcox.test(group1, group2)
  }
  pval <- test_res$p.value
  return(c(log2FC = log2fc, pvalue = pval))
}

# Apply the test to each KO using tpm_log
plant_samples <- rownames(colData)[colData$condition == "Plant"]
soil_samples <- rownames(colData)[colData$condition == "Soil"]

res_list <- apply(tpm_log, 1, calc_diff, group1_samples = soil_samples, group2_samples = plant_samples, method = "wilcox")  # Use "t" to select a t-test
res_df <- as.data.frame(t(res_list))  # Transpose
res_df$KO_ID <- rownames(res_df)
res_df$padj <- p.adjust(res_df$pvalue, method = "BH")  # Benjamini-Hochberg correction
res_df$baseMean <- rowMeans(tpm_df)  # Mean TPM for optional filtering

# Save results
write.csv(res_df, "wilcox_results_0612.csv", row.names = FALSE)  # Alternatively, t_test_results.csv

# ------------------------------------------------------------------
# 3. Prepare plotting data
# ------------------------------------------------------------------

# Define thresholds
fc_cutoff <- 0      # Log2 fold-change threshold
p_cutoff <- 0.05    # Adjusted p-value threshold

# Assign categories: Up (higher in Soil), Down (higher in Plant), NS
res_df <- res_df %>%
  mutate(group = case_when(
    padj < p_cutoff & log2FC > fc_cutoff ~ "Up",    # Significantly higher in Soil
    padj < p_cutoff & log2FC < -fc_cutoff ~ "Down", # Significantly higher in Plant
    TRUE ~ "NS"
  ))

res_df$group <- factor(res_df$group, levels = c("Up", "Down", "NS"))

# Select the top 10 Up/Down KOs ranked by adjusted p-value
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
# 3. Define and annotate stress-response KO categories
# ------------------------------------------------------------------

# DNA repair
dna_repair_kos <- c("K03648", "K10771", "K03660", "K10847", "K10844", "K10843", "K08735", "K04482", "K10866", "K10865", "K10867", "K10873", "K10901", "K10605", "K08775", "K10875", "K10879", "K10880", "K10843", "K10849", "K10846", "K10848", "K10747", "K02330")

# Oxidative stress
oxidative_stress_kos <- c("K04565", "K04564", "K03781", "K23856", "K24158", "K00799", "K03671", "K00384", "K00131", "K00434", "K00430", "K08232", "K00128")

# Carotenoids
carotenoids_kos <- c("K02291", "K02293", "K00514", "K09835", "K06443", "K22502", "K02294", "K09838", "K09839", "K09840", "K10027", "K01669", "K06976")

# Amino acids, including cysteine and glutamate
amino_acids_kos <- c("K00640", "K01738", "K12339", "K10150", "K01697", "K00261", "K00262", "K01425", "K01919", "K11205", "K14454", "K14455", "K01915")

# Terpenoids and precursors
terpenoids_kos <- c("K00626", "K01641", "K00021", "K00869", "K00869", "K01597", "K01823", "K01662", "K00099", "K00991", "K00919", "K01770", "K03526", "K13789", "K00787")

# Flavonoids
flavonoids_kos <- c("K10775", "K00660", "K01859", "K13406", "K23135", "K13083", "K05278", "K00091", "K13082", "K05277", "K13081", "K21102", "K05281", "K13064")

# Cold stress
cold_stress_kos <- c("K11574", "K11575", "K11576", "K11482", "K03704")

# Siderophore-related KOs
iron_carrier_kos <- c("K02361", "K00216", "K01252", "K02364", "K08225", "K19611", "K25290")

# Trehalose-related KOs
trehalose_kos <- c("K02361", "K01087", "K01194", "K01194", "K00130", "K00286")

# Vitamin B12 and ACC
vb12_acc_kos <- c("K01505", "K02303", "K05934", "K05936", "K00595", "K02224", "K02232")

# Colonization and biofilm formation
colonization_kos <- c("K07173", "K10909", "K10911", "K10912", "K13061", "K18099", "K19421", "K19422", "K19423")

# List all categories
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
# 4. Plot with ggplot2
# ------------------------------------------------------------------

p <- ggplot(res_df, aes(x = log2FC, y = -log10(padj), color = color_group)) +
  # 1. Draw all points in the background
  geom_point(alpha = 0.6, size = 1.5) +
  
  # 2. Overlay stress-response KOs using a data subset
  geom_point(data = subset(res_df, category != "Other"), alpha = 0.8, size = 2) +  # Increase opacity and point size for emphasis
  
  # 3. Define group and category colors
  scale_color_manual(values = c(
    "Up" = "#01B7ED",          # Higher in Soil (blue)
    "Down" = "#37B484",        # Higher in Plant (green)
    "NS" = "#999999",          # Not significant (grey)
    "DNA_repair" = "purple",   # Purple
    "Oxidative_stress" = "red",# Red
    "Carotenoids" = "orange",  # Orange
    "Amino_acids" = "brown",   # Brown (including Cys and Glu)
    "Terpenoids" = "darkgreen",# Dark green
    "Flavonoids" = "yellow",   # Yellow
    "Cold_stress" = "cyan",    # Cyan
    "Iron_carrier" = "magenta",# Magenta
    "Trehalose" = "lightblue", # Light blue
    "VB12_ACC" = "violet",     # Violet
    "Colonization" = "darkred" # Dark red
  ),
  labels = c("Enriched in Soil", "Enriched in Plant", "NS", names(categories)),
  drop = FALSE) +  # Keep all factor levels in the legend, including levels without observations
  
  # 4. Add threshold reference lines
  geom_vline(xintercept = c(-fc_cutoff, fc_cutoff), linetype = "dashed", color = "black", alpha=0.5) +
  geom_hline(yintercept = -log10(p_cutoff), linetype = "dashed", color = "black", alpha=0.5) +
  
  # 5. Label top genes with ggrepel
  geom_text_repel(
    data = top_genes,
    aes(label = KO_ID),
    size = 3,
    color = "black",      # Use black for all labels
    box.padding = 0.5,    # Padding around label boxes
    max.overlaps = 50,    # Maximum number of overlaps
    segment.color = "grey50" # Connector line color
  ) +
  
  # 6. Set axis labels and titles
  labs(x = expression(log[2] ~ "Fold Change"), 
       y = expression(-log[10] ~ "Adjusted P-value"),
       title = "Soil vs Plant (KO Abundance)",
       subtitle = paste0("Up: ", sum(res_df$group == "Up"), 
                         " | Down: ", sum(res_df$group == "Down")),
       color = "Group/Category") +  # Set the legend title
  
  # 7. Apply the publication-style classic theme
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
# 5. Save results
# ------------------------------------------------------------------

# Display the plot
print(p)


# Save the PDF
ggsave("Volcano_Plot_wilcox_12_0612.pdf", plot = p, width = 6, height = 5)
