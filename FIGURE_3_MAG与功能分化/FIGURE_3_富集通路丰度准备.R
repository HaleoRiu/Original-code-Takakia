#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(patchwork)
  library(scales)
})

base_dir <- "/Volumes/Haleo_DK01/WD_Takakia_F/metagenome-data"
out_dir <- file.path(base_dir, "Figure4")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

tpm_file <- file.path(base_dir, "tpm_renormalized.txt")
enrich_file <- file.path(base_dir, "kegg_compare_enrichment_0_fdr_0612.csv")

blue <- "#82B3D7"
red <- "#CF6C72"
group_cols <- c("Rhizomatousphere" = blue, "Endophytic" = red)
direction_cols <- c("Rhizomatousphere enriched" = blue, "Endophytic enriched" = red)

direction_label <- function(cluster) {
  dplyr::case_when(
    cluster == "Soil_Up" ~ "Rhizomatousphere enriched",
    cluster == "Plant_Up" ~ "Endophytic enriched",
    TRUE ~ cluster
  )
}

sig_symbol <- function(q) {
  dplyr::case_when(
    is.na(q) ~ "NA",
    q < 0.001 ~ "***",
    q < 0.01 ~ "**",
    q < 0.05 ~ "*",
    TRUE ~ "ns"
  )
}

sample_info <- function(samples) {
  tibble(sample = samples) %>%
    filter(grepl("^[PS][0-9]{2}_[0-9]+$", sample)) %>%
    mutate(
      Prefix = substr(sample, 1, 1),
      Group = if_else(Prefix == "S", "R", "E"),
      Habitat = if_else(Group == "R", "Rhizomatousphere", "Endophytic"),
      PairID = sub("^[PS]", "", sample),
      Elevation = as.integer(substr(sample, 2, 3)) * 100,
      Replicate = as.integer(sub("^[PS][0-9]{2}_", "", sample)),
      group_rank = if_else(Group == "R", 1L, 2L)
    ) %>%
    arrange(group_rank, Elevation, Replicate) %>%
    mutate(SampleIndex = row_number())
}

save_pair <- function(plot, stem, width, height) {
  ggsave(file.path(out_dir, paste0(stem, ".pdf")), plot, width = width, height = height,
         units = "in", device = "pdf", useDingbats = FALSE, bg = "white")
  ggsave(file.path(out_dir, paste0(stem, ".png")), plot, width = width, height = height,
         units = "in", dpi = 600, bg = "white")
}

theme_panel <- function(base_size = 9) {
  theme_classic(base_size = base_size, base_family = "Helvetica") +
    theme(
      axis.text = element_text(colour = "black"),
      axis.title = element_text(face = "bold", colour = "black"),
      axis.line = element_line(colour = "black", linewidth = 0.55),
      axis.ticks = element_line(colour = "black", linewidth = 0.5),
      panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.65),
      plot.title = element_text(face = "bold", hjust = 0.5, size = base_size + 1.5),
      legend.title = element_text(face = "bold"),
      plot.margin = margin(4, 2, 4, 2)
    )
}

tpm_raw <- read.delim(tpm_file, check.names = FALSE)
sample_cols0 <- names(tpm_raw)[grepl("^[PS][0-9]{2}_[0-9]+$", names(tpm_raw))]
meta <- sample_info(sample_cols0)
sample_cols <- meta$sample

tpm <- tpm_raw %>%
  rename(KO = `feature-id`) %>%
  mutate(KO = as.character(KO)) %>%
  select(KO, all_of(sample_cols))

enrich <- read.csv(enrich_file, check.names = FALSE) %>%
  mutate(
    Direction = direction_label(Cluster),
    Direction = factor(Direction, levels = c("Rhizomatousphere enriched", "Endophytic enriched")),
    Description = as.character(Description),
    GeneList = strsplit(as.character(geneID), "/")
  ) %>%
  arrange(Direction, qvalue, p.adjust)

pathway_order <- enrich$Description
display_levels <- rev(pathway_order)

pathway_abundance <- lapply(seq_len(nrow(enrich)), function(i) {
  genes <- unique(enrich$GeneList[[i]])
  matched <- genes[genes %in% tpm$KO]
  vals <- if (length(matched) == 0) {
    rep(0, length(sample_cols))
  } else {
    colSums(as.matrix(tpm[tpm$KO %in% matched, sample_cols, drop = FALSE]), na.rm = TRUE)
  }
  tibble(
    Cluster = enrich$Cluster[i],
    Direction = enrich$Direction[i],
    ID = enrich$ID[i],
    Description = enrich$Description[i],
    RichFactor = enrich$RichFactor[i],
    Count = enrich$Count[i],
    qvalue_enrichment = enrich$qvalue[i],
    sample = sample_cols,
    PathwayAbundance = as.numeric(vals),
    MatchedKOs = length(matched)
  )
}) %>%
  bind_rows() %>%
  left_join(meta, by = "sample") %>%
  mutate(
    Description = factor(Description, levels = display_levels),
    Direction = factor(Direction, levels = c("Rhizomatousphere enriched", "Endophytic enriched"))
  )

pathway_stats <- pathway_abundance %>%
  select(Description, Direction, ID, RichFactor, Count, MatchedKOs, PairID, Group, PathwayAbundance) %>%
  distinct() %>%
  pivot_wider(names_from = Group, values_from = PathwayAbundance) %>%
  group_by(Description, Direction, ID, RichFactor, Count, MatchedKOs) %>%
  summarise(
    mean_R = mean(R, na.rm = TRUE),
    mean_E = mean(E, na.rm = TRUE),
    se_R = sd(R, na.rm = TRUE) / sqrt(sum(!is.na(R))),
    se_E = sd(E, na.rm = TRUE) / sqrt(sum(!is.na(E))),
    log2FC_R_over_E = log2((mean_R + 1e-6) / (mean_E + 1e-6)),
    p_value = wilcox.test(R, E, paired = TRUE, exact = FALSE)$p.value,
    .groups = "drop"
  ) %>%
  mutate(
    q_value = p.adjust(p_value, method = "BH"),
    significance = sig_symbol(q_value),
    Description = factor(as.character(Description), levels = display_levels)
  )

write.csv(pathway_abundance, file.path(out_dir, "Figure4_integrated_pathway_abundance_by_sample.csv"), row.names = FALSE)
write.csv(pathway_stats, file.path(out_dir, "Figure4_integrated_pathway_stats.csv"), row.names = FALSE)

bar_df <- pathway_stats %>%
  select(Description, Direction, mean_R, mean_E, se_R, se_E, significance) %>%
  pivot_longer(c(mean_R, mean_E), names_to = "MeanType", values_to = "MeanTPM") %>%
  mutate(
    Habitat = if_else(MeanType == "mean_R", "Rhizomatousphere", "Endophytic"),
    SE = if_else(MeanType == "mean_R", se_R, se_E),
    Habitat = factor(Habitat, levels = c("Rhizomatousphere", "Endophytic"))
  )

label_df <- pathway_stats %>%
  mutate(x_label = pmax(mean_R, mean_E) + max(c(mean_R, mean_E), na.rm = TRUE) * 0.025)

heatmap_df <- pathway_abundance %>%
  group_by(Description) %>%
  mutate(
    log_abundance = log10(PathwayAbundance + 1),
    row_z = as.numeric(scale(log_abundance)),
    row_z = pmax(pmin(row_z, 2.5), -2.5)
  ) %>%
  ungroup() %>%
  mutate(
    SampleIndex = as.numeric(SampleIndex),
    Description = factor(as.character(Description), levels = display_levels)
  )

altitude_centres <- meta %>%
  group_by(Group, Elevation) %>%
  summarise(Centre = mean(SampleIndex), .groups = "drop") %>%
  arrange(if_else(Group == "R", 1L, 2L), Elevation) %>%
  mutate(Label = as.character(Elevation))

direction_boundary <- which(display_levels == "Aminoacyl-tRNA biosynthesis") + 0.5

p_names <- ggplot(pathway_stats, aes(y = Description)) +
  geom_text(aes(x = 1, label = Description), hjust = 1, family = "Helvetica", size = 3.0) +
  scale_y_discrete(drop = FALSE) +
  scale_x_continuous(limits = c(0, 1.02), expand = expansion(mult = c(0, 0))) +
  labs(title = "Pathway", x = NULL, y = NULL) +
  theme_void(base_family = "Helvetica") +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5, size = 10.5),
    plot.margin = margin(18, 0, 22, 2)
  )

p_bar <- ggplot(bar_df, aes(x = MeanTPM, y = Description, fill = Habitat, group = Habitat)) +
  geom_col(
    position = position_dodge(width = 0.82),
    width = 0.72,
    colour = "black",
    linewidth = 0.25
  ) +
  geom_errorbar(
    aes(xmin = pmax(MeanTPM - SE, 0), xmax = MeanTPM + SE),
    position = position_dodge(width = 0.82),
    width = 0.28,
    orientation = "y",
    linewidth = 0.35,
    colour = "black"
  ) +
  geom_text(
    data = label_df,
    aes(x = x_label, y = Description, label = significance),
    inherit.aes = FALSE,
    hjust = 0,
    family = "Helvetica",
    fontface = "bold",
    colour = "red",
    size = 3.0
  ) +
  scale_fill_manual(values = group_cols, name = NULL) +
  scale_y_discrete(drop = FALSE) +
  scale_x_continuous(labels = comma, expand = expansion(mult = c(0, 0.10))) +
  labs(title = "Pathway abundance", x = "TPM", y = NULL) +
  theme_panel(9) +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    legend.position = "top",
    legend.text = element_text(size = 8.2),
    legend.key.width = grid::unit(0.45, "cm")
  )

p_log2fc <- ggplot(pathway_stats, aes(x = "log2FC", y = Description, fill = log2FC_R_over_E)) +
  geom_tile(colour = "#A8A8A8", linewidth = 0.45, width = 0.95, height = 0.90) +
  geom_text(aes(label = sprintf("%.2f", log2FC_R_over_E)), family = "Helvetica", size = 3.0, colour = "black") +
  scale_fill_gradient2(low = blue, mid = "white", high = red, midpoint = 0, name = expression(log[2]*FC)) +
  scale_y_discrete(drop = FALSE) +
  labs(title = expression(log[2]*FC), x = NULL, y = NULL) +
  theme_void(base_family = "Helvetica") +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5, size = 10.5),
    legend.position = "bottom",
    legend.key.width = grid::unit(0.65, "cm"),
    legend.key.height = grid::unit(0.20, "cm"),
    plot.margin = margin(18, 2, 22, 2)
  )

p_heat <- ggplot(heatmap_df, aes(x = SampleIndex, y = Description, fill = row_z)) +
  geom_tile(width = 1, height = 0.92, colour = NA) +
  geom_vline(xintercept = c(6.5, 12.5, 18.5, 24.5, 30.5, 36.5, 42.5, 48.5, 54.5), colour = "#D0D0D0", linewidth = 0.25) +
  geom_vline(xintercept = 30.5, colour = "black", linewidth = 0.45) +
  scale_fill_gradient2(
    low = "#3B7FB6",
    mid = "white",
    high = "#C95760",
    midpoint = 0,
    limits = c(-2.5, 2.5),
    name = "Row z-score\nlog10(TPM+1)"
  ) +
  scale_x_continuous(
    breaks = altitude_centres$Centre,
    labels = altitude_centres$Label,
    expand = expansion(mult = c(0, 0), add = c(0.5, 0.5))
  ) +
  scale_y_discrete(drop = FALSE) +
  labs(title = "Samples", x = NULL, y = NULL) +
  theme_panel(9) +
  theme(
    axis.text.x = element_text(face = "bold", angle = 45, hjust = 1, size = 8.2),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    legend.position = "right"
  )

integrated <- p_names + p_bar + p_log2fc + p_heat +
  plot_layout(widths = c(1.55, 1.8, 0.72, 2.2), guides = "keep") +
  plot_annotation(tag_levels = list(c("(c)"))) &
  theme(plot.tag = element_text(face = "bold", size = 16, family = "Helvetica"))

save_pair(integrated, "Figure4c_integrated_pathway_abundance_log2FC_heatmap", width = 13.2, height = 6.1)

# KO contribution analysis ------------------------------------------------------
ko_contrib <- lapply(seq_len(nrow(enrich)), function(i) {
  genes <- unique(enrich$GeneList[[i]])
  matched <- genes[genes %in% tpm$KO]
  if (length(matched) == 0) return(NULL)
  long <- tpm %>%
    filter(KO %in% matched) %>%
    pivot_longer(all_of(sample_cols), names_to = "sample", values_to = "TPM") %>%
    left_join(meta, by = "sample")
  stats <- long %>%
    group_by(KO) %>%
    summarise(
      mean_R = mean(TPM[Group == "R"], na.rm = TRUE),
      mean_E = mean(TPM[Group == "E"], na.rm = TRUE),
      delta_R_minus_E = mean_R - mean_E,
      abs_delta = abs(delta_R_minus_E),
      .groups = "drop"
    )
  total_abs <- sum(stats$abs_delta, na.rm = TRUE)
  stats %>%
    mutate(
      Cluster = enrich$Cluster[i],
      Direction = enrich$Direction[i],
      ID = enrich$ID[i],
      Description = enrich$Description[i],
      log2FC_R_over_E_KO = log2((mean_R + 1e-6) / (mean_E + 1e-6)),
      contribution_percent = if (total_abs > 0) 100 * abs_delta / total_abs else 0
    )
}) %>%
  bind_rows() %>%
  group_by(Description) %>%
  arrange(desc(contribution_percent), .by_group = TRUE) %>%
  mutate(rank_in_pathway = row_number()) %>%
  ungroup()

write.csv(ko_contrib, file.path(out_dir, "Figure4_KO_contribution_to_pathway_change.csv"), row.names = FALSE)

top_contrib <- ko_contrib %>%
  group_by(Description) %>%
  slice_max(contribution_percent, n = 3, with_ties = FALSE) %>%
  ungroup() %>%
  mutate(
    PathwayShort = as.character(Description),
    KO = factor(KO, levels = rev(unique(KO[order(Description, contribution_percent)])))
  )

make_contrib_column <- function(df, title_text) {
  ggplot(df, aes(x = contribution_percent, y = reorder(KO, contribution_percent), fill = Direction)) +
    geom_col(width = 0.68, colour = "black", linewidth = 0.25) +
    facet_wrap(~ Description, scales = "free_y", ncol = 1) +
    scale_fill_manual(values = direction_cols, guide = "none") +
    labs(
      title = title_text,
      x = NULL,
      y = NULL
    ) +
    theme_panel(9) +
    theme(
      strip.background = element_rect(fill = "#E6E6E6", colour = "black", linewidth = 0.55),
      strip.text = element_text(face = "bold", size = 7.2),
      axis.text.y = element_text(size = 7.5),
      plot.title = element_text(face = "bold", hjust = 0.5, size = 10.5)
    )
}

p_contrib_R <- make_contrib_column(
  top_contrib %>% filter(Direction == "Rhizomatousphere enriched"),
  "Rhizomatousphere enriched"
)

p_contrib_E <- make_contrib_column(
  top_contrib %>% filter(Direction == "Endophytic enriched"),
  "Endophytic enriched"
)

p_contrib <- p_contrib_R + p_contrib_E +
  plot_layout(widths = c(1, 1)) +
  plot_annotation(
    title = "Top KO contributors to pathway abundance differences",
    caption = "Contribution to absolute pathway difference (%)"
  ) &
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5, size = 13, family = "Helvetica"),
    plot.caption = element_text(face = "bold", hjust = 0.5, size = 10, family = "Helvetica")
  )

save_pair(p_contrib, "Figure4_top_KO_contributors_by_pathway", width = 9.8, height = 12.2)

message("Done. Integrated Figure4 pathway draft written to: ", out_dir)
