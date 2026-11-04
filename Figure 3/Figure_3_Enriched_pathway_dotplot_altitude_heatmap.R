#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(scales)
})

base_dir <- "/Volumes/Haleo_DK01/WD_Takakia_F/metagenome-data"
out_dir <- file.path(base_dir, "Figure4")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

tpm_file <- file.path(base_dir, "tpm_renormalized.txt")
enrich_file <- file.path(base_dir, "kegg_compare_enrichment_0_fdr_0612.csv")

blue <- "#82B3D7"
red <- "#CF6C72"

direction_label <- function(cluster) {
  dplyr::case_when(
    cluster == "Soil_Up" ~ "Rhizomatousphere enriched",
    cluster == "Plant_Up" ~ "Endophytic enriched",
    TRUE ~ cluster
  )
}

sig_symbol <- function(q) {
  dplyr::case_when(
    is.na(q) ~ "",
    q < 0.001 ~ "***",
    q < 0.01 ~ "**",
    q < 0.05 ~ "*",
    TRUE ~ ""
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
    arrange(Elevation, group_rank, Replicate)
}

paired_wilcox_safe <- function(x, y) {
  ok <- is.finite(x) & is.finite(y)
  x <- x[ok]
  y <- y[ok]
  if (length(x) < 2 || length(unique(x - y)) < 2) return(NA_real_)
  suppressWarnings(wilcox.test(x, y, paired = TRUE, exact = FALSE)$p.value)
}

friedman_safe <- function(df, value_col = "PathwayAbundance") {
  wide <- df %>%
    select(Replicate, Elevation, value = all_of(value_col)) %>%
    pivot_wider(names_from = Elevation, values_from = value) %>%
    arrange(Replicate)
  altitude_cols <- as.character(sort(unique(df$Elevation)))
  complete <- wide %>% select(all_of(altitude_cols)) %>% tidyr::drop_na()
  if (nrow(complete) < 2) {
    return(tibble(statistic = NA_real_, df = length(altitude_cols) - 1, p_value = NA_real_))
  }
  test <- suppressWarnings(friedman.test(as.matrix(complete)))
  tibble(
    statistic = unname(test$statistic),
    df = unname(test$parameter),
    p_value = test$p.value
  )
}

save_pair <- function(plot, stem, width, height) {
  ggsave(file.path(out_dir, paste0(stem, ".pdf")), plot, width = width, height = height,
         units = "in", device = "pdf", useDingbats = FALSE, bg = "white")
  ggsave(file.path(out_dir, paste0(stem, ".png")), plot, width = width, height = height,
         units = "in", dpi = 600, bg = "white")
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

# One heatmap block per pathway x elevation x habitat. Value is the mean of
# six biological replicates in that elevation-habitat cell.
pathway_altitude_mean <- pathway_abundance %>%
  group_by(Cluster, Direction, ID, Description, RichFactor, Count, MatchedKOs, Elevation, Group, Habitat) %>%
  summarise(
    mean_TPM = mean(PathwayAbundance, na.rm = TRUE),
    se_TPM = sd(PathwayAbundance, na.rm = TRUE) / sqrt(sum(!is.na(PathwayAbundance))),
    n = sum(!is.na(PathwayAbundance)),
    .groups = "drop"
  ) %>%
  mutate(
    Column = paste0(Elevation, "_", Group),
    Column = factor(Column, levels = as.vector(t(outer(c(3800, 3900, 4000, 4100, 4200), c("R", "E"), paste, sep = "_")))),
    Description = factor(as.character(Description), levels = display_levels)
  ) %>%
  group_by(Description) %>%
  mutate(
    log_mean_TPM = log10(mean_TPM + 1),
    row_z = as.numeric(scale(log_mean_TPM)),
    row_z = if_else(is.na(row_z), 0, row_z),
    row_z = pmax(pmin(row_z, 2.5), -2.5)
  ) %>%
  ungroup()

write.csv(
  pathway_altitude_mean,
  file.path(out_dir, "Figure4_pathway_altitude_RE_mean_heatmap_values.csv"),
  row.names = FALSE
)

# R/E paired Wilcoxon test within each pathway x elevation, BH-corrected across
# all pathway-elevation tests.
pathway_altitude_RE_wilcox <- pathway_abundance %>%
  select(Description, Direction, ID, Elevation, PairID, Group, PathwayAbundance) %>%
  distinct() %>%
  pivot_wider(names_from = Group, values_from = PathwayAbundance) %>%
  group_by(Description, Direction, ID, Elevation) %>%
  summarise(
    mean_R = mean(R, na.rm = TRUE),
    mean_E = mean(E, na.rm = TRUE),
    log2FC_R_over_E = log2((mean_R + 1e-6) / (mean_E + 1e-6)),
    p_value = paired_wilcox_safe(R, E),
    .groups = "drop"
  ) %>%
  mutate(
    q_value = p.adjust(p_value, method = "BH"),
    significance = sig_symbol(q_value),
    Description = factor(as.character(Description), levels = display_levels)
  )

write.csv(
  pathway_altitude_RE_wilcox,
  file.path(out_dir, "Figure4_pathway_altitude_RE_paired_wilcoxon_BH.csv"),
  row.names = FALSE
)

# Friedman tests within R and E separately across five altitudes for each
# pathway, BH-corrected across all pathway-habitat tests.
pathway_altitude_friedman <- pathway_abundance %>%
  group_by(Description, Direction, ID, Group, Habitat) %>%
  group_modify(function(.x, .g) friedman_safe(.x)) %>%
  ungroup() %>%
  mutate(
    q_value = p.adjust(p_value, method = "BH"),
    Description = factor(as.character(Description), levels = display_levels)
  )

write.csv(
  pathway_altitude_friedman,
  file.path(out_dir, "Figure4_pathway_altitude_within_group_friedman_BH.csv"),
  row.names = FALSE
)

# Pairwise Wilcoxon signed-rank tests among five altitudes within R and E,
# BH-corrected within each pathway-habitat family.
pathway_altitude_pairwise <- pathway_abundance %>%
  group_by(Description, Direction, ID, Group, Habitat) %>%
  group_modify(function(.x, .g) {
    elevs <- sort(unique(.x$Elevation))
    pairs <- combn(elevs, 2, simplify = FALSE)
    out <- lapply(pairs, function(pair) {
      wide <- .x %>%
        filter(Elevation %in% pair) %>%
        select(Replicate, Elevation, PathwayAbundance) %>%
        pivot_wider(names_from = Elevation, values_from = PathwayAbundance) %>%
        arrange(Replicate)
      p <- paired_wilcox_safe(wide[[as.character(pair[1])]], wide[[as.character(pair[2])]])
      tibble(Level1 = pair[1], Level2 = pair[2], p_value = p)
    })
    bind_rows(out) %>% mutate(q_value = p.adjust(p_value, method = "BH"))
  }) %>%
  ungroup()

write.csv(
  pathway_altitude_pairwise,
  file.path(out_dir, "Figure4_pathway_altitude_within_group_pairwise_wilcoxon_BH.csv"),
  row.names = FALSE
)

# New Figure4 exploratory version:
# Rhizomatousphere-enriched pathways are shown only for R across altitudes;
# Endophytic-enriched pathways are shown only for E across altitudes.
selected_group_mean <- pathway_altitude_mean %>%
  mutate(
    SelectedGroup = case_when(
      Direction == "Rhizomatousphere enriched" ~ "R",
      Direction == "Endophytic enriched" ~ "E",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(Group == SelectedGroup) %>%
  mutate(
    Altitude = factor(as.character(Elevation), levels = c("3800", "3900", "4000", "4100", "4200")),
    SelectedHabitat = if_else(SelectedGroup == "R", "Rhizomatousphere", "Endophytic")
  ) %>%
  group_by(Description) %>%
  mutate(
    log_mean_TPM = log10(mean_TPM + 1),
    row_z_selected = as.numeric(scale(log_mean_TPM)),
    row_z_selected = if_else(is.na(row_z_selected), 0, row_z_selected),
    row_z_selected = pmax(pmin(row_z_selected, 2.5), -2.5)
  ) %>%
  ungroup()

write.csv(
  selected_group_mean,
  file.path(out_dir, "Figure4_pathway_enriched_group_altitude_mean_heatmap_values.csv"),
  row.names = FALSE
)

selected_friedman <- pathway_altitude_friedman %>%
  mutate(
    SelectedGroup = case_when(
      Direction == "Rhizomatousphere enriched" ~ "R",
      Direction == "Endophytic enriched" ~ "E",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(Group == SelectedGroup) %>%
  group_by(Direction) %>%
  mutate(q_value_selected_direction = p.adjust(p_value, method = "BH")) %>%
  ungroup()

write.csv(
  selected_friedman,
  file.path(out_dir, "Figure4_pathway_enriched_group_altitude_friedman_BH.csv"),
  row.names = FALSE
)

selected_pairwise <- pathway_altitude_pairwise %>%
  mutate(
    SelectedGroup = case_when(
      Direction == "Rhizomatousphere enriched" ~ "R",
      Direction == "Endophytic enriched" ~ "E",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(Group == SelectedGroup) %>%
  group_by(Description, Direction, ID, Group, Habitat) %>%
  mutate(q_value_selected_pathway = p.adjust(p_value, method = "BH")) %>%
  ungroup()

write.csv(
  selected_pairwise,
  file.path(out_dir, "Figure4_pathway_enriched_group_altitude_pairwise_wilcoxon_BH.csv"),
  row.names = FALSE
)

friedman_label_df <- selected_friedman %>%
  mutate(
    Description = factor(as.character(Description), levels = display_levels),
    label = sig_symbol(q_value_selected_direction),
    label = if_else(label == "", "ns", label),
    Altitude = factor("4200", levels = c("3800", "3900", "4000", "4100", "4200"))
  )

p_heat <- ggplot(selected_group_mean, aes(x = Altitude, y = Description, fill = row_z_selected)) +
  geom_tile(width = 0.98, height = 0.92, colour = "white", linewidth = 0.18) +
  geom_hline(yintercept = 9.5, colour = "black", linewidth = 0.45) +
  geom_text(
    data = friedman_label_df,
    aes(x = Altitude, y = Description, label = label),
    inherit.aes = FALSE,
    family = "Helvetica",
    fontface = "bold",
    colour = "black",
    size = 2.8,
    nudge_x = 0.36
  ) +
  scale_fill_gradient2(
    low = "#3B7FB6",
    mid = "white",
    high = "#C95760",
    midpoint = 0,
    limits = c(-2.5, 2.5),
    name = "Row z-score\nlog10(mean TPM+1)"
  ) +
  scale_x_discrete(expand = expansion(mult = c(0, 0.06))) +
  scale_y_discrete(drop = FALSE) +
  labs(
    title = "Enriched pathway abundance across altitude",
    subtitle = "Blocks: mean TPM of six samples from the enriched habitat; right labels: Friedman BH.",
    x = "Altitude (m)",
    y = NULL
  ) +
  coord_cartesian(clip = "off") +
  theme_classic(base_size = 9, base_family = "Helvetica") +
  theme(
    axis.text.x = element_text(colour = "black", face = "bold", size = 8.5),
    axis.text.y = element_text(colour = "black", size = 8.2),
    axis.ticks = element_line(colour = "black", linewidth = 0.35),
    axis.line = element_line(colour = "black", linewidth = 0.45),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.55),
    plot.title = element_text(face = "bold", hjust = 0.5, size = 12),
    plot.subtitle = element_text(hjust = 0.5, size = 8.2),
    legend.title = element_text(face = "bold"),
    plot.margin = margin(18, 6, 6, 6)
  )

save_pair(p_heat, "Figure4_pathway_enriched_group_altitude_heatmap", width = 6.9, height = 6.4)

message("Done. Wrote altitude mean pathway heatmap and statistics to: ", out_dir)
