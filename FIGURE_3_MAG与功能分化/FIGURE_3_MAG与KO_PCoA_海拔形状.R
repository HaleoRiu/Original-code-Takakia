#!/usr/bin/env Rscript

# Resizable MAG PCoA and functional-gene PCoA panels.
# Key point: this script deliberately avoids theme(aspect.ratio = ...)
# and coord_fixed(), so the plotting panel follows the width/height
# specified in save_pair().

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(vegan)
})

base_dir <- "/Volumes/Haleo_DK01/WD_Takakia_F/metagenome-data"
out_dir <- file.path(base_dir, "metagenome_redraw")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

mags_file <- file.path(base_dir, "MAGs_with_taxonomy_70.txt")
tpm_file <- file.path(base_dir, "tpm_renormalized.txt")

# -----------------------------------------------------------------------------
# Change the width/height here. The panel will follow these settings.
# -----------------------------------------------------------------------------
MAG_PCOA_WIDTH <- 3.12
MAG_PCOA_HEIGHT <- 2.9
FUNC_PCA_WIDTH <- 3.2
FUNC_PCA_HEIGHT <- 2.9

group_full_colours <- c("Rhizomatousphere" = "#82B3D7", "Endophytic" = "#CF6C72")
group_full_labels <- c(
  "Rhizomatousphere" = "Rhizomatousphere, R",
  "Endophytic" = "Endophytic, E"
)

altitude_shapes <- c(
  "3800" = 24,
  "3900" = 22,
  "4000" = 21,
  "4100" = 23,
  "4200" = 25
)

format_p <- function(p) {
  ifelse(p < 0.001, "0.001", sprintf("%.3f", p))
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
      AltitudeShape = factor(Elevation, levels = c(3800, 3900, 4000, 4100, 4200)),
      Replicate = as.integer(sub("^[PS][0-9]{2}_", "", sample)),
      Group = factor(Group, levels = c("R", "E")),
      Habitat = factor(Habitat, levels = c("Rhizomatousphere", "Endophytic"))
    )
}

save_pair <- function(plot, stem, width, height) {
  png_file <- file.path(out_dir, paste0(stem, ".png"))
  pdf_file <- file.path(out_dir, paste0(stem, ".pdf"))

  ggsave(
    png_file,
    plot = plot,
    width = width,
    height = height,
    units = "in",
    dpi = 600,
    bg = "white"
  )
  ggsave(
    pdf_file,
    plot = plot,
    width = width,
    height = height,
    units = "in",
    device = grDevices::pdf,
    bg = "white"
  )
}

ordination_theme_resizable <- function(base_size = 10) {
  theme_classic(base_size = base_size, base_family = "Helvetica") +
    theme(
      panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.8),
      axis.line = element_blank(),
      axis.ticks = element_line(colour = "black", linewidth = 0.55),
      axis.ticks.length = grid::unit(0.11, "cm"),
      axis.text = element_text(colour = "black", size = base_size - 1),
      axis.title = element_text(face = "bold", size = base_size + 1),
      plot.title = element_blank(),
      plot.margin = margin(4, 4, 4, 4),
      legend.position = "none"
    )
}

make_ordination_plot_resizable <- function(df, x_col, y_col, pc, adonis_stats,
                                           x_prefix = "PCoA1", y_prefix = "PCoA2",
                                           base_size = 10) {
  x <- df[[x_col]]
  y <- df[[y_col]]
  x_range <- diff(range(x, na.rm = TRUE))
  y_range <- diff(range(y, na.rm = TRUE))
  x_left <- min(x, na.rm = TRUE) + x_range * 0.04
  y_top <- max(y, na.rm = TRUE) - y_range * 0.04

  ggplot(df, aes(x = .data[[x_col]], y = .data[[y_col]])) +
    stat_ellipse(
      aes(fill = Habitat),
      type = "t",
      level = 0.95,
      geom = "polygon",
      alpha = 0.22,
      colour = NA,
      show.legend = FALSE
    ) +
    stat_ellipse(
      aes(colour = Habitat),
      type = "t",
      level = 0.95,
      linewidth = 0.9,
      linetype = "solid",
      show.legend = FALSE
    ) +
    geom_point(
      aes(fill = Habitat, colour = Habitat, shape = AltitudeShape),
      size = 2.0,
      stroke = 0.30,
      alpha = 0.68,
      show.legend = FALSE
    ) +
    annotate(
      "text",
      x = x_left,
      y = y_top,
      label = paste0("Adonis ", "R\u00B2", " = ", sprintf("%.3f", adonis_stats$R2)),
      hjust = 0,
      vjust = 1,
      family = "Helvetica",
      size = base_size * 0.35
    ) +
    annotate(
      "text",
      x = x_left,
      y = y_top - y_range * 0.10,
      label = "p =",
      hjust = 0,
      vjust = 1,
      family = "Helvetica",
      fontface = "italic",
      size = base_size * 0.35
    ) +
    annotate(
      "text",
      x = x_left + x_range * 0.13,
      y = y_top - y_range * 0.10,
      label = format_p(adonis_stats$p_value),
      hjust = 0,
      vjust = 1,
      family = "Helvetica",
      colour = "red",
      size = base_size * 0.35
    ) +
    scale_fill_manual(values = group_full_colours, breaks = names(group_full_labels), drop = FALSE) +
    scale_colour_manual(values = group_full_colours, breaks = names(group_full_labels), drop = FALSE) +
    scale_shape_manual(values = altitude_shapes, breaks = names(altitude_shapes), drop = FALSE) +
    labs(
      x = paste0(x_prefix, " (", sprintf("%.2f", pc[1]), "%)"),
      y = paste0(y_prefix, " (", sprintf("%.2f", pc[2]), "%)")
    ) +
    coord_cartesian(clip = "on") +
    ordination_theme_resizable(base_size = base_size)
}

# -----------------------------------------------------------------------------
# MAG-level PCoA
# -----------------------------------------------------------------------------
mags <- read.delim(mags_file, check.names = FALSE)
mag_sample_cols <- names(mags)[grepl("^[PS][0-9]{2}_[0-9]+$", names(mags))]
mag_meta <- sample_info(mag_sample_cols)

mag_mat <- mags %>%
  select(Genome, all_of(mag_sample_cols)) %>%
  tibble::column_to_rownames("Genome") %>%
  as.matrix()
storage.mode(mag_mat) <- "numeric"
mag_sample_mat <- t(mag_mat)
mag_sample_mat <- mag_sample_mat[mag_meta$sample, , drop = FALSE]

mag_dist <- vegdist(mag_sample_mat, method = "bray")
mag_pcoa <- cmdscale(mag_dist, eig = TRUE, k = 2)
mag_eig <- mag_pcoa$eig
mag_var <- mag_eig[seq_len(2)] / sum(mag_eig[mag_eig > 0]) * 100
mag_scores <- as.data.frame(mag_pcoa$points)
colnames(mag_scores) <- c("PCoA1", "PCoA2")
mag_scores$sample <- rownames(mag_scores)
mag_scores <- left_join(mag_scores, mag_meta, by = "sample")

set.seed(20260625)
mag_adonis <- adonis2(mag_dist ~ Group, data = mag_meta, permutations = 999)
mag_adonis_stats <- tibble(R2 = mag_adonis$R2[1], p_value = mag_adonis$`Pr(>F)`[1])

p_mag_pcoa_resizable <- make_ordination_plot_resizable(
  df = mag_scores,
  x_col = "PCoA1",
  y_col = "PCoA2",
  pc = mag_var,
  adonis_stats = mag_adonis_stats,
  x_prefix = "PCoA1",
  y_prefix = "PCoA2",
  base_size = 10
)

save_pair(
  p_mag_pcoa_resizable,
  "02_MAG_level_PCoA_R_E_resizable_wide_test_altitude_shape",
  width = MAG_PCOA_WIDTH,
  height = MAG_PCOA_HEIGHT
)

# -----------------------------------------------------------------------------
# Functional-gene PCoA based on Bray-Curtis distance
# -----------------------------------------------------------------------------
tpm <- read.delim(tpm_file, check.names = FALSE)
tpm_sample_cols <- names(tpm)[grepl("^[PS][0-9]{2}_[0-9]+$", names(tpm))]
tpm_meta <- sample_info(tpm_sample_cols)

tpm_mat <- tpm %>%
  rename(feature = `feature-id`) %>%
  select(feature, all_of(tpm_sample_cols)) %>%
  tibble::column_to_rownames("feature") %>%
  as.matrix()
storage.mode(tpm_mat) <- "numeric"
tpm_sample_mat <- t(tpm_mat)
tpm_sample_mat <- tpm_sample_mat[tpm_meta$sample, , drop = FALSE]
keep <- apply(tpm_sample_mat, 2, sd, na.rm = TRUE) > 0
tpm_sample_mat <- tpm_sample_mat[, keep, drop = FALSE]

func_dist <- vegdist(tpm_sample_mat, method = "bray")
func_pcoa <- cmdscale(func_dist, eig = TRUE, k = 2)
func_eig <- func_pcoa$eig
func_var <- func_eig[seq_len(2)] / sum(func_eig[func_eig > 0]) * 100
func_scores <- as.data.frame(func_pcoa$points)
colnames(func_scores) <- c("PCoA1", "PCoA2")
func_scores$sample <- rownames(func_scores)
func_scores <- left_join(func_scores, tpm_meta, by = "sample")

set.seed(20260625)
func_adonis <- adonis2(func_dist ~ Group, data = tpm_meta, permutations = 999)
func_adonis_stats <- tibble(R2 = func_adonis$R2[1], p_value = func_adonis$`Pr(>F)`[1])

p_func_pcoa_resizable <- make_ordination_plot_resizable(
  df = func_scores,
  x_col = "PCoA1",
  y_col = "PCoA2",
  pc = func_var,
  adonis_stats = func_adonis_stats,
  x_prefix = "PCoA1",
  y_prefix = "PCoA2",
  base_size = 10
)

save_pair(
  p_func_pcoa_resizable,
  "04_functional_gene_PCoA_Bray_R_E_resizable_wide_test_altitude_shape",
  width = FUNC_PCA_WIDTH,
  height = FUNC_PCA_HEIGHT
)

message("Done. Resizable ordination outputs written to: ", out_dir)
