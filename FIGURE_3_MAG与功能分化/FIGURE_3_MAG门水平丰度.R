#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(patchwork)
  library(vegan)
  library(scales)
})

base_dir <- "/Volumes/Haleo_DK01/WD_Takakia_F/metagenome-data"
out_dir <- file.path(base_dir, "metagenome_redraw")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

mags_file <- file.path(base_dir, "MAGs_with_taxonomy_70.txt")
div_file <- file.path(base_dir, "div.txt")
tpm_file <- file.path(base_dir, "tpm_renormalized.txt")

group_colours <- c("R" = "#82B3D7", "E" = "#CF6C72")
group_outline <- c("R" = "#6EA7D0", "E" = "#C95760")
group_full_colours <- c("Rhizomatousphere" = "#82B3D7", "Endophytic" = "#CF6C72")
group_full_labels <- c(
  "Rhizomatousphere" = "Rhizomatousphere, R",
  "Endophytic" = "Endophytic, E"
)

phylum_colours <- c(
  "Pseudomonadota" = "#F2B8B5",
  "Myxococcota" = "#BFD8F0",
  "Acidobacteriota" = "#BFECC4",
  "Bacteroidota" = "#E9D2B7",
  "Gemmatimonadota" = "#CBBCE7",
  "Actinomycetota" = "#EBC2DA",
  "Patescibacteria" = "#B7D9E8",
  "Chloroflexota" = "#F3CCCC",
  "Eremiobacterota" = "#CFEAF0",
  "Others" = "#666666"
)
phylum_order <- names(phylum_colours)

format_p <- function(p) {
  ifelse(p < 0.001, "0.001", sprintf("%.3f", p))
}

sig_symbol <- function(q) {
  case_when(
    q < 0.001 ~ "***",
    q < 0.01 ~ "**",
    q < 0.05 ~ "*",
    TRUE ~ "NS"
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
      Group = factor(Group, levels = c("R", "E")),
      Habitat = factor(Habitat, levels = c("Rhizomatousphere", "Endophytic"))
    )
}

clean_phylum <- function(x) {
  x <- sub("^p__", "", as.character(x))
  x[is.na(x) | x == "" | x == "NA"] <- "Others"
  ifelse(x %in% setdiff(phylum_order, "Others"), x, "Others")
}

save_pair <- function(plot, stem, width, height) {
  ggsave(
    file.path(out_dir, paste0(stem, ".pdf")),
    plot,
    width = width,
    height = height,
    units = "in",
    device = "pdf",
    useDingbats = FALSE,
    bg = "white"
  )
  ggsave(
    file.path(out_dir, paste0(stem, ".png")),
    plot,
    width = width,
    height = height,
    units = "in",
    dpi = 600,
    bg = "white"
  )
}

# -----------------------------------------------------------------------------
# Shared styles copied/adapted from the amplicon redraws
# -----------------------------------------------------------------------------
pcoa_theme <- function(base_size = 11) {
  theme_bw(base_size = base_size, base_family = "Helvetica") +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.95),
      axis.title = element_text(face = "bold", size = base_size + 1),
      axis.text = element_text(colour = "black", size = base_size - 1),
      plot.title = element_text(face = "bold", hjust = 0.5, size = base_size + 2),
      legend.position = "bottom",
      legend.title = element_blank(),
      legend.text = element_text(size = base_size),
      legend.key.width = grid::unit(0.55, "cm"),
      plot.margin = margin(5, 5, 5, 5),
      panel.ontop = FALSE
    )
}

make_ordination_plot <- function(df, x_col, y_col, pc, adonis_stats, title,
                                 x_prefix = "PCoA1", y_prefix = "PCoA2",
                                 show_legend = TRUE, base_size = 11) {
  x <- df[[x_col]]
  y <- df[[y_col]]
  x_range <- diff(range(x, na.rm = TRUE))
  y_range <- diff(range(y, na.rm = TRUE))
  x_left <- min(x, na.rm = TRUE) + x_range * 0.03
  y_top <- max(y, na.rm = TRUE) - y_range * 0.03

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
      linewidth = 1.05,
      linetype = "solid",
      show.legend = FALSE
    ) +
    geom_point(
      aes(fill = Habitat, colour = Habitat),
      shape = 21,
      size = 2.35,
      stroke = 0.35,
      alpha = 0.95
    ) +
    annotate(
      "rect",
      xmin = x_left - x_range * 0.01,
      xmax = x_left + x_range * 0.42,
      ymin = y_top - y_range * 0.20,
      ymax = y_top + y_range * 0.03,
      fill = "white",
      alpha = 0.72,
      colour = NA
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
    scale_fill_manual(
      values = group_full_colours,
      breaks = names(group_full_labels),
      labels = group_full_labels,
      drop = FALSE
    ) +
    scale_colour_manual(
      values = group_full_colours,
      breaks = names(group_full_labels),
      labels = group_full_labels,
      drop = FALSE
    ) +
    guides(
      fill = guide_legend(
        override.aes = list(shape = 22, size = 5, colour = "black", alpha = 1)
      ),
      colour = "none"
    ) +
    labs(
      title = title,
      x = paste0(x_prefix, " (", sprintf("%.2f", pc[1]), "%)"),
      y = paste0(y_prefix, " (", sprintf("%.2f", pc[2]), "%)")
    ) +
    coord_cartesian(clip = "on") +
    pcoa_theme(base_size) +
    theme(
      legend.position = if (show_legend) "bottom" else "none",
      aspect.ratio = 0.82
    )
}

make_title_strip <- function(title, left_margin = 26, text_size = 5.0) {
  ggplot() +
    annotate(
      "rect",
      xmin = 0,
      xmax = 1,
      ymin = 0,
      ymax = 1,
      fill = "#E6E6E6",
      colour = "black",
      linewidth = 0.9
    ) +
    annotate(
      "text",
      x = 0.5,
      y = 0.5,
      label = title,
      family = "Helvetica",
      fontface = "bold",
      size = text_size
    ) +
    coord_cartesian(xlim = c(0, 1), ylim = c(0, 1), expand = FALSE) +
    theme_void() +
    theme(plot.margin = margin(0, 0, 0, left_margin))
}

make_alpha_body <- function(alpha, metric_name, box_width = 0.34, dodge_width = 0.52) {
  elevation_levels <- c("All", paste0(c(3800, 3900, 4000, 4100, 4200), " m"))
  elevation_positions <- c("All" = 1, "3800 m" = 2, "3900 m" = 3,
                           "4000 m" = 4, "4100 m" = 5, "4200 m" = 6)
  df <- alpha$data %>%
    filter(as.character(Metric) == metric_name) %>%
    mutate(
      ElevationGroup = factor(ElevationGroup, levels = elevation_levels),
      Section = if_else(as.character(ElevationGroup) == "All", "All", "Elevation"),
      Section = factor(Section, levels = c("All", "Elevation")),
      ElevationX = unname(elevation_positions[as.character(ElevationGroup)])
    )

  y_extent <- df %>%
    summarise(y_min = min(Value, na.rm = TRUE), y_max = max(Value, na.rm = TRUE)) %>%
    mutate(
      y_range = y_max - y_min,
      y_range = if_else(y_range == 0, y_max, y_range),
      y_lower = y_min - y_range * 0.10,
      y_upper = y_max + y_range * 0.30
    )

  ranges <- df %>%
    group_by(Metric, ElevationGroup) %>%
    summarise(
      y_min = min(Value, na.rm = TRUE),
      y_max = max(Value, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      y_range = y_max - y_min,
      y_range = if_else(y_range == 0, y_max, y_range),
      y_lower = y_min - y_range * 0.10,
      y_upper = y_max + y_range * 0.30,
      bracket_y = y_max + y_range * 0.11,
      tick_y = y_max + y_range * 0.06,
      label_y = y_max + y_range * 0.16,
      Section = if_else(as.character(ElevationGroup) == "All", "All", "Elevation"),
      Section = factor(Section, levels = c("All", "Elevation")),
      x_center = unname(elevation_positions[as.character(ElevationGroup)]),
      x_left = x_center - dodge_width / 4.6,
      x_right = x_center + dodge_width / 4.6
    ) %>%
    left_join(alpha$stats, by = c("Metric", "ElevationGroup"))

  span_df <- tibble(
    Section = factor(c("All", "All", "Elevation", "Elevation"), levels = c("All", "Elevation")),
    x = c(0.50, 1.50, 1.50, 6.50),
    y = c(y_extent$y_lower, y_extent$y_upper, y_extent$y_lower, y_extent$y_upper)
  )

  ggplot(df, aes(x = ElevationX, y = Value, fill = Group)) +
    geom_blank(data = span_df, aes(x = x, y = y), inherit.aes = FALSE) +
    geom_boxplot(
      aes(group = interaction(ElevationX, Group)),
      position = position_dodge(width = dodge_width),
      width = box_width,
      outlier.shape = 21,
      outlier.size = 1.45,
      outlier.stroke = 0.25,
      outlier.colour = "black",
      colour = "black",
      linewidth = 0.72,
      alpha = 0.90,
      show.legend = FALSE
    ) +
    geom_segment(
      data = ranges,
      aes(x = x_left, xend = x_right, y = bracket_y, yend = bracket_y),
      inherit.aes = FALSE,
      colour = "#9E9E9E",
      linewidth = 0.7
    ) +
    geom_segment(
      data = ranges,
      aes(x = x_left, xend = x_left, y = tick_y, yend = bracket_y),
      inherit.aes = FALSE,
      colour = "#9E9E9E",
      linewidth = 0.7
    ) +
    geom_segment(
      data = ranges,
      aes(x = x_right, xend = x_right, y = tick_y, yend = bracket_y),
      inherit.aes = FALSE,
      colour = "#9E9E9E",
      linewidth = 0.7
    ) +
    geom_text(
      data = ranges,
      aes(x = x_center, y = label_y, label = label, colour = label != "NS"),
      inherit.aes = FALSE,
      family = "Helvetica",
      fontface = "bold",
      size = 3.35
    ) +
    scale_fill_manual(values = group_colours, drop = FALSE) +
    scale_colour_manual(values = c("TRUE" = "red", "FALSE" = "black"), guide = "none") +
    scale_x_continuous(
      breaks = 1:6,
      labels = c("All", "3800", "3900", "4000", "4100", "4200"),
      expand = expansion(mult = c(0, 0))
    ) +
    scale_y_continuous(limits = c(y_extent$y_lower, y_extent$y_upper), expand = expansion(mult = c(0, 0))) +
    labs(x = NULL, y = NULL) +
    facet_grid(. ~ Section, scales = "free_x", space = "free_x") +
    theme_classic(base_size = 12, base_family = "Helvetica") +
    theme(
      panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.9),
      axis.line = element_blank(),
      axis.ticks.x = element_line(colour = "black", linewidth = 0.65),
      axis.ticks.y = element_line(colour = "black", linewidth = 0.75),
      axis.ticks.length = grid::unit(0.12, "cm"),
      axis.text.x = element_text(colour = "black", face = "bold", size = 7.6, angle = 0),
      axis.text.y = element_text(
        colour = "black",
        size = 8.7,
        angle = 90,
        hjust = 0.5,
        vjust = 0.5
      ),
      axis.title = element_blank(),
      plot.margin = margin(0, 2, 0, 18),
      legend.position = "none",
      strip.background = element_blank(),
      strip.text = element_blank(),
      panel.spacing.x = grid::unit(0.30, "cm")
    )
}

make_alpha_figure <- function(alpha, title = "Functional richness/diversity") {
  richness_panel <- wrap_plots(
    make_title_strip("Functional gene richness", text_size = 4.1),
    make_alpha_body(alpha, "Richness"),
    ncol = 1,
    heights = c(0.20, 1)
  )
  shannon_panel <- wrap_plots(
    make_title_strip("Shannon index", text_size = 4.1),
    make_alpha_body(alpha, "Shannon"),
    ncol = 1,
    heights = c(0.20, 1)
  )

  wrap_plots(richness_panel, shannon_panel, ncol = 1, heights = c(1, 1))
}

bar_theme <- function(base_size = 11) {
  theme_classic(base_size = base_size, base_family = "Helvetica") +
    theme(
      panel.grid.major.x = element_line(colour = "#B8B8B8", linewidth = 0.45),
      panel.grid.major.y = element_blank(),
      panel.grid.minor = element_blank(),
      axis.line = element_line(colour = "black", linewidth = 0.8),
      axis.ticks = element_line(colour = "black", linewidth = 0.7),
      axis.text = element_text(colour = "black"),
      axis.text.y = element_text(face = "bold", size = base_size + 1),
      axis.title.x = element_text(face = "bold", size = base_size + 1),
      axis.title.y = element_blank(),
      panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.9),
      legend.title = element_text(face = "bold", size = base_size),
      legend.text = element_text(size = base_size - 1),
      legend.key = element_rect(colour = NA),
      plot.margin = margin(5, 5, 5, 5)
    )
}

sample_sort_key <- function(sample) {
  group_rank <- if_else(grepl("^S", sample), 1L, 2L)
  altitude <- as.integer(sub("^[SP]([0-9]+)_.*$", "\\1", sample))
  replicate <- as.integer(sub("^[SP][0-9]+_([0-9]+)$", "\\1", sample))
  tibble(group_rank = group_rank, altitude = altitude, replicate = replicate)
}

mag_sample_plot_theme <- function(base_size = 13) {
  theme_classic(base_size = base_size, base_family = "Helvetica") +
    theme(
      panel.grid.major.y = element_line(colour = "#B8B8B8", linewidth = 0.55),
      panel.grid.major.x = element_blank(),
      panel.grid.minor = element_blank(),
      axis.line = element_line(colour = "black", linewidth = 0.9),
      axis.ticks = element_line(colour = "black", linewidth = 0.8),
      axis.ticks.length = grid::unit(0.16, "cm"),
      axis.text = element_text(colour = "black"),
      axis.text.x = element_text(face = "bold", size = base_size - 1.5, angle = 45, hjust = 1),
      axis.text.y = element_text(size = base_size - 1),
      axis.title.y = element_text(face = "bold", size = base_size + 1),
      axis.title.x = element_blank(),
      plot.title = element_text(face = "bold", size = base_size + 2.5, hjust = 0.5),
      panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.9),
      plot.margin = margin(7, 4, 7, 7)
    )
}

prepare_mag_sample_abundance <- function(mag_long) {
  samples <- mag_long %>%
    distinct(sample, Group) %>%
    bind_cols(sample_sort_key(.$sample)) %>%
    arrange(group_rank, altitude, replicate) %>%
    mutate(SampleIndex = row_number()) %>%
    select(sample, Group, Altitude = altitude, Replicate = replicate, SampleIndex)

  if (nrow(samples) != 60 ||
      sum(samples$Group == "R") != 30 ||
      sum(samples$Group == "E") != 30) {
    stop("MAG table must contain 30 R samples and 30 E samples named like S38_1/P38_1.")
  }

  plot_values <- mag_long %>%
    group_by(sample, Phylum_clean) %>%
    summarise(Abundance = sum(Abundance, na.rm = TRUE), .groups = "drop") %>%
    right_join(
      crossing(samples, Phylum_clean = phylum_order),
      by = c("sample", "Phylum_clean")
    ) %>%
    mutate(
      Abundance = replace_na(Abundance, 0),
      Phylum_clean = factor(Phylum_clean, levels = phylum_order)
    ) %>%
    arrange(SampleIndex, Phylum_clean) %>%
    group_by(SampleIndex) %>%
    mutate(
      ymin = lag(cumsum(Abundance), default = 0),
      ymax = cumsum(Abundance),
      xmin = SampleIndex - 0.5,
      xmax = SampleIndex + 0.5
    ) %>%
    ungroup()

  sample_totals <- plot_values %>%
    group_by(sample, SampleIndex) %>%
    summarise(Total = max(ymax, na.rm = TRUE), .groups = "drop")

  max_total <- max(sample_totals$Total, na.rm = TRUE)
  y_breaks <- pretty(c(0, max_total), n = 5)
  y_breaks <- y_breaks[y_breaks >= 0]
  y_max <- max(y_breaks)
  if (y_max < max_total) {
    y_max <- ceiling(max_total)
    y_breaks <- pretty(c(0, y_max), n = 5)
    y_breaks <- y_breaks[y_breaks >= 0 & y_breaks <= y_max]
  }
  strip_h <- y_max * 0.14
  strip_top <- y_max + strip_h

  list(
    plot_values = plot_values,
    sample_totals = sample_totals,
    sample_count = nrow(samples),
    y_breaks = y_breaks,
    y_max = y_max,
    strip_top = strip_top,
    altitude_centres = samples %>%
      group_by(Group, Altitude) %>%
      summarise(Centre = mean(SampleIndex), .groups = "drop") %>%
      arrange(factor(Group, levels = c("R", "E")), Altitude) %>%
      mutate(Label = as.character(Altitude * 100)),
    altitude_lines = tibble(
      x = c(6.5, 12.5, 18.5, 24.5, 36.5, 42.5, 48.5, 54.5),
      y = 0,
      yend = y_max
    ),
    top_strip = tibble(
      xmin = c(0.5, 30.5),
      xmax = c(30.5, 60.5),
      ymin = c(y_max, y_max),
      ymax = c(strip_top, strip_top),
      Label = c("Rhizomatousphere", "Endophytic"),
      x = c(15.5, 45.5),
      y = y_max + strip_h * 0.5
    )
  )
}

make_mag_sample_panel <- function(prepared, show_legend = TRUE, show_y_title = TRUE) {
  ggplot(prepared$plot_values, aes(fill = Phylum_clean)) +
    geom_rect(
      data = prepared$top_strip,
      aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
      inherit.aes = FALSE,
      fill = "#E6E6E6",
      colour = "black",
      linewidth = 0.75
    ) +
    geom_text(
      data = prepared$top_strip,
      aes(x = x, y = y, label = Label),
      inherit.aes = FALSE,
      family = "Helvetica",
      fontface = "bold",
      size = 4.1
    ) +
    geom_rect(
      aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
      colour = NA,
      linewidth = 0
    ) +
    geom_segment(
      data = prepared$altitude_lines,
      aes(x = x, xend = x, y = y, yend = yend),
      inherit.aes = FALSE,
      colour = "#C8C8C8",
      linewidth = 0.32
    ) +
    annotate(
      "segment",
      x = 30.5,
      xend = 30.5,
      y = 0,
      yend = prepared$strip_top,
      colour = "black",
      linewidth = 0.55
    ) +
    scale_fill_manual(values = phylum_colours, breaks = phylum_order, drop = FALSE) +
    scale_x_continuous(
      breaks = prepared$altitude_centres$Centre,
      labels = prepared$altitude_centres$Label,
      expand = expansion(mult = c(0, 0), add = c(0, 0))
    ) +
    scale_y_continuous(
      breaks = prepared$y_breaks,
      expand = expansion(mult = c(0, 0), add = c(0, 0))
    ) +
    coord_cartesian(
      xlim = c(0.5, prepared$sample_count + 0.5),
      ylim = c(0, prepared$strip_top),
      clip = "on"
    ) +
    labs(
      title = "MAGs",
      x = NULL,
      y = if (show_y_title) "Relative abundance (%)" else NULL,
      fill = "Phylum"
    ) +
    guides(
      fill = guide_legend(
        title.position = "top",
        title.hjust = 0,
        ncol = 1,
        byrow = TRUE,
        keyheight = grid::unit(0.30, "cm"),
        keywidth = grid::unit(0.48, "cm")
      )
    ) +
    mag_sample_plot_theme() +
    theme(
      legend.position = if (show_legend) "right" else "none",
      legend.title = element_text(face = "bold", size = 12),
      legend.text = element_text(size = 10.5),
      legend.spacing.y = grid::unit(0.015, "cm"),
      legend.key = element_rect(colour = NA),
      legend.margin = margin(0, 0, 0, 0)
    )
}

# -----------------------------------------------------------------------------
# MAG table: phylum abundance and MAG-level PCoA
# -----------------------------------------------------------------------------
mags <- read.delim(mags_file, check.names = FALSE)
mag_sample_cols <- names(mags)[grepl("^[PS][0-9]{2}_[0-9]+$", names(mags))]
mag_meta <- sample_info(mag_sample_cols)

mag_long <- mags %>%
  mutate(Phylum_clean = clean_phylum(Phylum)) %>%
  select(Genome, Phylum_clean, all_of(mag_sample_cols)) %>%
  pivot_longer(all_of(mag_sample_cols), names_to = "sample", values_to = "Abundance") %>%
  left_join(mag_meta, by = "sample") %>%
  mutate(Abundance = as.numeric(Abundance))

mag_sample_abundance <- prepare_mag_sample_abundance(mag_long)
write.csv(
  mag_sample_abundance$plot_values %>%
    select(sample, Group, Altitude, Replicate, SampleIndex, Phylum_clean, Abundance, ymin, ymax),
  file.path(out_dir, "MAG_phylum_abundance_60_samples.csv"),
  row.names = FALSE
)

p_mag_phylum <- make_mag_sample_panel(mag_sample_abundance, show_legend = TRUE)
p_mag_phylum_nolegend <- make_mag_sample_panel(mag_sample_abundance, show_legend = FALSE) +
  theme(plot.margin = margin(7, 4, 7, 24))

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
write.csv(as.data.frame(mag_adonis), file.path(out_dir, "MAG_PCoA_adonis.csv"))

p_mag_pcoa <- make_ordination_plot(
  df = mag_scores,
  x_col = "PCoA1",
  y_col = "PCoA2",
  pc = mag_var,
  adonis_stats = mag_adonis_stats,
  title = "MAGs",
  x_prefix = "PCoA1",
  y_prefix = "PCoA2",
  show_legend = TRUE,
  base_size = 11
)
p_mag_pcoa_nolegend <- make_ordination_plot(
  df = mag_scores,
  x_col = "PCoA1",
  y_col = "PCoA2",
  pc = mag_var,
  adonis_stats = mag_adonis_stats,
  title = "MAGs",
  x_prefix = "PCoA1",
  y_prefix = "PCoA2",
  show_legend = FALSE,
  base_size = 10
)

# -----------------------------------------------------------------------------
# Functional genes: alpha diversity and PCA
# -----------------------------------------------------------------------------
div_df <- read.delim(div_file, check.names = FALSE) %>%
  rename(sample = `feature-id`) %>%
  left_join(sample_info(.$sample), by = "sample") %>%
  mutate(Group = factor(Group, levels = c("R", "E")))

alpha_long <- div_df %>%
  select(sample, PairID, Group, Elevation, Replicate, Richness, Shannon) %>%
  pivot_longer(c(Richness, Shannon), names_to = "Metric", values_to = "Value") %>%
  mutate(Metric = factor(Metric, levels = c("Richness", "Shannon")))

pair_check <- alpha_long %>%
  distinct(PairID, Group) %>%
  count(PairID, Group) %>%
  pivot_wider(names_from = Group, values_from = n, values_fill = 0)
if (nrow(pair_check) != 30 ||
    !all(c("R", "E") %in% names(pair_check)) ||
    any(pair_check$R != 1) ||
    any(pair_check$E != 1)) {
  stop("div.txt must contain 30 complete R/E pairs for paired Wilcoxon.")
}

alpha_stats <- alpha_long %>%
  mutate(ElevationGroup = "All") %>%
  bind_rows(
    alpha_long %>%
      mutate(ElevationGroup = paste0(Elevation, " m"))
  ) %>%
  mutate(ElevationGroup = factor(ElevationGroup, levels = c("All", paste0(c(3800, 3900, 4000, 4100, 4200), " m")))) %>%
  select(PairID, Group, Metric, ElevationGroup, Value) %>%
  pivot_wider(names_from = Group, values_from = Value) %>%
  group_by(Metric, ElevationGroup) %>%
  summarise(
    p_value = wilcox.test(R, E, paired = TRUE, exact = FALSE)$p.value,
    .groups = "drop"
  ) %>%
  group_by(Metric) %>%
  mutate(
    q_value = p.adjust(p_value, method = "BH"),
    label = sig_symbol(q_value)
  ) %>%
  ungroup()
write.csv(alpha_stats, file.path(out_dir, "functional_gene_alpha_elevation_paired_wilcoxon_BH.csv"), row.names = FALSE)

alpha_plot_long <- alpha_long %>%
  mutate(ElevationGroup = "All") %>%
  bind_rows(
    alpha_long %>%
      mutate(ElevationGroup = paste0(Elevation, " m"))
  ) %>%
  mutate(ElevationGroup = factor(ElevationGroup, levels = c("All", paste0(c(3800, 3900, 4000, 4100, 4200), " m"))))

alpha_obj <- list(data = alpha_plot_long, stats = alpha_stats)
p_alpha <- make_alpha_figure(alpha_obj, "Functional gene richness and diversity")

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
tpm_sample_mat <- log1p(tpm_sample_mat)
keep <- apply(tpm_sample_mat, 2, sd, na.rm = TRUE) > 0
tpm_sample_mat <- tpm_sample_mat[, keep, drop = FALSE]

func_pca <- prcomp(tpm_sample_mat, center = TRUE, scale. = FALSE)
func_var <- (func_pca$sdev^2 / sum(func_pca$sdev^2)) * 100
func_scores <- as.data.frame(func_pca$x[, 1:2, drop = FALSE])
func_scores$sample <- rownames(func_scores)
func_scores <- left_join(func_scores, tpm_meta, by = "sample")

func_dist <- dist(tpm_sample_mat, method = "euclidean")
set.seed(20260625)
func_adonis <- adonis2(func_dist ~ Group, data = tpm_meta, permutations = 999)
func_adonis_stats <- tibble(R2 = func_adonis$R2[1], p_value = func_adonis$`Pr(>F)`[1])
write.csv(as.data.frame(func_adonis), file.path(out_dir, "functional_gene_PCA_adonis.csv"))

p_func_pca <- make_ordination_plot(
  df = func_scores,
  x_col = "PC1",
  y_col = "PC2",
  pc = func_var,
  adonis_stats = func_adonis_stats,
  title = "Functional genes",
  x_prefix = "PC1",
  y_prefix = "PC2",
  show_legend = TRUE,
  base_size = 11
)
p_func_pca_nolegend <- make_ordination_plot(
  df = func_scores,
  x_col = "PC1",
  y_col = "PC2",
  pc = func_var,
  adonis_stats = func_adonis_stats,
  title = "Functional genes",
  x_prefix = "PC1",
  y_prefix = "PC2",
  show_legend = FALSE,
  base_size = 10
)

save_pair(p_mag_phylum, "01_MAG_phylum_abundance_R_E", width = 8.6, height = 4.2)
save_pair(p_mag_pcoa, "02_MAG_level_PCoA_R_E", width = 3.9, height = 3.35)
save_pair(p_alpha, "03_functional_gene_alpha_diversity_R_E", width = 5.2, height = 4.15)
save_pair(p_func_pca, "04_functional_gene_PCA_R_E", width = 3.9, height = 3.35)

combined_right <- (wrap_elements(full = p_mag_phylum_nolegend) + wrap_elements(full = p_alpha)) /
  (wrap_elements(full = p_mag_pcoa_nolegend) + wrap_elements(full = p_func_pca_nolegend)) +
  plot_layout(widths = c(1.28, 1), heights = c(0.78, 1.05)) +
  plot_annotation(tag_levels = list(c("(b)", "(d)", "(c)", "(e)"))) &
  theme(plot.tag = element_text(face = "bold", size = 16, family = "Helvetica"))

save_pair(combined_right, "metagenome_right_four_panels_combined", width = 9.2, height = 6.2)

vertical_right <- wrap_elements(
  full = p_mag_phylum_nolegend +
    theme(plot.margin = margin(7, 4, 7, 28))
) /
  wrap_elements(
    full = p_mag_pcoa_nolegend +
      theme(aspect.ratio = NULL, plot.margin = margin(7, 4, 7, 28))
  ) /
  wrap_elements(
    full = p_alpha +
      theme(plot.margin = margin(7, 4, 7, 0))
  ) /
  wrap_elements(
    full = p_func_pca_nolegend +
      theme(aspect.ratio = NULL, plot.margin = margin(7, 4, 7, 28))
  ) +
  plot_layout(ncol = 1, heights = c(1.00, 1.08, 1.00, 1.08)) +
  plot_annotation(tag_levels = list(c("(b)", "(c)", "(d)", "(e)"))) &
  theme(plot.tag = element_text(face = "bold", size = 16, family = "Helvetica"))

save_pair(vertical_right, "metagenome_right_four_panels_vertical", width = 5.2, height = 12.8)

message("Done. Outputs written to: ", out_dir)
