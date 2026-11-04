#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(cowplot)
})

args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
script_path <- if (length(file_arg)) {
  normalizePath(sub("^--file=", "", file_arg[1]))
} else {
  normalizePath("alpha_diversity_violin.R")
}
data_dir <- dirname(script_path)

bacteria_file <- file.path(data_dir, "div_bacteria.txt")
fungi_file <- file.path(data_dir, "div_fungi.txt")

required_files <- c(bacteria_file, fungi_file)
if (!all(file.exists(required_files))) {
  stop(
    "Missing input file(s): ",
    paste(required_files[!file.exists(required_files)], collapse = ", ")
  )
}

group_colours <- c(
  "R" = "#82B3D7",
  "E" = "#CF6C72"
)

metric_labels <- c(
  "Richness" = "Richness",
  "Shannon" = "Shannon index"
)

metric_bottom_labels <- c(
  "Richness" = "Richness",
  "Shannon" = "Shannon"
)

sig_symbol <- function(q) {
  case_when(
    q < 0.001 ~ "***",
    q < 0.01 ~ "**",
    q < 0.05 ~ "*",
    TRUE ~ "NS"
  )
}

read_alpha <- function(path) {
  raw <- read.delim(
    path,
    header = TRUE,
    sep = "\t",
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  names(raw) <- tolower(names(raw))
  if (!all(c("sample", "sobs", "shannon") %in% names(raw))) {
    stop(basename(path), " must contain sample, sobs, and shannon columns.")
  }

  long <- raw %>%
    transmute(
      Sample = sample,
      Group = case_when(
        grepl("^S", Sample) ~ "R",
        grepl("^P", Sample) ~ "E",
        TRUE ~ NA_character_
      ),
      PairID = sub("^[SP]", "", Sample),
      Richness = as.numeric(sobs),
      Shannon = as.numeric(shannon)
    ) %>%
    filter(!is.na(Group)) %>%
    pivot_longer(
      cols = c(Richness, Shannon),
      names_to = "Metric",
      values_to = "Value"
    ) %>%
    mutate(
      Group = factor(Group, levels = c("R", "E")),
      Metric = factor(Metric, levels = c("Richness", "Shannon"))
    )

  pair_check <- long %>%
    distinct(PairID, Group) %>%
    count(PairID, Group) %>%
    pivot_wider(names_from = Group, values_from = n, values_fill = 0)

  if (nrow(pair_check) != 30 ||
      !all(c("R", "E") %in% names(pair_check)) ||
      any(pair_check$R != 1) ||
      any(pair_check$E != 1)) {
    stop(basename(path), " must contain 30 complete R/E pairs.")
  }

  stats <- long %>%
    select(PairID, Group, Metric, Value) %>%
    pivot_wider(names_from = Group, values_from = Value) %>%
    group_by(Metric) %>%
    summarise(
      p_value = wilcox.test(R, E, paired = TRUE, exact = FALSE)$p.value,
      .groups = "drop"
    ) %>%
    mutate(
      q_value = p.adjust(p_value, method = "BH"),
      label = sig_symbol(q_value)
    )

  list(data = long, stats = stats)
}

panel_theme <- function() {
  theme_classic(base_size = 12, base_family = "Helvetica") +
    theme(
      panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.9),
      axis.line = element_blank(),
      axis.ticks = element_line(colour = "black", linewidth = 0.75),
      axis.ticks.length = grid::unit(0.14, "cm"),
      axis.text = element_text(colour = "black", size = 10.5),
      axis.text.x = element_text(face = "bold", size = 12),
      axis.title.y = element_text(face = "bold", size = 12.5),
      axis.title.x = element_blank(),
      plot.margin = margin(0, 2, 5, 2),
      legend.position = "none"
    )
}

make_alpha_body <- function(alpha) {
  df <- alpha$data %>%
    mutate(
      GroupX = if_else(as.character(Group) == "R", 1, 1.58),
      MetricLabel = recode(as.character(Metric), !!!metric_bottom_labels),
      MetricLabel = factor(MetricLabel, levels = unname(metric_bottom_labels))
    )

  ranges <- df %>%
    group_by(Metric, MetricLabel) %>%
    summarise(
      y_min = min(Value, na.rm = TRUE),
      y_max = max(Value, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      y_range = y_max - y_min,
      y_range = if_else(y_range == 0, y_max, y_range),
      y_lower = y_min - y_range * 0.06,
      y_upper = y_max + y_range * 0.23,
      bracket_y = y_max + y_range * 0.09,
      tick_y = y_max + y_range * 0.06,
      label_y = y_max + y_range * 0.13
    ) %>%
    left_join(alpha$stats, by = "Metric")

  ggplot(df, aes(x = GroupX, y = Value, fill = Group, group = Group)) +
    geom_blank(
      data = ranges,
      aes(x = 1.29, y = y_lower),
      inherit.aes = FALSE
    ) +
    geom_blank(
      data = ranges,
      aes(x = 1.29, y = y_upper),
      inherit.aes = FALSE
    ) +
    geom_violin(
      width = 0.35,
      trim = TRUE,
      alpha = 0.86,
      colour = "black",
      linewidth = 0.85,
      scale = "width"
    ) +
    geom_point(
      position = position_jitter(width = 0.045, height = 0, seed = 20260625),
      shape = 21,
      size = 1.65,
      stroke = 0.35,
      colour = "#3A3A3A",
      alpha = 0.82,
      show.legend = FALSE
    ) +
    geom_segment(
      data = ranges,
      aes(x = 1, xend = 1.58, y = bracket_y, yend = bracket_y),
      inherit.aes = FALSE,
      colour = "#9E9E9E",
      linewidth = 0.7
    ) +
    geom_segment(
      data = ranges,
      aes(x = 1, xend = 1, y = tick_y, yend = bracket_y),
      inherit.aes = FALSE,
      colour = "#9E9E9E",
      linewidth = 0.7
    ) +
    geom_segment(
      data = ranges,
      aes(x = 1.58, xend = 1.58, y = tick_y, yend = bracket_y),
      inherit.aes = FALSE,
      colour = "#9E9E9E",
      linewidth = 0.7
    ) +
    geom_text(
      data = ranges,
      aes(x = 1.29, y = label_y, label = label, colour = label != "NS"),
      inherit.aes = FALSE,
      family = "Helvetica",
      fontface = "bold",
      size = 5
    ) +
    scale_fill_manual(values = group_colours, drop = FALSE) +
    scale_colour_manual(
      values = c("TRUE" = "red", "FALSE" = "black"),
      guide = "none"
    ) +
    scale_x_continuous(
      breaks = NULL,
      limits = c(0.72, 1.86),
      expand = expansion(mult = c(0, 0))
    ) +
    scale_y_continuous(expand = expansion(mult = c(0, 0))) +
    facet_wrap(
      ~ MetricLabel,
      nrow = 1,
      scales = "free_y",
      strip.position = "bottom"
    ) +
    labs(x = NULL, y = NULL) +
    theme_classic(base_size = 12, base_family = "Helvetica") +
    theme(
      panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.9),
      axis.line = element_blank(),
      axis.ticks.x = element_blank(),
      axis.ticks.y = element_line(colour = "black", linewidth = 0.75),
      axis.ticks.length = grid::unit(0.12, "cm"),
      axis.text.x = element_blank(),
      axis.text.y = element_text(
        colour = "black",
        size = 10.5,
        angle = 90,
        hjust = 0.5,
        vjust = 0.5
      ),
      axis.title = element_blank(),
      strip.placement = "outside",
      strip.background = element_blank(),
      strip.text.x.bottom = element_text(
        face = "bold",
        size = 13,
        margin = margin(t = 4, b = 0)
      ),
      panel.spacing.x = grid::unit(0.10, "cm"),
      plot.margin = margin(0, 0, 0, 26),
      legend.position = "none"
    )
}

make_title_strip <- function(title, left_margin = 26) {
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
      size = 5.2
    ) +
    coord_cartesian(xlim = c(0, 1), ylim = c(0, 1), expand = FALSE) +
    theme_void() +
    theme(plot.margin = margin(0, 0, 0, left_margin))
}

make_alpha_figure <- function(alpha, title = "Richness and diversity") {
  body <- make_alpha_body(alpha)

  plot_grid(
    make_title_strip(title),
    body,
    ncol = 1,
    rel_heights = c(0.13, 1),
    align = "v",
    axis = "lr"
  )
}

bacteria <- read_alpha(bacteria_file)
fungi <- read_alpha(fungi_file)

bacteria_figure <- make_alpha_figure(bacteria, "Bacterial richness and diversity")
fungi_figure <- make_alpha_figure(fungi, "Fungal richness and diversity")

combined_figure <- plot_grid(
  bacteria_figure,
  fungi_figure,
  nrow = 1,
  rel_widths = c(1, 1),
  align = "h",
  axis = "tb"
)

save_plot_pair <- function(filename_stem, plot, width, height) {
  ggsave(
    file.path(data_dir, paste0(filename_stem, ".pdf")),
    plot,
    width = width,
    height = height,
    units = "in",
    device = "pdf",
    useDingbats = FALSE,
    bg = "white"
  )
  ggsave(
    file.path(data_dir, paste0(filename_stem, ".png")),
    plot,
    width = width,
    height = height,
    units = "in",
    dpi = 600,
    bg = "white"
  )
}

write.csv(
  bacteria$stats,
  file.path(data_dir, "bacteria_alpha_diversity_paired_wilcoxon_BH_FDR.csv"),
  row.names = FALSE
)

write.csv(
  fungi$stats,
  file.path(data_dir, "fungi_alpha_diversity_paired_wilcoxon_BH_FDR.csv"),
  row.names = FALSE
)

save_plot_pair("bacteria_alpha_diversity_violin", bacteria_figure, 4.4, 3.075)
save_plot_pair("fungi_alpha_diversity_violin", fungi_figure, 4.4, 3.075)
save_plot_pair("bacteria_fungi_alpha_diversity_violin_combined", combined_figure, 8.8, 3.1875)

message("Created alpha diversity violin plots.")
