#!/usr/bin/env Rscript

# Draw only the functional-gene alpha diversity figure.
# Output size follows the supplied reference image:
# reference chart/image size = 715 x 716 px -> 7.15 x 7.16 inch at 100 dpi.

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(patchwork)
  library(png)
})

base_dir <- "/Volumes/Haleo_DK01/WD_Takakia_F/metagenome-data"
out_dir <- file.path(base_dir, "metagenome_redraw")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

div_file <- file.path(base_dir, "div.txt")
ref_image <- "/Users/haleoriu/Library/Containers/com.tencent.qq/Data/Library/Application Support/QQ/nt_qq_4bae52e983238d99afb091a3553917d4/nt_data/Pic/2026-07/Ori/cb2d028117c7a986b8c7247b8e06e3e0.png"

ref <- png::readPNG(ref_image)
REF_WIDTH_PX <- dim(ref)[2]
REF_HEIGHT_PX <- dim(ref)[1]
OUT_DPI <- 100
OUT_WIDTH <- REF_WIDTH_PX / OUT_DPI
OUT_HEIGHT <- REF_HEIGHT_PX / OUT_DPI

group_colours <- c("R" = "#82B3D7", "E" = "#CF6C72")
group_labels <- c("R" = "Rhizomatousphere", "E" = "Endophytic")

pt_to_lw <- function(pt) pt / 2.845276
lw_box <- pt_to_lw(0.5)
lw_axis <- pt_to_lw(0.8)
lw_strip <- pt_to_lw(0.9)
lw_sig <- pt_to_lw(0.8)

sample_info <- function(samples) {
  tibble(sample = samples) %>%
    filter(grepl("^[PS][0-9]{2}_[0-9]+$", sample)) %>%
    mutate(
      Prefix = substr(sample, 1, 1),
      Group = if_else(Prefix == "S", "R", "E"),
      PairID = sub("^[PS]", "", sample),
      Elevation = as.integer(substr(sample, 2, 3)) * 100,
      Replicate = as.integer(sub("^[PS][0-9]{2}_", "", sample)),
      Group = factor(Group, levels = c("R", "E"))
    )
}

sig_symbol <- function(q) {
  case_when(
    q < 0.001 ~ "***",
    q < 0.01 ~ "**",
    q < 0.05 ~ "*",
    TRUE ~ "NS"
  )
}

div_df <- read.delim(div_file, check.names = FALSE) %>%
  rename(sample = `feature-id`) %>%
  left_join(sample_info(.$sample), by = "sample") %>%
  mutate(Group = factor(Group, levels = c("R", "E")))

alpha_long <- div_df %>%
  select(sample, PairID, Group, Elevation, Replicate, Richness, Shannon) %>%
  pivot_longer(c(Richness, Shannon), names_to = "Metric", values_to = "Value") %>%
  mutate(Metric = factor(Metric, levels = c("Richness", "Shannon")))

plot_levels <- c("All samples", "3800", "3900", "4000", "4100", "4200")

alpha_plot <- bind_rows(
  alpha_long %>%
    mutate(XGroup = "All samples"),
  alpha_long %>%
    mutate(XGroup = as.character(Elevation))
) %>%
  mutate(
    XGroup = factor(XGroup, levels = plot_levels),
    X = as.numeric(XGroup)
  )

alpha_stats <- alpha_plot %>%
  select(PairID, Group, Metric, XGroup, Value) %>%
  pivot_wider(names_from = Group, values_from = Value) %>%
  group_by(Metric, XGroup) %>%
  summarise(
    p_value = wilcox.test(R, E, paired = TRUE, exact = FALSE)$p.value,
    .groups = "drop"
  ) %>%
  group_by(Metric) %>%
  mutate(
    q_value = p.adjust(p_value, method = "BH"),
    label = sig_symbol(q_value),
    X = as.numeric(factor(XGroup, levels = plot_levels))
  ) %>%
  ungroup()

write.csv(
  alpha_stats,
  file.path(out_dir, "functional_gene_alpha_reference_size_stats.csv"),
  row.names = FALSE
)

make_title_strip <- function(title, text_size = 8.3) {
  ggplot() +
    annotate("rect", xmin = 0, xmax = 1, ymin = 0, ymax = 1,
             fill = "#E6E6E6", colour = "black", linewidth = lw_strip) +
    annotate("text", x = 0.5, y = 0.5, label = title,
             family = "Helvetica", fontface = "bold", size = text_size) +
    coord_cartesian(xlim = c(0, 1), ylim = c(0, 1), expand = FALSE) +
    theme_void() +
    theme(plot.margin = margin(0, 0, 0, 0))
}

make_alpha_panel <- function(metric_name) {
  df <- alpha_plot %>%
    filter(as.character(Metric) == metric_name)

  y_extent <- df %>%
    summarise(y_min = min(Value, na.rm = TRUE), y_max = max(Value, na.rm = TRUE)) %>%
    mutate(
      y_range = y_max - y_min,
      y_lower = y_min - y_range * 0.14,
      y_upper = y_max + y_range * 0.24
    )

  stat_df <- df %>%
    group_by(Metric, XGroup, X) %>%
    summarise(
      y_min = min(Value, na.rm = TRUE),
      y_max = max(Value, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      y_range = y_max - y_min,
      y_range = if_else(y_range == 0, y_extent$y_range * 0.08, y_range),
      bracket_y = y_max + y_range * 0.16,
      tick_y = y_max + y_range * 0.10,
      label_y = y_max + y_range * 0.21,
      x_left = X - 0.17,
      x_right = X + 0.17
    ) %>%
    left_join(
      alpha_stats %>% select(Metric, XGroup, label),
      by = c("Metric", "XGroup")
    )

  ggplot(df, aes(x = X, y = Value, fill = Group)) +
    geom_boxplot(
      aes(group = interaction(X, Group)),
      position = position_dodge(width = 0.60),
      width = 0.52,
      linewidth = lw_box,
      colour = "black",
      outlier.shape = 21,
      outlier.size = 0.65,
      outlier.stroke = 0.12,
      outlier.colour = "black",
      alpha = 0.92,
      show.legend = metric_name == "Richness"
    ) +
    geom_segment(
      data = stat_df,
      aes(x = x_left, xend = x_right, y = bracket_y, yend = bracket_y),
      inherit.aes = FALSE,
      colour = "#9E9E9E",
      linewidth = lw_sig
    ) +
    geom_segment(
      data = stat_df,
      aes(x = x_left, xend = x_left, y = tick_y, yend = bracket_y),
      inherit.aes = FALSE,
      colour = "#9E9E9E",
      linewidth = lw_sig
    ) +
    geom_segment(
      data = stat_df,
      aes(x = x_right, xend = x_right, y = tick_y, yend = bracket_y),
      inherit.aes = FALSE,
      colour = "#9E9E9E",
      linewidth = lw_sig
    ) +
    geom_text(
      data = stat_df,
      aes(x = X, y = label_y, label = label, colour = label != "NS"),
      inherit.aes = FALSE,
      family = "Helvetica",
      fontface = "bold",
      size = 4.1
    ) +
    scale_fill_manual(values = group_colours, labels = group_labels, drop = FALSE) +
    scale_colour_manual(values = c("TRUE" = "red", "FALSE" = "black"), guide = "none") +
    scale_x_continuous(
      breaks = seq_along(plot_levels),
      labels = plot_levels,
      limits = c(0.50, length(plot_levels) + 0.50),
      expand = expansion(mult = c(0, 0))
    ) +
    scale_y_continuous(
      limits = c(y_extent$y_lower, y_extent$y_upper),
      expand = expansion(mult = c(0, 0))
    ) +
    labs(x = NULL, y = NULL, fill = NULL) +
    guides(fill = guide_legend(nrow = 1, byrow = TRUE, override.aes = list(size = 5))) +
    theme_classic(base_size = 13, base_family = "Helvetica") +
    theme(
      panel.border = element_rect(colour = "black", fill = NA, linewidth = lw_axis),
      axis.line = element_blank(),
      axis.ticks = element_line(colour = "black", linewidth = lw_axis),
      axis.ticks.length = grid::unit(0.11, "cm"),
      axis.text.x = element_text(colour = "black", size = 12),
      axis.text.y = element_text(colour = "black", size = 12, angle = 90, hjust = 0.5, vjust = 0.5),
      legend.position = if (metric_name == "Richness") "bottom" else "none",
      legend.text = element_text(size = 11),
      legend.key.size = grid::unit(0.42, "cm"),
      legend.margin = margin(-2, 0, -1, 0),
      legend.spacing.x = grid::unit(0.16, "cm"),
      plot.margin = margin(0, 0, 0, 0)
    )
}

p_richness <- wrap_plots(
  make_title_strip("Functional gene richness"),
  make_alpha_panel("Richness"),
  ncol = 1,
  heights = c(0.18, 1)
)

p_shannon <- wrap_plots(
  make_title_strip("Functional gene shannon index"),
  make_alpha_panel("Shannon"),
  ncol = 1,
  heights = c(0.18, 1)
)

p_final <- wrap_plots(
  p_richness,
  p_shannon,
  ncol = 1,
  heights = c(1, 1)
) &
  theme(plot.margin = margin(0, 0, 0, 0))

p_final <- p_final +
  plot_annotation(
    theme = theme(
      plot.margin = margin(0, 0, 0, 0)
    )
  ) &
  theme(
    plot.background = element_rect(fill = "white", colour = NA)
  )

out_stem <- "functional_gene_alpha_reference_size_only"
ggsave(
  file.path(out_dir, paste0(out_stem, ".png")),
  p_final,
  width = OUT_WIDTH,
  height = OUT_HEIGHT,
  units = "in",
  dpi = OUT_DPI,
  bg = "white"
)
ggsave(
  file.path(out_dir, paste0(out_stem, ".pdf")),
  p_final,
  width = OUT_WIDTH,
  height = OUT_HEIGHT,
  units = "in",
  device = grDevices::pdf,
  bg = "white"
)

message(
  "Done. Reference image = ", REF_WIDTH_PX, " x ", REF_HEIGHT_PX,
  " px; output = ", sprintf("%.2f", OUT_WIDTH), " x ",
  sprintf("%.2f", OUT_HEIGHT), " inches at ", OUT_DPI, " dpi."
)
