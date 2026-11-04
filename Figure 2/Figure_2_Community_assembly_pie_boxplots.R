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
  normalizePath("assembly_process_pie_boxplot.R")
}
data_dir <- dirname(script_path)

bacteria_file <- file.path(data_dir, "AP_b.txt")
fungi_file <- file.path(data_dir, "AP_f.txt")

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

process_levels_box <- c("HoS", "DL", "DR")
process_levels_pie <- c("HeS", "HoS", "DL", "HD", "DR")

process_colours <- c(
  "HeS" = "#8FBBD9",
  "HoS" = "#A8CFA0",
  "DL" = "#E7C07A",
  "HD" = "#B9A7CF",
  "DR" = "#D98F8F"
)

pie_values <- tibble(
  Kingdom = rep(c("Bacteria", "Fungi"), each = 10),
  Group = rep(rep(c("R", "E"), each = 5), times = 2),
  Process = rep(process_levels_pie, times = 4),
  Percentage = c(
    0.42571648, 48.96315710, 27.94114940, 1.74018391, 20.92979310,
    2.37826820, 39.71314940, 10.69718010, 0.61995402, 46.59144830,
    2.29996853, 11.80347580, 52.31548120, 0.78045558, 32.80061890,
    0.16582212, 8.39516603, 53.97217120, 0.75719084, 36.70964980
  )
) %>%
  mutate(
    Group = factor(Group, levels = c("R", "E")),
    Process = factor(Process, levels = process_levels_pie)
  )

sig_symbol <- function(q) {
  case_when(
    q < 0.001 ~ "***",
    q < 0.01 ~ "**",
    q < 0.05 ~ "*",
    TRUE ~ "NS"
  )
}

read_process <- function(path) {
  raw <- read.delim(
    path,
    header = TRUE,
    sep = "\t",
    check.names = FALSE,
    stringsAsFactors = FALSE
  )

  required_cols <- as.vector(outer(process_levels_box, c("R", "E"), paste, sep = "_"))
  if (!all(required_cols %in% colnames(raw))) {
    stop(
      basename(path),
      " must contain columns: ",
      paste(required_cols, collapse = ", ")
    )
  }

  long <- raw %>%
    mutate(RowID = row_number()) %>%
    pivot_longer(
      cols = all_of(required_cols),
      names_to = c("Process", "Group"),
      names_sep = "_",
      values_to = "Value"
    ) %>%
    mutate(
      Value = as.numeric(Value),
      Group = factor(Group, levels = c("R", "E")),
      Process = factor(Process, levels = process_levels_box)
    )

  stats <- long %>%
    group_by(Process) %>%
    summarise(
      p_value = wilcox.test(
        Value[Group == "R"],
        Value[Group == "E"],
        paired = FALSE,
        exact = FALSE
      )$p.value,
      .groups = "drop"
    ) %>%
    mutate(
      q_value = p.adjust(p_value, method = "BH"),
      label = sig_symbol(q_value)
    )

  list(data = long, stats = stats)
}

make_process_body <- function(process_data) {
  df <- process_data$data %>%
    mutate(
      GroupX = if_else(as.character(Group) == "R", 1, 1.52),
      ProcessLabel = factor(as.character(Process), levels = process_levels_box)
    )

  ranges <- df %>%
    group_by(Process, ProcessLabel) %>%
    summarise(
      y_min = min(Value, na.rm = TRUE),
      y_max = max(Value, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      y_range = y_max - y_min,
      y_range = if_else(y_range == 0, y_max, y_range),
      y_lower = y_min - y_range * 0.06,
      y_upper = y_max + y_range * 0.24,
      bracket_y = y_max + y_range * 0.09,
      tick_y = y_max + y_range * 0.06,
      label_y = y_max + y_range * 0.14
    ) %>%
    left_join(process_data$stats, by = "Process")

  ggplot(df, aes(x = GroupX, y = Value, fill = Group, group = Group)) +
    geom_blank(
      data = ranges,
      aes(x = 1.26, y = y_lower),
      inherit.aes = FALSE
    ) +
    geom_blank(
      data = ranges,
      aes(x = 1.26, y = y_upper),
      inherit.aes = FALSE
    ) +
    geom_boxplot(
      width = 0.28,
      linewidth = 0.75,
      alpha = 0.82,
      outlier.shape = 21,
      outlier.size = 1.2,
      outlier.stroke = 0.25,
      outlier.alpha = 0.65
    ) +
    geom_segment(
      data = ranges,
      aes(x = 1, xend = 1.52, y = bracket_y, yend = bracket_y),
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
      aes(x = 1.52, xend = 1.52, y = tick_y, yend = bracket_y),
      inherit.aes = FALSE,
      colour = "#9E9E9E",
      linewidth = 0.7
    ) +
    geom_text(
      data = ranges,
      aes(x = 1.26, y = label_y, label = label, colour = label != "NS"),
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
      limits = c(0.76, 1.76),
      expand = expansion(mult = c(0, 0))
    ) +
    scale_y_continuous(expand = expansion(mult = c(0, 0))) +
    facet_wrap(
      ~ ProcessLabel,
      nrow = 1,
      scales = "free_y",
      strip.position = "bottom"
    ) +
    labs(x = NULL, y = "Percentage (%)") +
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
      axis.title.y = element_text(face = "bold", size = 12.5, margin = margin(r = 6)),
      axis.title.x = element_blank(),
      strip.placement = "outside",
      strip.background = element_blank(),
      strip.text.x.bottom = element_text(
        face = "bold",
        size = 13,
        margin = margin(t = 4, b = 0)
      ),
      panel.spacing.x = grid::unit(0.10, "cm"),
      plot.margin = margin(0, 0, 0, 18),
      legend.position = "none"
    )
}

make_donut <- function(df, label) {
  ggplot(df, aes(x = 2, y = Percentage, fill = Process)) +
    geom_col(
      width = 1,
      colour = "white",
      linewidth = 0.35
    ) +
    coord_polar(theta = "y", start = pi / 2) +
    xlim(0.65, 2.55) +
    scale_fill_manual(
      values = process_colours,
      breaks = process_levels_pie,
      drop = FALSE
    ) +
    annotate(
      "text",
      x = 0.65,
      y = 0,
      label = label,
      family = "Helvetica",
      fontface = "bold",
      size = 5
    ) +
    theme_void(base_family = "Helvetica") +
    theme(
      plot.margin = margin(0, 0, 0, 0),
      legend.position = "none"
    )
}

make_pie_column <- function(kingdom) {
  r_pie <- make_donut(
    pie_values %>% filter(Kingdom == kingdom, Group == "R"),
    "R"
  )
  e_pie <- make_donut(
    pie_values %>% filter(Kingdom == kingdom, Group == "E"),
    "E"
  )

  plot_grid(
    r_pie,
    e_pie,
    ncol = 1,
    rel_heights = c(1, 1),
    align = "v"
  )
}

make_process_legend <- function() {
  legend_plot <- ggplot(
    pie_values,
    aes(x = Process, y = Percentage, fill = Process)
  ) +
    geom_col() +
    scale_fill_manual(
      values = process_colours,
      breaks = process_levels_pie,
      drop = FALSE
    ) +
    guides(
      fill = guide_legend(
        title = NULL,
        nrow = 1,
        byrow = TRUE,
        keywidth = grid::unit(0.42, "cm"),
        keyheight = grid::unit(0.30, "cm")
      )
    ) +
    theme_void(base_family = "Helvetica") +
    theme(
      legend.position = "bottom",
      legend.text = element_text(size = 10),
      legend.margin = margin(0, 0, 0, 0)
    )

  get_legend(legend_plot)
}

make_title_strip <- function(title) {
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
    theme(plot.margin = margin(0, 0, 0, 0))
}

make_process_figure <- function(process_data, kingdom, title) {
  body <- plot_grid(
    make_pie_column(kingdom),
    make_process_body(process_data),
    nrow = 1,
    rel_widths = c(0.70, 3.30),
    align = "h",
    axis = "tb"
  )

  plot_grid(
    make_title_strip(title),
    body,
    ncol = 1,
    rel_heights = c(0.13, 1),
    align = "v",
    axis = "lr"
  )
}

bacteria <- read_process(bacteria_file)
fungi <- read_process(fungi_file)

bacteria_figure <- make_process_figure(
  bacteria,
  "Bacteria",
  "Bacterial community assembly"
)

fungi_figure <- make_process_figure(
  fungi,
  "Fungi",
  "Fungal community assembly"
)

legend <- make_process_legend()

combined_body <- plot_grid(
  bacteria_figure,
  fungi_figure,
  nrow = 1,
  rel_widths = c(1, 1),
  align = "h",
  axis = "tb"
)

combined_figure <- plot_grid(
  combined_body,
  legend,
  ncol = 1,
  rel_heights = c(1, 0.12)
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
  file.path(data_dir, "bacteria_assembly_process_pie_boxplot_wilcoxon_BH_FDR.csv"),
  row.names = FALSE
)

write.csv(
  fungi$stats,
  file.path(data_dir, "fungi_assembly_process_pie_boxplot_wilcoxon_BH_FDR.csv"),
  row.names = FALSE
)

save_plot_pair("bacteria_assembly_process_pie_boxplot", bacteria_figure, 7.6, 3.075)
save_plot_pair("fungi_assembly_process_pie_boxplot", fungi_figure, 7.6, 3.075)
save_plot_pair("bacteria_fungi_assembly_process_pie_boxplot_combined", combined_figure, 15.2, 3.45)

message("Created assembly process pie plus boxplot figures.")
