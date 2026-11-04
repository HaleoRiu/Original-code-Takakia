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
  normalizePath("relative_abundance_60_samples.R")
}
data_dir <- dirname(script_path)

bacteria_file <- file.path(data_dir, "phylum_bacteria.txt")
fungi_file <- file.path(data_dir, "phylum_fungi.txt")

required_files <- c(bacteria_file, fungi_file)
if (!all(file.exists(required_files))) {
  stop(
    "Missing input file(s): ",
    paste(required_files[!file.exists(required_files)], collapse = ", ")
  )
}

bacteria_phyla <- c(
  "Pseudomonadota", "Acidobacteriota", "Actinomycetota",
  "Chloroflexi", "Bacteroidota", "Firmicutes",
  "Verrucomicrobiota", "WPS-2", "Planctomycetota", "Myxococcota"
)

bacteria_colours <- c(
  "Pseudomonadota" = "#FF7F27",
  "Acidobacteriota" = "#9BE9E6",
  "Actinomycetota" = "#F7658A",
  "Chloroflexi" = "#64E85B",
  "Bacteroidota" = "#000080",
  "Firmicutes" = "#FFF125",
  "Verrucomicrobiota" = "#F52525",
  "WPS-2" = "#9EF0CF",
  "Planctomycetota" = "#F43FC2",
  "Myxococcota" = "#C97BE1",
  "Others" = "#4DC4D0"
)

fungi_phyla <- c(
  "Ascomycota", "Basidiomycota", "Mortierellomycota",
  "Unclassified", "Rozellomycota"
)

fungi_colours <- c(
  "Ascomycota" = "#75B5D4",
  "Basidiomycota" = "#A979DD",
  "Mortierellomycota" = "#155A9C",
  "Unclassified" = "#62A94D",
  "Rozellomycota" = "#7770A6",
  "Others" = "#858585"
)

canonicalize_phylum <- function(x) {
  recode(
    x,
    "Proteobacteria" = "Pseudomonadota",
    "Actinobacteriota" = "Actinomycetota",
    .default = x
  )
}

sample_sort_key <- function(sample) {
  group_rank <- if_else(grepl("^S", sample), 1L, 2L)
  altitude <- as.integer(sub("^[SP]([0-9]+)_.*$", "\\1", sample))
  replicate <- as.integer(sub("^[SP][0-9]+_([0-9]+)$", "\\1", sample))
  tibble(group_rank = group_rank, altitude = altitude, replicate = replicate)
}

read_sample_abundance <- function(path, selected_phyla) {
  raw <- read.delim(
    path,
    header = TRUE,
    sep = "\t",
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  names(raw)[1] <- "Phylum"
  raw$Phylum <- canonicalize_phylum(raw$Phylum)

  long <- raw %>%
    pivot_longer(-Phylum, names_to = "Sample", values_to = "Abundance") %>%
    mutate(
      Abundance = as.numeric(Abundance),
      Group = case_when(
        grepl("^S", Sample) ~ "R",
        grepl("^P", Sample) ~ "E",
        TRUE ~ NA_character_
      )
    ) %>%
    filter(!is.na(Group)) %>%
    group_by(Sample, Group, Phylum) %>%
    summarise(Abundance = sum(Abundance), .groups = "drop") %>%
    group_by(Sample) %>%
    mutate(RelativeAbundance = 100 * Abundance / sum(Abundance)) %>%
    ungroup()

  samples <- long %>%
    distinct(Sample, Group) %>%
    bind_cols(sample_sort_key(.$Sample)) %>%
    arrange(group_rank, altitude, replicate) %>%
    mutate(SampleIndex = row_number()) %>%
    select(Sample, Group, Altitude = altitude, SampleIndex)

  if (nrow(samples) != 60 ||
      sum(samples$Group == "R") != 30 ||
      sum(samples$Group == "E") != 30) {
    stop(basename(path), " must contain 30 R samples and 30 E samples.")
  }

  plot_values <- long %>%
    mutate(
      DisplayPhylum = if_else(Phylum %in% selected_phyla, Phylum, "Others")
    ) %>%
    group_by(Sample, Group, DisplayPhylum) %>%
    summarise(RelativeAbundance = sum(RelativeAbundance), .groups = "drop") %>%
    right_join(
      crossing(samples, DisplayPhylum = c(selected_phyla, "Others")),
      by = c("Sample", "Group", "DisplayPhylum")
    ) %>%
    mutate(
      RelativeAbundance = replace_na(RelativeAbundance, 0),
      DisplayPhylum = factor(DisplayPhylum, levels = c(selected_phyla, "Others"))
    ) %>%
    arrange(SampleIndex, DisplayPhylum) %>%
    group_by(SampleIndex) %>%
    mutate(
      ymin = lag(cumsum(RelativeAbundance), default = 0),
      ymax = cumsum(RelativeAbundance),
      xmin = SampleIndex - 0.515,
      xmax = SampleIndex + 0.515
    ) %>%
    ungroup()

  list(
    plot_values = plot_values,
    sample_count = nrow(samples),
    altitude_centres = samples %>%
      group_by(Group, Altitude) %>%
      summarise(Centre = mean(SampleIndex), .groups = "drop") %>%
      arrange(factor(Group, levels = c("R", "E")), Altitude) %>%
      mutate(Label = as.character(Altitude * 100)),
    altitude_gaps = tibble(
      x = c(6.5, 12.5, 18.5, 24.5, 36.5, 42.5, 48.5, 54.5),
      xmin = x - 0.04,
      xmax = x + 0.04,
      ymin = 0,
      ymax = 100
    ),
    top_strip = tibble(
      xmin = c(0.5, 30.5),
      xmax = c(30.5, 60.5),
      ymin = c(100, 100),
      ymax = c(114, 114),
      Label = c("Rhizomatousphere", "Endophytic"),
      x = c(15.5, 45.5),
      y = c(107, 107)
    )
  )
}

sample_plot_theme <- function() {
  theme_classic(base_size = 13, base_family = "Helvetica") +
    theme(
      panel.grid.major.y = element_line(colour = "#B8B8B8", linewidth = 0.55),
      panel.grid.minor = element_blank(),
      axis.line = element_line(colour = "black", linewidth = 0.9),
      axis.ticks = element_line(colour = "black", linewidth = 0.8),
      axis.ticks.length = grid::unit(0.16, "cm"),
      axis.text = element_text(colour = "black"),
      axis.text.x = element_text(face = "bold", size = 11.5, angle = 45, hjust = 1),
      axis.text.y = element_text(size = 12),
      axis.title.y = element_text(face = "bold", size = 14),
      plot.title = element_text(face = "bold", size = 15.5, hjust = 0.5),
      panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.9),
      plot.margin = margin(7, 4, 7, 7)
    )
}

make_sample_panel <- function(
    prepared,
    title,
    colours,
    selected_phyla,
    show_legend = FALSE,
    show_y_title = TRUE) {
  ggplot(
    prepared$plot_values,
    aes(fill = DisplayPhylum)
  ) +
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
    geom_rect(
      data = prepared$altitude_gaps,
      aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
      inherit.aes = FALSE,
      fill = "white",
      colour = NA,
      linewidth = 0
    ) +
    annotate(
      "segment",
      x = 30.5,
      xend = 30.5,
      y = 0,
      yend = 114,
      colour = "black",
      linewidth = 0.55
    ) +
    scale_fill_manual(
      values = colours,
      breaks = c(selected_phyla, "Others"),
      drop = FALSE
    ) +
    scale_x_continuous(
      breaks = prepared$altitude_centres$Centre,
      labels = prepared$altitude_centres$Label,
      expand = expansion(mult = c(0, 0), add = c(0, 0))
    ) +
    scale_y_continuous(
      breaks = seq(0, 100, 20),
      expand = expansion(mult = c(0, 0), add = c(0, 0))
    ) +
    coord_cartesian(
      xlim = c(0.5, prepared$sample_count + 0.5),
      ylim = c(0, 114),
      clip = "on"
    ) +
    labs(
      title = title,
      x = NULL,
      y = if (show_y_title) "Relative abundance (%)" else NULL,
      fill = NULL
    ) +
    guides(
      fill = guide_legend(
        ncol = 1,
        byrow = TRUE,
        keyheight = grid::unit(0.34, "cm"),
        keywidth = grid::unit(0.52, "cm"),
        title.position = "top"
      )
    ) +
    sample_plot_theme() +
    theme(
      legend.position = if (show_legend) "right" else "none",
      legend.title = element_text(face = "bold", size = 12),
      legend.text = element_text(size = 11),
      legend.spacing.y = grid::unit(0.02, "cm"),
      legend.key = element_rect(colour = NA),
      legend.margin = margin(0, 0, 0, 0)
    )
}

make_legend <- function(prepared, title, colours, selected_phyla, heading) {
  legend_plot <- make_sample_panel(
    prepared = prepared,
    title = title,
    colours = colours,
    selected_phyla = selected_phyla,
    show_legend = TRUE
  ) +
    guides(
      fill = guide_legend(
        title = heading,
        title.position = "top",
        title.hjust = 0,
        ncol = 1,
        byrow = TRUE,
        keyheight = grid::unit(0.34, "cm"),
        keywidth = grid::unit(0.52, "cm")
      )
    ) +
    theme(
      legend.title = element_text(face = "bold", size = 12),
      legend.text = element_text(size = 11)
    )
  get_legend(legend_plot)
}

make_single_figure <- function(prepared, title, colours, selected_phyla) {
  panel <- make_sample_panel(
    prepared,
    title,
    colours,
    selected_phyla,
    show_legend = FALSE
  )
  legend <- make_legend(prepared, title, colours, selected_phyla, "Phylum")

  plot_grid(
    panel,
    legend,
    nrow = 1,
    rel_widths = c(1.35, 0.72),
    align = "h",
    axis = "tb"
  )
}

bacteria <- read_sample_abundance(bacteria_file, bacteria_phyla)
fungi <- read_sample_abundance(fungi_file, fungi_phyla)

bacteria_panel <- make_sample_panel(
  bacteria,
  "Bacteria",
  bacteria_colours,
  bacteria_phyla,
  show_legend = FALSE
)

fungi_panel <- make_sample_panel(
  fungi,
  "Fungi",
  fungi_colours,
  fungi_phyla,
  show_legend = FALSE
)

bacteria_figure <- make_single_figure(
  bacteria,
  "Bacteria",
  bacteria_colours,
  bacteria_phyla
)

fungi_figure <- make_single_figure(
  fungi,
  "Fungi",
  fungi_colours,
  fungi_phyla
)

combined_legend <- plot_grid(
  make_legend(bacteria, "Bacteria", bacteria_colours, bacteria_phyla, "Bacteria phyla"),
  make_legend(fungi, "Fungi", fungi_colours, fungi_phyla, "Fungal phyla"),
  ncol = 1,
  rel_heights = c(1.00, 0.58)
)

combined_panels <- plot_grid(
  bacteria_panel,
  fungi_panel,
  nrow = 1,
  rel_widths = c(1, 1),
  labels = c("(a)", "(b)"),
  label_size = 15,
  label_fontface = "plain",
  label_x = c(0.01, 0.01),
  label_y = c(0.99, 0.99),
  hjust = 0,
  vjust = 1
)

combined_figure <- plot_grid(
  combined_panels,
  combined_legend,
  nrow = 1,
  rel_widths = c(1.15, 0.48),
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

save_plot_pair("bacteria_60_samples_relative_abundance", bacteria_figure, 8.6, 4.2)
save_plot_pair("fungi_60_samples_relative_abundance", fungi_figure, 7.8, 4.2)
save_plot_pair("bacteria_fungi_60_samples_relative_abundance_combined", combined_figure, 13.2, 4.45)

message("Created 60-sample relative abundance stacked bar plots without significance tests.")
