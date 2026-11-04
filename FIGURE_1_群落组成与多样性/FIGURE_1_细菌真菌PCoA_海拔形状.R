#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(vegan)
  library(ggplot2)
  library(dplyr)
  library(cowplot)
})

args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
script_path <- if (length(file_arg)) {
  normalizePath(sub("^--file=", "", file_arg[1]))
} else {
  normalizePath("pcoa_95ci.R")
}
data_dir <- "/Volumes/Haleo_DK01/WD_Takakia_F/amplicon-data"

bacteria_file <- file.path(data_dir, "asv_bacteria.txt")
fungi_file <- file.path(data_dir, "asv_fungi.txt")

required_files <- c(bacteria_file, fungi_file)
if (!all(file.exists(required_files))) {
  stop(
    "Missing input file(s): ",
    paste(required_files[!file.exists(required_files)], collapse = ", ")
  )
}

group_colours <- c(
  "Rhizomatousphere" = "#82B3D7",
  "Endophytic" = "#CF6C72"
)

group_labels <- c(
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

read_asv <- function(path) {
  raw <- read.delim(
    path,
    header = TRUE,
    sep = "\t",
    row.names = 1,
    check.names = FALSE,
    stringsAsFactors = FALSE,
    comment.char = ""
  )

  sample_cols <- grep("^[SP][0-9]+_[0-9]+$", colnames(raw), value = TRUE)
  if (length(sample_cols) != 60) {
    stop(basename(path), " must contain 60 sample columns named like P38_1/S38_1.")
  }

  comm <- t(raw[, sample_cols, drop = FALSE])
  storage.mode(comm) <- "numeric"

  meta <- tibble(
    Sample = rownames(comm),
    Group = case_when(
      grepl("^S", Sample) ~ "Rhizomatousphere",
      grepl("^P", Sample) ~ "Endophytic",
      TRUE ~ NA_character_
    ),
    Altitude = factor(as.integer(substr(Sample, 2, 3)) * 100,
                      levels = c(3800, 3900, 4000, 4100, 4200))
  ) %>%
    filter(!is.na(Group)) %>%
    mutate(Group = factor(Group, levels = c("Rhizomatousphere", "Endophytic")))

  comm <- comm[meta$Sample, , drop = FALSE]

  list(comm = comm, meta = meta)
}

run_pcoa <- function(path) {
  prepared <- read_asv(path)
  comm <- prepared$comm
  meta <- prepared$meta

  bray <- vegdist(comm, method = "bray")
  pcoa <- cmdscale(bray, k = 2, eig = TRUE)
  eig <- pcoa$eig
  pc <- round(100 * eig / sum(eig[eig > 0]), 2)

  scores <- as.data.frame(pcoa$points[, 1:2, drop = FALSE])
  colnames(scores) <- c("PCoA1", "PCoA2")
  scores$Sample <- rownames(scores)

  plot_df <- left_join(scores, meta, by = "Sample")

  set.seed(20260625)
  ad <- adonis2(comm ~ Group, data = meta, permutations = 999, method = "bray")
  ad_stats <- tibble(
    R2 = ad$R2[1],
    p_value = ad$`Pr(>F)`[1]
  )

  list(plot_df = plot_df, pc = pc, adonis = ad_stats)
}

format_p <- function(p) {
  ifelse(p < 0.001, "0.001", sprintf("%.3f", p))
}

pcoa_theme <- function() {
  theme_bw(base_size = 12, base_family = "Helvetica") +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.95),
      axis.title = element_text(face = "bold", size = 13),
      axis.text = element_text(colour = "black", size = 11),
      plot.title = element_text(face = "bold", hjust = 0.5, size = 14),
      legend.position = "bottom",
      legend.title = element_blank(),
      legend.text = element_text(size = 9.2),
      legend.key.width = grid::unit(0.40, "cm"),
      legend.key.height = grid::unit(0.33, "cm"),
      legend.box = "vertical",
      legend.margin = margin(0, 0, 0, 0),
      legend.spacing.y = grid::unit(0.01, "cm"),
      legend.spacing.x = grid::unit(0.12, "cm"),
      plot.margin = margin(5, 5, 5, 5),
      panel.ontop = FALSE
    )
}

make_pcoa_plot <- function(result, title, show_legend = TRUE) {
  df <- result$plot_df
  ad <- result$adonis

  x_range <- diff(range(df$PCoA1, na.rm = TRUE))
  y_range <- diff(range(df$PCoA2, na.rm = TRUE))
  x_left <- min(df$PCoA1, na.rm = TRUE) + x_range * 0.03
  y_top <- max(df$PCoA2, na.rm = TRUE) - y_range * 0.03

  ggplot(df, aes(x = PCoA1, y = PCoA2)) +
    stat_ellipse(
      aes(fill = Group),
      type = "t",
      level = 0.95,
      geom = "polygon",
      alpha = 0.22,
      colour = NA,
      show.legend = FALSE
    ) +
    stat_ellipse(
      aes(colour = Group),
      type = "t",
      level = 0.95,
      linewidth = 1.05,
      linetype = "solid",
      show.legend = FALSE
    ) +
    geom_point(
      aes(fill = Group, colour = Group, shape = Altitude),
      size = 2.7,
      stroke = 0.35,
      alpha = 0.68
    ) +
    annotate(
      "rect",
      xmin = x_left - x_range * 0.01,
      xmax = x_left + x_range * 0.34,
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
      label = paste0("Adonis ", "R\u00B2", " = ", sprintf("%.3f", ad$R2)),
      hjust = 0,
      vjust = 1,
      family = "Helvetica",
      size = 4.2
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
      size = 4.2
    ) +
    annotate(
      "text",
      x = x_left + x_range * 0.13,
      y = y_top - y_range * 0.10,
      label = format_p(ad$p_value),
      hjust = 0,
      vjust = 1,
      family = "Helvetica",
      colour = "red",
      size = 4.2
    ) +
    scale_fill_manual(
      values = group_colours,
      breaks = names(group_labels),
      labels = group_labels,
      drop = FALSE
    ) +
    scale_colour_manual(
      values = group_colours,
      breaks = names(group_labels),
      labels = group_labels,
      drop = FALSE
    ) +
    scale_shape_manual(
      values = altitude_shapes,
      breaks = names(altitude_shapes),
      drop = FALSE
    ) +
    guides(
      fill = guide_legend(
        order = 1,
        nrow = 1,
        override.aes = list(shape = 22, size = 4, colour = "black", alpha = 1)
      ),
      colour = "none",
      shape = guide_legend(
        order = 2,
        nrow = 1,
        override.aes = list(size = 3.3, colour = "black", fill = "grey75", alpha = 1)
      )
    ) +
    labs(
      title = title,
      x = paste0("PCoA1 (", result$pc[1], "%)"),
      y = paste0("PCoA2 (", result$pc[2], "%)")
    ) +
    coord_cartesian(clip = "on") +
    pcoa_theme() +
    theme(
      legend.position = if (show_legend) "bottom" else "none",
      aspect.ratio = 0.82
    )
}

bacteria <- run_pcoa(bacteria_file)
fungi <- run_pcoa(fungi_file)

bacteria_plot <- make_pcoa_plot(bacteria, "Bacteria", show_legend = TRUE)
fungi_plot <- make_pcoa_plot(fungi, "Fungi", show_legend = TRUE)

shared_legend <- get_legend(
  make_pcoa_plot(bacteria, "Bacteria", show_legend = TRUE)
)

combined_panels <- plot_grid(
  make_pcoa_plot(bacteria, "Bacteria", show_legend = FALSE),
  make_pcoa_plot(fungi, "Fungi", show_legend = FALSE),
  ncol = 1,
  labels = c("(a)", "(b)"),
  label_size = 15,
  label_fontface = "plain",
  label_x = 0.01,
  label_y = 0.99,
  hjust = 0,
  vjust = 1,
  align = "v",
  axis = "lr"
)

combined_plot <- plot_grid(
  combined_panels,
  shared_legend,
  ncol = 1,
  rel_heights = c(1, 0.14)
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
  bacteria$adonis,
  file.path(data_dir, "bacteria_pcoa_adonis.csv"),
  row.names = FALSE
)

write.csv(
  fungi$adonis,
  file.path(data_dir, "fungi_pcoa_adonis.csv"),
  row.names = FALSE
)

save_plot_pair("bacteria_pcoa_95ci_nogrid_altitude_shape", bacteria_plot, 4.65, 4.35)
save_plot_pair("fungi_pcoa_95ci_nogrid_altitude_shape", fungi_plot, 4.65, 4.35)
save_plot_pair("bacteria_fungi_pcoa_95ci_nogrid_combined_altitude_shape", combined_plot, 4.75, 7.95)

message("Created PCoA plots with 95% CI ellipses, no internal grid/reference lines, and altitude-specific point shapes.")
