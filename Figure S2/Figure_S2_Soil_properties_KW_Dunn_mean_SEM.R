#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(patchwork)
  library(rstatix)
  library(multcompView)
})

input_file <- "/Volumes/Haleo_DK01/WD_Takakia_F/amplicon-data/psychem.txt"
project_dir <- "/Volumes/Haleo_DK01/WD_Takakia_F/amplicon-data/physicochemical_altitude_2026-08-02_KW_Dunn"
output_dir <- file.path(project_dir, "output")
table_dir <- file.path(project_dir, "tables")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)

altitude_levels <- c("3800", "3900", "4000", "4100", "4200")
palette_altitude <- c(
  "3800" = "#F29AC2",
  "3900" = "#96B8F1",
  "4000" = "#F7D99B",
  "4100" = "#96EF91",
  "4200" = "#CCCCCC"
)

variable_info <- tibble::tribble(
  ~variable,   ~axis_label,
  "SOC",       "SOC (g/kg)",
  "TN",        "TN (g/kg)",
  "TP",        "TP (g/kg)",
  "AN",        "AN (mg/kg)",
  "AP",        "AP (mg/kg)",
  "AvCa",      "AvCa (mg/kg)",
  "AvK",       "AvK (mg/kg)",
  "AvMg",      "AvMg (mg/kg)",
  "AvFe",      "AvFe (mg/kg)",
  "TFe",       "TFe (mg/kg)",
  "pH",        "pH",
  "Moisture",  "Moisture (%)"
)

raw <- read.delim(
  input_file,
  header = TRUE,
  sep = "\t",
  check.names = FALSE,
  stringsAsFactors = FALSE
)

required_columns <- c("Sample", variable_info$variable, "elev")
missing_columns <- setdiff(required_columns, names(raw))
if (length(missing_columns) > 0) {
  stop("Missing required columns: ", paste(missing_columns, collapse = ", "))
}
if (anyDuplicated(raw$Sample)) stop("Duplicated Sample identifiers detected.")
if (anyNA(raw[, required_columns])) stop("Missing values detected in required columns.")
if (!setequal(as.character(unique(raw$elev)), altitude_levels)) {
  stop("Unexpected altitude groups: ", paste(sort(unique(raw$elev)), collapse = ", "))
}

raw$elev <- factor(as.character(raw$elev), levels = altitude_levels, ordered = TRUE)

group_counts <- raw %>% count(elev, name = "n")
write.table(
  group_counts,
  file.path(table_dir, "physicochemical_group_counts.tsv"),
  sep = "\t", row.names = FALSE, quote = FALSE
)

long <- raw %>%
  select(Sample, elev, all_of(variable_info$variable)) %>%
  pivot_longer(
    cols = all_of(variable_info$variable),
    names_to = "variable",
    values_to = "value"
  ) %>%
  left_join(variable_info, by = "variable")

summary_stats <- long %>%
  group_by(variable, axis_label, elev) %>%
  summarise(
    n = n(),
    mean = mean(value),
    sd = sd(value),
    sem = sd / sqrt(n),
    median = median(value),
    q1 = quantile(value, 0.25),
    q3 = quantile(value, 0.75),
    .groups = "drop"
  )

diagnostics <- lapply(variable_info$variable, function(v) {
  fit <- aov(raw[[v]] ~ raw$elev)
  shapiro <- shapiro.test(residuals(fit))
  fligner <- fligner.test(raw[[v]] ~ raw$elev)
  tibble(
    variable = v,
    shapiro_residual_W = unname(shapiro$statistic),
    shapiro_residual_p = shapiro$p.value,
    fligner_killeen_chi2 = unname(fligner$statistic),
    fligner_killeen_p = fligner$p.value
  )
}) %>% bind_rows()

kruskal_results <- lapply(variable_info$variable, function(v) {
  test <- kruskal.test(raw[[v]] ~ raw$elev)
  tibble(
    variable = v,
    n_total = sum(!is.na(raw[[v]])),
    groups = nlevels(raw$elev),
    df = unname(test$parameter),
    kruskal_wallis_chi2 = unname(test$statistic),
    p_raw = test$p.value,
    significance = case_when(
      test$p.value < 0.001 ~ "***",
      test$p.value < 0.01 ~ "**",
      test$p.value < 0.05 ~ "*",
      TRUE ~ "ns"
    )
  )
}) %>% bind_rows()

dunn_all <- list()
letter_all <- list()

for (v in variable_info$variable) {
  means_desc <- summary_stats %>%
    filter(variable == v) %>%
    arrange(desc(mean)) %>%
    pull(elev) %>%
    as.character()

  test_data <- raw %>%
    transmute(value = .data[[v]], altitude_for_test = factor(elev, levels = means_desc))

  dunn <- rstatix::dunn_test(
    test_data,
    value ~ altitude_for_test,
    p.adjust.method = "BH"
  ) %>%
    transmute(
      variable = v,
      group1,
      group2,
      n1,
      n2,
      statistic,
      p_raw = p,
      p_BH_within_variable = p.adj
    )

  p_vector <- setNames(
    dunn$p_BH_within_variable,
    paste(dunn$group1, dunn$group2, sep = "-")
  )
  compact_letters <- multcompView::multcompLetters(
    p_vector,
    threshold = 0.05,
    compare = "<"
  )$Letters

  dunn_all[[v]] <- dunn
  letter_all[[v]] <- tibble(
    variable = v,
    elev = factor(names(compact_letters), levels = altitude_levels, ordered = TRUE),
    letter = unname(compact_letters)
  )
}

dunn_results <- bind_rows(dunn_all) %>%
  arrange(match(variable, variable_info$variable), group1, group2)

letters <- bind_rows(letter_all)

summary_stats <- summary_stats %>%
  left_join(letters, by = c("variable", "elev")) %>%
  group_by(variable) %>%
  mutate(
    data_max = max(mean + sem),
    data_min = min(0, min(mean - sem)),
    data_range = max(data_max - data_min, abs(data_max) * 0.10, 1e-9),
    letter_y = mean + sem + 0.045 * data_range
  ) %>%
  ungroup()

write.table(
  raw,
  file.path(table_dir, "physicochemical_source_data.tsv"),
  sep = "\t", row.names = FALSE, quote = FALSE
)
write.table(
  summary_stats %>% select(-data_max, -data_min, -data_range, -letter_y),
  file.path(table_dir, "physicochemical_mean_SEM_and_letters.tsv"),
  sep = "\t", row.names = FALSE, quote = FALSE
)
write.table(
  diagnostics,
  file.path(table_dir, "physicochemical_distribution_variance_diagnostics.tsv"),
  sep = "\t", row.names = FALSE, quote = FALSE
)
write.table(
  kruskal_results,
  file.path(table_dir, "physicochemical_Kruskal_Wallis_omnibus.tsv"),
  sep = "\t", row.names = FALSE, quote = FALSE
)
write.table(
  dunn_results,
  file.path(table_dir, "physicochemical_Dunn_pairwise_BH_within_variable.tsv"),
  sep = "\t", row.names = FALSE, quote = FALSE
)

make_panel <- function(v) {
  panel_data <- summary_stats %>% filter(variable == v)
  label <- variable_info$axis_label[match(v, variable_info$variable)]

  ggplot(panel_data, aes(x = elev, y = mean, fill = elev)) +
    geom_col(
      width = 0.58,
      colour = "black",
      linewidth = 0.45,
      show.legend = FALSE
    ) +
    geom_errorbar(
      aes(ymin = pmax(0, mean - sem), ymax = mean + sem),
      width = 0.18,
      linewidth = 0.50
    ) +
    geom_text(
      aes(y = letter_y, label = letter),
      family = "Helvetica",
      size = 4.1,
      vjust = 0
    ) +
    scale_fill_manual(values = palette_altitude, drop = FALSE) +
    scale_y_continuous(
      limits = c(0, NA),
      expand = expansion(mult = c(0, 0.15))
    ) +
    labs(x = NULL, y = label) +
    theme_classic(base_family = "Helvetica", base_size = 11) +
    theme(
      axis.title.y = element_text(face = "bold", size = 12, margin = margin(r = 7)),
      axis.text.x = element_text(size = 10.5, colour = "black"),
      axis.text.y = element_text(size = 10.5, colour = "black"),
      axis.ticks = element_line(linewidth = 0.45, colour = "black"),
      axis.line = element_line(linewidth = 0.50, colour = "black"),
      plot.margin = margin(8, 8, 8, 8)
    )
}

panels <- lapply(variable_info$variable, make_panel)
combined_plot <- wrap_plots(panels, ncol = 4, nrow = 3, byrow = TRUE) +
  plot_layout(guides = "collect") +
  plot_annotation(
    caption = "Bars show mean +/- SEM (n = 6 per altitude). Kruskal-Wallis omnibus tests were followed by Dunn pairwise tests; pairwise P values were BH-adjusted within each variable (q < 0.05).",
    theme = theme(
      plot.caption = element_text(
        family = "Helvetica",
        size = 9.2,
        colour = "black",
        hjust = 0.5,
        margin = margin(t = 8)
      )
    )
  )

pdf_file <- file.path(output_dir, "physicochemical_altitude_mean_SEM_KW_Dunn_BH_within_variable.pdf")
png_file <- file.path(output_dir, "physicochemical_altitude_mean_SEM_KW_Dunn_BH_within_variable.png")
tiff_file <- file.path(output_dir, "physicochemical_altitude_mean_SEM_KW_Dunn_BH_within_variable.tiff")

ggsave(
  pdf_file,
  combined_plot,
  width = 14.8,
  height = 10.2,
  units = "in",
  device = grDevices::pdf,
  useDingbats = FALSE,
  bg = "white"
)
ggsave(
  png_file,
  combined_plot,
  width = 14.8,
  height = 10.2,
  units = "in",
  dpi = 600,
  device = ragg::agg_png,
  bg = "white"
)
ggsave(
  tiff_file,
  combined_plot,
  width = 14.8,
  height = 10.2,
  units = "in",
  dpi = 600,
  device = ragg::agg_tiff,
  compression = "lzw",
  bg = "white"
)

session_info <- capture.output(sessionInfo())
writeLines(session_info, file.path(project_dir, "R_sessionInfo.txt"))

readme <- c(
  "Physicochemical properties across altitude",
  "==========================================",
  "",
  paste0("Input: ", input_file),
  "Biological replication: n = 6 independent plots per altitude; 5 altitude groups.",
  "Display: arithmetic mean +/- SEM, matching the supplied reference layout.",
  "Inferential analysis: a Kruskal-Wallis omnibus test is performed separately for each variable, followed by Dunn rank-based pairwise comparisons among the five altitudes.",
  "Multiplicity: Benjamini-Hochberg correction is applied separately within each variable (10 pairwise comparisons per variable).",
  "The 12 Kruskal-Wallis omnibus P values are reported as raw P values; no correction is applied across physicochemical variables.",
  "Letters in the figure are determined only from the within-variable BH-adjusted pairwise p-values.",
  "Reason for the rank-based workflow: n = 6 per group and residual normality/variance diagnostics were not satisfied for several variables (notably AvCa, AvMg, AvFe and TFe).",
  "The plotted mean +/- SEM is descriptive; Dunn inference is rank-based.",
  "",
  "Outputs:",
  paste0("- ", basename(pdf_file)),
  paste0("- ", basename(png_file), " (600 dpi)"),
  paste0("- ", basename(tiff_file), " (600 dpi, LZW)"),
  "- tables/physicochemical_source_data.tsv",
  "- tables/physicochemical_mean_SEM_and_letters.tsv",
  "- tables/physicochemical_Kruskal_Wallis_omnibus.tsv",
  "- tables/physicochemical_Dunn_pairwise_BH_within_variable.tsv",
  "- tables/physicochemical_distribution_variance_diagnostics.tsv",
  "- R_sessionInfo.txt"
)
writeLines(readme, file.path(project_dir, "README_methods.txt"))

message("Completed. Outputs written to: ", project_dir)
