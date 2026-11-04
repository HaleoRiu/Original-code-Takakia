#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(grid)
})

args <- commandArgs(trailingOnly = TRUE)
prefix <- ifelse(length(args) >= 1, args[[1]], "Figure5_antistress_GSHaa_main")

base_dir <- "/Volumes/Haleo_DK01/WD_Takakia_F/metagenome-data"
out_dir <- file.path(base_dir, "Figure5")
input_file <- file.path(out_dir, paste0(prefix, "_category_TPM_by_sample.tsv"))
summary_file <- file.path(out_dir, paste0(prefix, "_category_summary.tsv"))
out_stem <- paste0(prefix, "_boxplots_3x4_pairedWilcoxBH_KW")
stats_file <- file.path(out_dir, paste0(out_stem, "_statistics_raw_p_and_BH_q.csv"))

OUT_WIDTH <- 10.8
OUT_HEIGHT <- 9.8
OUT_DPI <- 300

group_colours <- c("R" = "#82B3D7", "E" = "#CF6C72")
group_labels <- c("R" = "Rhizomatousphere, R", "E" = "Endophytic, E")
plot_levels_all <- c("All samples", "3800", "3900", "4000", "4100", "4200")
plot_labels_all <- c("All", "3800", "3900", "4000", "4100", "4200")

category_levels <- c(
  "All anti-stress genes",
  "DNA repair",
  "Glutathione and related amino acids",
  "Cold stress",
  "Flavonoids",
  "Terpenoids",
  "VB12",
  "Carotenoids",
  "Iron carrier",
  "Colonization",
  "Oxidative stress",
  "Trehalose"
)

category_labels <- c(
  "All anti-stress genes" = "All anti-stress functional genes",
  "DNA repair" = "UV/oxidative DNA repair",
  "Glutathione and related amino acids" = "Glutathione and related\namino acids",
  "Cold stress" = "Cold shock response",
  "Flavonoids" = "Flavonoid biosynthesis",
  "Terpenoids" = "Terpenoid backbone biosynthesis",
  "VB12" = "Cobalamin (VB12) biosynthesis",
  "Carotenoids" = "Carotenoid biosynthesis",
  "Iron carrier" = "Siderophore biosynthesis",
  "Colonization" = "Quorum sensing-mediated\ncolonization",
  "Oxidative stress" = "Oxidative stress",
  "Trehalose" = "Trehalose metabolism"
)

pt_to_lw <- function(pt) pt / 2.845276
lw_box <- pt_to_lw(0.50)
lw_axis <- pt_to_lw(0.80)
lw_sig <- pt_to_lw(0.65)

sig_symbol <- function(q) {
  dplyr::case_when(
    is.na(q) ~ "ns",
    q < 0.001 ~ "***",
    q < 0.01 ~ "**",
    q < 0.05 ~ "*",
    TRUE ~ "ns"
  )
}

sig_bool <- function(x) {
  dplyr::case_when(
    is.na(x) ~ FALSE,
    x < 0.05 ~ TRUE,
    TRUE ~ FALSE
  )
}

format_p <- function(p) {
  ifelse(is.na(p), "NA", ifelse(p < 0.0001, "<0.0001", sprintf("%.4f", p)))
}

category_long <- read.delim(input_file, check.names = FALSE) %>%
  mutate(
    Group = factor(Group, levels = c("R", "E")),
    category = factor(category, levels = category_levels),
    Elevation = as.integer(Elevation),
    PairID = paste0(Elevation / 100, "_", Replicate)
  )

nko_tbl <- read.delim(summary_file, check.names = FALSE) %>%
  mutate(category = factor(category, levels = category_levels)) %>%
  select(category, n_KO)

plot_df <- bind_rows(
  category_long %>% mutate(XGroup = "All samples"),
  category_long %>% mutate(XGroup = as.character(Elevation))
) %>%
  mutate(
    XGroup = factor(XGroup, levels = plot_levels_all),
    X = as.numeric(XGroup),
    category = factor(category, levels = category_levels)
  ) %>%
  left_join(nko_tbl, by = "category")

# R/E paired Wilcoxon tests: All samples and each altitude.
re_all_stats <- category_long %>%
  select(category, PairID, Group, TPM) %>%
  pivot_wider(names_from = Group, values_from = TPM) %>%
  group_by(category) %>%
  summarise(
    test_type = "R_vs_E_All_paired_Wilcoxon",
    group = "R_vs_E",
    comparison = "All samples",
    p_value = tryCatch(
      wilcox.test(R, E, paired = TRUE, exact = FALSE)$p.value,
      error = function(e) NA_real_
    ),
    .groups = "drop"
  ) %>%
  mutate(
    adjustment_family = "R_vs_E_All_across_categories",
    q_value = p.adjust(p_value, method = "BH"),
    q_global = q_value,
    q_significance = sig_symbol(q_value),
    q_is_significant = sig_bool(q_value)
  )

re_altitude_stats <- category_long %>%
  mutate(comparison = as.character(Elevation)) %>%
  select(category, comparison, Replicate, Group, TPM) %>%
  pivot_wider(names_from = Group, values_from = TPM) %>%
  group_by(category, comparison) %>%
  summarise(
    test_type = "R_vs_E_each_altitude_paired_Wilcoxon",
    group = "R_vs_E",
    p_value = tryCatch(
      wilcox.test(R, E, paired = TRUE, exact = FALSE)$p.value,
      error = function(e) NA_real_
    ),
    .groups = "drop"
  ) %>%
  group_by(category) %>%
  mutate(
    adjustment_family = paste0("R_vs_E_each_altitude_within_category_", as.character(category)),
    q_value = p.adjust(p_value, method = "BH"),
    q_significance = sig_symbol(q_value),
    q_is_significant = sig_bool(q_value)
  ) %>%
  ungroup() %>%
  mutate(
    q_global = p.adjust(p_value, method = "BH")
  )

# R/E separate altitude effects: Kruskal-Wallis tests across elevations.
kw_stats <- category_long %>%
  mutate(ElevationFactor = factor(Elevation, levels = c(3800, 3900, 4000, 4100, 4200))) %>%
  group_by(category, Group) %>%
  summarise(
    test_type = "Altitude_Kruskal_Wallis",
    group = as.character(first(Group)),
    comparison = "3800_vs_3900_vs_4000_vs_4100_vs_4200",
    statistic = tryCatch(
      unname(kruskal.test(TPM ~ ElevationFactor)$statistic),
      error = function(e) NA_real_
    ),
    df = tryCatch(
      unname(kruskal.test(TPM ~ ElevationFactor)$parameter),
      error = function(e) NA_real_
    ),
    p_value = tryCatch(
      kruskal.test(TPM ~ ElevationFactor)$p.value,
      error = function(e) NA_real_
    ),
    .groups = "drop"
  ) %>%
  group_by(group) %>%
  mutate(
    adjustment_family = paste0("Kruskal_altitude_", group, "_across_categories"),
    q_value = p.adjust(p_value, method = "BH"),
    q_global = q_value,
    p_significance = sig_symbol(p_value),
    q_significance = sig_symbol(q_value),
    p_is_significant = sig_bool(p_value),
    q_is_significant = sig_bool(q_value)
  ) %>%
  ungroup()

stat_table <- bind_rows(
  re_all_stats %>%
    mutate(statistic = NA_real_, df = NA_real_, p_significance = sig_symbol(p_value),
           p_is_significant = sig_bool(p_value)) %>%
    select(category, test_type, group, comparison, statistic, df, p_value,
           adjustment_family, q_value, q_global, p_significance, q_significance,
           p_is_significant, q_is_significant),
  re_altitude_stats %>%
    mutate(statistic = NA_real_, df = NA_real_, p_significance = sig_symbol(p_value),
           p_is_significant = sig_bool(p_value)) %>%
    select(category, test_type, group, comparison, statistic, df, p_value,
           adjustment_family, q_value, q_global, p_significance, q_significance,
           p_is_significant, q_is_significant),
  kw_stats %>%
    select(category, test_type, group, comparison, statistic, df, p_value,
           adjustment_family, q_value, q_global, p_significance, q_significance,
           p_is_significant, q_is_significant)
) %>%
  mutate(category = factor(category, levels = category_levels)) %>%
  arrange(test_type, category, group, comparison)

write.csv(stat_table, stats_file, row.names = FALSE)

y_extent <- plot_df %>%
  group_by(category) %>%
  summarise(
    y_min = min(TPM, na.rm = TRUE),
    y_max = max(TPM, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    y_range = if_else(y_max > y_min, y_max - y_min, pmax(abs(y_max), 1) * 0.10),
    y_lower = y_min - y_range * 0.12,
    y_upper = y_max + y_range * 0.42
  )

make_re_annot <- function(source_stats, all_samples = FALSE) {
  base <- plot_df %>%
    filter(if (all_samples) XGroup == "All samples" else XGroup != "All samples") %>%
    group_by(category, XGroup, X) %>%
    summarise(y_max = max(TPM, na.rm = TRUE), .groups = "drop") %>%
    left_join(y_extent %>% select(category, y_range), by = "category") %>%
    mutate(
      bracket_y = y_max + y_range * 0.075,
      tick_y = y_max + y_range * 0.045,
      label_y = y_max + y_range * 0.115,
      x_left = X - 0.18,
      x_right = X + 0.18
    )
  if (all_samples) {
    base %>%
      left_join(
        source_stats %>%
          mutate(label = q_significance) %>%
          select(category, p_value, q_value, q_global, label),
        by = "category"
      )
  } else {
    base %>%
      mutate(comparison = as.character(XGroup)) %>%
      left_join(
        source_stats %>%
          mutate(label = q_significance) %>%
          select(category, comparison, p_value, q_value, q_global, label),
        by = c("category", "comparison")
      )
  }
}

re_annot <- bind_rows(
  make_re_annot(re_all_stats, all_samples = TRUE),
  make_re_annot(re_altitude_stats, all_samples = FALSE)
) %>%
  filter(!is.na(label), label != "ns")

kw_annot <- kw_stats %>%
  left_join(y_extent, by = "category") %>%
  mutate(
    x = if_else(group == "R", 2.82, 6.30),
    y = if_else(group == "R", y_upper - y_range * 0.08, y_upper - y_range * 0.16),
    hjust_value = if_else(group == "R", 0, 1),
    label = paste0(group, "~Kruskal~italic(q)~' = ", format_p(q_value), "'")
  )

plot_theme <- theme_classic(base_size = 11, base_family = "Helvetica") +
  theme(
    panel.border = element_rect(colour = "black", fill = NA, linewidth = lw_axis),
    axis.line = element_blank(),
    axis.ticks = element_line(colour = "black", linewidth = lw_axis),
    axis.ticks.length = unit(0.08, "cm"),
    axis.text.x = element_text(colour = "black", size = 9),
    axis.text.y = element_text(colour = "black", size = 7.6, angle = 90, hjust = 0.5, vjust = 0.5),
    axis.title.x = element_text(colour = "black", size = 16, face = "bold", margin = margin(t = 2)),
    axis.title.y = element_text(colour = "black", size = 16, face = "bold", margin = margin(r = 3)),
    strip.background = element_blank(),
    strip.text = element_text(colour = "black", size = 10.2, face = "bold", margin = margin(2, 0, 3, 0)),
    legend.position = "top",
    legend.title = element_blank(),
    legend.text = element_text(size = 12),
    legend.key.size = unit(0.35, "cm"),
    legend.spacing.x = unit(0.14, "cm"),
    panel.spacing.x = unit(0.28, "cm"),
    panel.spacing.y = unit(0.30, "cm"),
    plot.margin = margin(2, 8, 2, 4),
    plot.background = element_rect(fill = "white", colour = NA)
  )

p <- ggplot(plot_df, aes(x = X, y = TPM, fill = Group)) +
  geom_boxplot(
    aes(group = interaction(X, Group)),
    position = position_dodge(width = 0.62),
    width = 0.48,
    linewidth = lw_box,
    colour = "black",
    outlier.shape = 16,
    outlier.size = 0.48,
    outlier.stroke = 0,
    alpha = 0.92
  ) +
  geom_vline(xintercept = 1.5, linetype = "dashed", colour = "#9E9E9E", linewidth = pt_to_lw(0.7)) +
  geom_segment(
    data = re_annot,
    aes(x = x_left, xend = x_right, y = bracket_y, yend = bracket_y),
    inherit.aes = FALSE,
    colour = "#9E9E9E",
    linewidth = lw_sig
  ) +
  geom_segment(
    data = re_annot,
    aes(x = x_left, xend = x_left, y = tick_y, yend = bracket_y),
    inherit.aes = FALSE,
    colour = "#9E9E9E",
    linewidth = lw_sig
  ) +
  geom_segment(
    data = re_annot,
    aes(x = x_right, xend = x_right, y = tick_y, yend = bracket_y),
    inherit.aes = FALSE,
    colour = "#9E9E9E",
    linewidth = lw_sig
  ) +
  geom_text(
    data = re_annot,
    aes(x = X, y = label_y, label = label),
    inherit.aes = FALSE,
    family = "Helvetica",
    fontface = "bold",
    size = 2.55,
    colour = "red"
  ) +
  geom_text(
    data = kw_annot,
    aes(x = x, y = y, label = label, hjust = hjust_value),
    inherit.aes = FALSE,
    family = "Helvetica",
    size = 2.50,
    colour = "black",
    parse = TRUE
  ) +
  scale_fill_manual(values = group_colours, labels = group_labels, drop = FALSE) +
  scale_x_continuous(
    breaks = seq_along(plot_levels_all),
    labels = plot_labels_all,
    limits = c(0.45, length(plot_levels_all) + 0.45),
    expand = expansion(mult = c(0, 0))
  ) +
  scale_y_continuous(
    breaks = function(x) pretty(x, n = 3),
    labels = function(x) format(x, trim = TRUE, scientific = FALSE),
    expand = expansion(mult = c(0.035, 0.10))
  ) +
  facet_wrap(~ category, ncol = 3, scales = "free_y", labeller = labeller(category = category_labels)) +
  labs(x = "Altitude (m)", y = "TPM") +
  guides(fill = guide_legend(nrow = 1, byrow = TRUE, override.aes = list(size = 4))) +
  plot_theme

ggsave(file.path(out_dir, paste0(out_stem, ".png")), p,
       width = OUT_WIDTH, height = OUT_HEIGHT, units = "in", dpi = OUT_DPI, bg = "white")
ggsave(file.path(out_dir, paste0(out_stem, ".pdf")), p,
       width = OUT_WIDTH, height = OUT_HEIGHT, units = "in", device = grDevices::pdf, bg = "white")
ggsave(file.path(out_dir, paste0(out_stem, ".svg")), p,
       width = OUT_WIDTH, height = OUT_HEIGHT, units = "in", device = grDevices::svg, bg = "white")

message("Done: ", file.path(out_dir, paste0(out_stem, ".png")))
message("Stats: ", stats_file)
