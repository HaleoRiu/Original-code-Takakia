#!/usr/bin/env Rscript

# Re-run significance tests for finished FIGURE 1-3 with a non-paired policy.
# Exclusions requested by the user:
#   - community assembly tests are not rerun;
#   - volcano plot tests are not rerun.
#
# Main replacement policy:
#   - R/E comparisons: unpaired Wilcoxon rank-sum tests;
#   - altitude effects: Kruskal-Wallis tests;
#   - altitude post-hoc: pairwise unpaired Wilcoxon tests with BH correction.
# Adonis/PERMANOVA rows are included as "unchanged" for bookkeeping only.

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
})

root <- "/Volumes/Haleo_DK01/WD_Takakia_F"
amp_dir <- file.path(root, "amplicon-data")
mg_dir <- file.path(root, "metagenome-data")
fig4_dir <- file.path(mg_dir, "Figure4")
out_dir <- file.path(root, "manuscript", "V3", "Final_version")
out_csv <- file.path(out_dir, "Figure1-3_nonpaired_statistical_recheck.csv")

sig_symbol <- function(q) {
  dplyr::case_when(
    is.na(q) ~ NA_character_,
    q < 0.001 ~ "***",
    q < 0.01 ~ "**",
    q < 0.05 ~ "*",
    TRUE ~ "ns"
  )
}

safe_wilcox <- function(x, y, paired = FALSE) {
  x <- x[is.finite(x)]
  y <- y[is.finite(y)]
  if (length(x) < 2 || length(y) < 2) return(NA_real_)
  tryCatch(wilcox.test(x, y, paired = paired, exact = FALSE)$p.value,
           error = function(e) NA_real_)
}

safe_kruskal <- function(value, group) {
  keep <- is.finite(value) & !is.na(group)
  value <- value[keep]
  group <- factor(group[keep])
  if (length(unique(group)) < 2 || length(unique(value)) <= 1) return(c(p = NA_real_, statistic = NA_real_, df = NA_real_))
  test <- tryCatch(kruskal.test(value ~ group), error = function(e) NULL)
  if (is.null(test)) return(c(p = NA_real_, statistic = NA_real_, df = NA_real_))
  c(p = test$p.value, statistic = unname(test$statistic), df = unname(test$parameter))
}

parse_sample_meta <- function(samples) {
  tibble(sample = samples) %>%
    mutate(
      Prefix = substr(sample, 1, 1),
      Group = recode(Prefix, "S" = "R", "P" = "E"),
      Habitat = recode(Group, "R" = "Rhizomatousphere", "E" = "Endophytic"),
      Elevation = as.integer(substr(sample, 2, 3)) * 100,
      Replicate = as.integer(sub(".*_", "", sample)),
      PairID = paste0(Elevation / 100, "_", Replicate)
    )
}

sample_cols_from_df <- function(df) {
  grep("^[SP][0-9]{2}_[0-9]+$", names(df), value = TRUE)
}

read_feature_table_long <- function(path, feature_col, value_name = "value") {
  dat <- read.delim(path, check.names = FALSE)
  sample_cols <- sample_cols_from_df(dat)
  meta <- parse_sample_meta(sample_cols)
  dat %>%
    pivot_longer(all_of(sample_cols), names_to = "sample", values_to = value_name) %>%
    left_join(meta, by = "sample") %>%
    rename(feature = all_of(feature_col))
}

add_common_cols <- function(df, figure, panel, data_source, previous_method, new_method, notes = "") {
  df %>%
    mutate(
      figure = figure,
      panel = panel,
      data_source = data_source,
      previous_method_in_plot = previous_method,
      rerun_method = new_method,
      notes = notes,
      significance = sig_symbol(q_value)
    ) %>%
    select(
      figure, panel, data_source, analysis_family, feature_type, feature_id, feature_name,
      metric, group, comparison, previous_method_in_plot, rerun_method,
      statistic, df, p_value, q_value, significance, adjustment_family, notes
    )
}

re_test_by_feature <- function(df, value_col, feature_cols, metric_col = NULL,
                               figure, panel, data_source, previous_method,
                               q_family_label, feature_type) {
  group_cols <- c(feature_cols, metric_col)
  res <- df %>%
    group_by(across(all_of(group_cols))) %>%
    summarise(
      p_value = safe_wilcox(.data[[value_col]][.data$Group == "R"], .data[[value_col]][.data$Group == "E"], paired = FALSE),
      R_median = median(.data[[value_col]][.data$Group == "R"], na.rm = TRUE),
      E_median = median(.data[[value_col]][.data$Group == "E"], na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      q_value = p.adjust(p_value, method = "BH"),
      statistic = NA_real_,
      df = NA_real_,
      analysis_family = q_family_label,
      feature_type = feature_type,
      feature_id = if ("ID" %in% names(.)) as.character(ID) else if ("feature" %in% names(.)) as.character(feature) else NA_character_,
      feature_name = if ("Description" %in% names(.)) as.character(Description) else if ("feature" %in% names(.)) as.character(feature) else NA_character_,
      metric = if (!is.null(metric_col)) as.character(.data[[metric_col]]) else NA_character_,
      group = "R_vs_E",
      comparison = "All samples",
      adjustment_family = q_family_label
    )
  add_common_cols(res, figure, panel, data_source, previous_method,
                  "unpaired Wilcoxon rank-sum + BH", "")
}

re_test_by_feature_xgroup <- function(df, value_col, feature_cols, metric_col,
                                      figure, panel, data_source, previous_method,
                                      q_family_label, feature_type) {
  res <- df %>%
    group_by(across(all_of(c(feature_cols, metric_col, "XGroup")))) %>%
    summarise(
      p_value = safe_wilcox(.data[[value_col]][.data$Group == "R"], .data[[value_col]][.data$Group == "E"], paired = FALSE),
      .groups = "drop"
    ) %>%
    group_by(across(all_of(c(feature_cols, metric_col)))) %>%
    mutate(q_value = p.adjust(p_value, method = "BH")) %>%
    ungroup() %>%
    mutate(
      statistic = NA_real_,
      df = NA_real_,
      analysis_family = q_family_label,
      feature_type = feature_type,
      feature_id = if ("ID" %in% names(.)) as.character(ID) else if ("feature" %in% names(.)) as.character(feature) else NA_character_,
      feature_name = if ("Description" %in% names(.)) as.character(Description) else if ("feature" %in% names(.)) as.character(feature) else NA_character_,
      metric = as.character(.data[[metric_col]]),
      group = "R_vs_E",
      comparison = as.character(XGroup),
      adjustment_family = paste0(q_family_label, "_within_feature_metric")
    )
  add_common_cols(res, figure, panel, data_source, previous_method,
                  "unpaired Wilcoxon rank-sum + BH within feature/metric", "")
}

kw_by_group <- function(df, value_col, feature_cols, metric_col = NULL,
                        figure, panel, data_source, previous_method,
                        q_family_label, feature_type, selected_only = FALSE) {
  value_sym <- rlang::sym(value_col)
  group_cols <- c(feature_cols, metric_col, "Group")
  res <- df %>%
    group_by(across(all_of(group_cols))) %>%
    group_modify(function(.x, .g) {
      out <- safe_kruskal(.x[[value_col]], .x$Elevation)
      tibble(p_value = out[["p"]], statistic = out[["statistic"]], df = out[["df"]])
    }) %>%
    ungroup() %>%
    mutate(q_value = p.adjust(p_value, method = "BH")) %>%
    mutate(
      analysis_family = q_family_label,
      feature_type = feature_type,
      feature_id = if ("ID" %in% names(.)) as.character(ID) else if ("feature" %in% names(.)) as.character(feature) else NA_character_,
      feature_name = if ("Description" %in% names(.)) as.character(Description) else if ("feature" %in% names(.)) as.character(feature) else NA_character_,
      metric = if (!is.null(metric_col)) as.character(.data[[metric_col]]) else NA_character_,
      group = as.character(Group),
      comparison = "Altitude: 3800/3900/4000/4100/4200",
      adjustment_family = q_family_label
    )
  add_common_cols(res, figure, panel, data_source, previous_method,
                  "Kruskal-Wallis + BH", if (selected_only) "Only the enriched/selected habitat was tested." else "")
}

pairwise_altitude <- function(df, value_col, feature_cols, metric_col = NULL,
                              figure, panel, data_source, previous_method,
                              q_family_label, feature_type) {
  pairs <- combn(sort(unique(df$Elevation)), 2, simplify = FALSE)
  group_cols <- c(feature_cols, metric_col, "Group")
  res <- df %>%
    group_by(across(all_of(group_cols))) %>%
    group_modify(function(.x, .g) {
      rows <- lapply(pairs, function(pp) {
        tibble(
          Level1 = pp[1],
          Level2 = pp[2],
          p_value = safe_wilcox(.x[[value_col]][.x$Elevation == pp[1]], .x[[value_col]][.x$Elevation == pp[2]], paired = FALSE)
        )
      })
      bind_rows(rows) %>% mutate(q_value = p.adjust(p_value, method = "BH"))
    }) %>%
    ungroup() %>%
    mutate(
      statistic = NA_real_,
      df = NA_real_,
      analysis_family = q_family_label,
      feature_type = feature_type,
      feature_id = if ("ID" %in% names(.)) as.character(ID) else if ("feature" %in% names(.)) as.character(feature) else NA_character_,
      feature_name = if ("Description" %in% names(.)) as.character(Description) else if ("feature" %in% names(.)) as.character(feature) else NA_character_,
      metric = if (!is.null(metric_col)) as.character(.data[[metric_col]]) else NA_character_,
      group = as.character(Group),
      comparison = paste0(Level1, "_vs_", Level2),
      adjustment_family = paste0(q_family_label, "_within_feature_metric_group")
    )
  add_common_cols(res, figure, panel, data_source, previous_method,
                  "pairwise unpaired Wilcoxon + BH within feature/metric/group", "")
}

results <- list()

# FIGURE 1: amplicon phylum abundance, alpha diversity, and PCoA bookkeeping.
for (domain in c("bacteria", "fungi")) {
  phylum_file <- file.path(amp_dir, paste0("phylum_", domain, ".txt"))
  phylum_long <- read_feature_table_long(phylum_file, "phylum", "Count") %>%
    group_by(sample) %>%
    mutate(RelativeAbundance = Count / sum(Count, na.rm = TRUE) * 100) %>%
    ungroup()
  results[[paste0("fig1_", domain, "_phylum")]] <- re_test_by_feature(
    phylum_long, "RelativeAbundance", "feature", NULL,
    "FIGURE 1", paste0(domain, " phylum abundance"),
    phylum_file,
    "paired Wilcoxon + BH in previous phylum table",
    paste0("FIGURE1_", domain, "_phylum_R_vs_E_all"),
    "Phylum"
  )

  alpha_file <- file.path(amp_dir, paste0("div_", domain, ".txt"))
  alpha <- read.delim(alpha_file, check.names = FALSE) %>%
    rename(sample = sample) %>%
    left_join(parse_sample_meta(.$sample), by = "sample") %>%
    transmute(sample, Group, Habitat, Elevation, Replicate, PairID,
              Richness = sobs, Shannon = shannon) %>%
    pivot_longer(c(Richness, Shannon), names_to = "Metric", values_to = "Value")
  results[[paste0("fig1_", domain, "_alpha")]] <- re_test_by_feature(
    alpha, "Value", character(0), "Metric",
    "FIGURE 1", paste0(domain, " alpha diversity"),
    alpha_file,
    "paired Wilcoxon + BH in previous alpha table",
    paste0("FIGURE1_", domain, "_alpha_R_vs_E_all"),
    "Alpha diversity"
  )

  adonis_file <- file.path(amp_dir, paste0(domain, "_pcoa_adonis.csv"))
  if (file.exists(adonis_file)) {
    ad <- read.csv(adonis_file, check.names = FALSE)
    results[[paste0("fig1_", domain, "_adonis")]] <- tibble(
      figure = "FIGURE 1",
      panel = paste0(domain, " PCoA"),
      data_source = adonis_file,
      analysis_family = paste0("FIGURE1_", domain, "_PCoA_adonis_unchanged"),
      feature_type = "Ordination",
      feature_id = NA_character_,
      feature_name = paste0(domain, " PCoA"),
      metric = "Bray-Curtis",
      group = "R_vs_E",
      comparison = "PERMANOVA/Adonis group effect",
      previous_method_in_plot = "Adonis/PERMANOVA; not a paired Wilcoxon test",
      rerun_method = "unchanged; copied from existing adonis output",
      statistic = NA_real_,
      df = NA_real_,
      p_value = ad$p_value[1],
      q_value = ad$p_value[1],
      significance = sig_symbol(ad$p_value[1]),
      adjustment_family = "none",
      notes = "No paired-test replacement was applied."
    )
  }
}

# FIGURE 2: network plus community assembly; community assembly deliberately excluded.
results[["fig2_note"]] <- tibble(
  figure = "FIGURE 2",
  panel = "network and community assembly",
  data_source = file.path(amp_dir, "assembly_process_boxplot.R"),
  analysis_family = "FIGURE2_excluded_by_user_request",
  feature_type = "Community assembly",
  feature_id = NA_character_,
  feature_name = "Community assembly processes",
  metric = "Percentage",
  group = NA_character_,
  comparison = NA_character_,
  previous_method_in_plot = "Mann-Whitney/Wilcoxon for assembly process",
  rerun_method = "not rerun",
  statistic = NA_real_,
  df = NA_real_,
  p_value = NA_real_,
  q_value = NA_real_,
  significance = NA_character_,
  adjustment_family = NA_character_,
  notes = "Excluded because the user requested all tests except community assembly."
)

# FIGURE 3: metagenome MAG/functional panels.
func_alpha_file <- file.path(mg_dir, "div.txt")
func_alpha <- read.delim(func_alpha_file, check.names = FALSE) %>%
  rename(sample = `feature-id`) %>%
  left_join(parse_sample_meta(.$sample), by = "sample") %>%
  pivot_longer(c(Richness, Shannon), names_to = "Metric", values_to = "Value")
func_alpha_x <- bind_rows(
  func_alpha %>% mutate(XGroup = "All samples"),
  func_alpha %>% mutate(XGroup = as.character(Elevation))
)
results[["fig3_func_alpha_RE"]] <- re_test_by_feature_xgroup(
  func_alpha_x, "Value", character(0), "Metric",
  "FIGURE 3", "functional gene alpha diversity R/E by all samples and altitude",
  func_alpha_file,
  "paired Wilcoxon + BH in previous functional alpha figure",
  "FIGURE3_functional_alpha_R_vs_E_all_and_altitudes",
  "Functional gene alpha"
)
results[["fig3_func_alpha_KW"]] <- kw_by_group(
  func_alpha, "Value", character(0), "Metric",
  "FIGURE 3", "functional gene alpha diversity altitude effect",
  func_alpha_file,
  "Friedman across altitudes in previous altitude-letter versions",
  "FIGURE3_functional_alpha_altitude_Kruskal_Wallis",
  "Functional gene alpha"
)
results[["fig3_func_alpha_pairwise_alt"]] <- pairwise_altitude(
  func_alpha, "Value", character(0), "Metric",
  "FIGURE 3", "functional gene alpha diversity altitude post-hoc",
  func_alpha_file,
  "paired altitude post-hoc in previous versions",
  "FIGURE3_functional_alpha_altitude_pairwise_Wilcoxon",
  "Functional gene alpha"
)

for (item in list(
  list(file = file.path(mg_dir, "metagenome_redraw", "MAG_PCoA_adonis.csv"),
       panel = "MAG-level PCoA", metric = "MAG abundance Bray-Curtis"),
  list(file = file.path(mg_dir, "metagenome_redraw", "functional_gene_PCA_adonis.csv"),
       panel = "functional gene ordination", metric = "functional gene ordination")
)) {
  if (file.exists(item$file)) {
    ad <- read.csv(item$file, check.names = FALSE)
    model <- ad[1, ]
    results[[paste0("fig3_adonis_", gsub("[^A-Za-z]", "_", item$panel))]] <- tibble(
      figure = "FIGURE 3",
      panel = item$panel,
      data_source = item$file,
      analysis_family = paste0("FIGURE3_", item$panel, "_adonis_unchanged"),
      feature_type = "Ordination",
      feature_id = NA_character_,
      feature_name = item$panel,
      metric = item$metric,
      group = "R_vs_E",
      comparison = "PERMANOVA/Adonis group effect",
      previous_method_in_plot = "Adonis/PERMANOVA; not a paired Wilcoxon test",
      rerun_method = "unchanged; copied from existing adonis output",
      statistic = suppressWarnings(as.numeric(model$F)),
      df = suppressWarnings(as.numeric(model$Df)),
      p_value = suppressWarnings(as.numeric(model$`Pr(>F)`)),
      q_value = suppressWarnings(as.numeric(model$`Pr(>F)`)),
      significance = sig_symbol(suppressWarnings(as.numeric(model$`Pr(>F)`))),
      adjustment_family = "none",
      notes = "No paired-test replacement was applied."
    )
  }
}

# FIGURE 3 also contains the former Figure 4 pathway panels.
pathway_sample_file <- file.path(fig4_dir, "Figure4_integrated_pathway_abundance_by_sample.csv")
if (file.exists(pathway_sample_file)) {
  pathway <- read.csv(pathway_sample_file, check.names = FALSE) %>%
    mutate(
      SelectedGroup = if_else(grepl("Rhizomatousphere", Direction), "R", "E"),
      SelectedHabitat = if_else(SelectedGroup == "R", "Rhizomatousphere", "Endophytic"),
      Metric = "PathwayAbundance"
    )

  pathway_x <- bind_rows(
    pathway %>% mutate(XGroup = "All samples"),
    pathway %>% mutate(XGroup = as.character(Elevation))
  )
  results[["fig3_pathway_RE"]] <- re_test_by_feature_xgroup(
    pathway_x, "PathwayAbundance", c("ID", "Description", "Direction"), "Metric",
    "FIGURE 3", "comparative enrichment pathway abundance R/E",
    pathway_sample_file,
    "paired Wilcoxon in previous pathway abundance checks",
    "FIGURE3_pathway_R_vs_E_all_and_altitudes",
    "KEGG pathway"
  )

  pathway_selected <- pathway %>% filter(Group == SelectedGroup)
  results[["fig3_pathway_selected_KW"]] <- kw_by_group(
    pathway_selected, "PathwayAbundance", c("ID", "Description", "Direction"), "Metric",
    "FIGURE 3", "enriched-group pathway abundance altitude effect",
    pathway_sample_file,
    "Friedman across altitudes in previous enriched-group heatmap",
    "FIGURE3_pathway_selected_group_altitude_Kruskal_Wallis",
    "KEGG pathway",
    selected_only = TRUE
  )
  results[["fig3_pathway_selected_pairwise_alt"]] <- pairwise_altitude(
    pathway_selected, "PathwayAbundance", c("ID", "Description", "Direction"), "Metric",
    "FIGURE 3", "enriched-group pathway abundance altitude post-hoc",
    pathway_sample_file,
    "paired Wilcoxon post-hoc in previous enriched-group heatmap",
    "FIGURE3_pathway_selected_group_altitude_pairwise_Wilcoxon",
    "KEGG pathway"
  )
}

out <- bind_rows(results) %>%
  mutate(
    p_value = as.numeric(p_value),
    q_value = as.numeric(q_value),
    p_value_fmt = if_else(is.na(p_value), NA_character_, if_else(p_value < 1e-4, "<0.0001", sprintf("%.6f", p_value))),
    q_value_fmt = if_else(is.na(q_value), NA_character_, if_else(q_value < 1e-4, "<0.0001", sprintf("%.6f", q_value)))
  ) %>%
  arrange(figure, panel, analysis_family, feature_name, metric, group, comparison)

write.csv(out, out_csv, row.names = FALSE)
message("Wrote: ", out_csv)
message("Rows: ", nrow(out))
