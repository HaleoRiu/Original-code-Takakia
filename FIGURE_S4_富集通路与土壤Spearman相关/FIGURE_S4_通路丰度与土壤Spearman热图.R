#!/usr/bin/env Rscript

# R/E-enriched KEGG pathway abundance (sum of matched KO TPM values) versus
# soil physicochemical properties. Ordinary Spearman correlation; altitude is
# deliberately not controlled, following the revised figure specification.

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(cowplot)
})

out_dir <- "/Users/haleoriu/Documents/Takakia绘图和文章修改/enriched_pathway_soil_partial_spearman_2026-07-30"
pathway_file <- "/Volumes/Haleo_DK01/WD_Takakia_F/metagenome-data/Figure4/Figure4_integrated_pathway_abundance_by_sample.csv"
soil_file <- "/Volumes/Haleo_DK01/WD_Takakia_F/amplicon-data/psychem.txt"

soil_vars <- c("SOC", "TN", "TP", "AN", "AP", "AvCa", "AvK", "AvMg", "AvFe", "TFe", "pH", "Moisture")
soil_labels <- c(AvCa = "Ca_ex", AvK = "K_ex", AvMg = "Mg_ex", AvFe = "Fe_ex")
# Top-to-bottom pathway order in Figure 3g (R pathways first, then E pathways).
figure3_order <- c(
  "ko01210", "ko03010", "ko00020", "ko00785", "ko00190", "ko00720", "ko01200", "ko03018", "ko00541", "ko00970",
  "ko00860", "ko00571", "ko00984", "ko00572", "ko02040", "ko00074", "ko02010", "ko01220", "ko00624"
)

pathways <- read.csv(pathway_file, check.names = FALSE, stringsAsFactors = FALSE)
soil <- read.delim(soil_file, check.names = FALSE, stringsAsFactors = FALSE) %>%
  mutate(PairID = sub("^PS", "", Sample))
joined <- pathways %>%
  select(Direction, ID, Description, sample, PathwayAbundance, PairID, Group) %>%
  inner_join(soil %>% select(PairID, all_of(soil_vars)), by = "PairID") %>%
  # Each pathway exists in both habitats in the source data.  Correlate an
  # R-enriched pathway only in R samples, and an E-enriched pathway only in E
  # samples, matching the corresponding comparison-enrichment heatmap cohort.
  filter((Direction == "Rhizomatousphere enriched" & Group == "R") |
         (Direction == "Endophytic enriched" & Group == "E"))
stopifnot(all(table(joined$Direction) %% 30 == 0), all(table(paste(joined$Direction, joined$ID)) == 30))

catalogue <- joined %>% distinct(Direction, ID, Description) %>%
  mutate(Figure3_order = match(ID, figure3_order)) %>%
  arrange(Figure3_order) %>%
  select(-Figure3_order)
results <- bind_rows(lapply(seq_len(nrow(catalogue)), function(i) {
  info <- catalogue[i, ]
  d <- joined %>% filter(Direction == info$Direction, ID == info$ID) %>% arrange(PairID)
  bind_rows(lapply(soil_vars, function(v) {
    z <- suppressWarnings(cor.test(d$PathwayAbundance, d[[v]], method = "spearman", exact = FALSE))
    data.frame(Direction = info$Direction, ID = info$ID, Description = info$Description,
               Soil_variable = v, rho = unname(z$estimate), p = z$p.value)
  }))
})) %>%
  group_by(Direction, ID) %>%
  mutate(p_BH_within_pathway = p.adjust(p, method = "BH"),
         stars_BH = case_when(p_BH_within_pathway < 0.001 ~ "***", p_BH_within_pathway < 0.01 ~ "**",
                               p_BH_within_pathway < 0.05 ~ "*", TRUE ~ ""),
         abs_rho = abs(rho),
         Soil_label = dplyr::recode(Soil_variable, !!!soil_labels)) %>%
  ungroup()

labels <- catalogue %>% mutate(Pathway_label = vapply(Description, function(x) paste(strwrap(x, width = 29), collapse = "\n"), character(1)))
results <- results %>% left_join(labels %>% select(Direction, ID, Pathway_label), by = c("Direction", "ID"))

make_heatmap <- function(direction) {
  d <- results %>% filter(Direction == direction)
  ordered_labels <- labels %>% filter(Direction == direction) %>% pull(Pathway_label)
  d$Pathway_label <- factor(d$Pathway_label, levels = rev(ordered_labels))
  d$Soil_label <- factor(d$Soil_label, levels = c("SOC", "TN", "TP", "AN", "AP", "Ca_ex", "K_ex", "Mg_ex", "Fe_ex", "TFe", "pH", "Moisture"))
  title <- if (direction == "Rhizomatousphere enriched") "R-enriched pathways" else "E-enriched pathways"
  ggplot(d, aes(Soil_label, Pathway_label, fill = rho)) +
    geom_tile(colour = "white", linewidth = 0.55) +
    geom_text(aes(label = stars_BH), family = "Helvetica", size = 4.2, fontface = "bold") +
    scale_fill_gradient2(low = "#762A83", mid = "#F7F7F7", high = "#1B7837", midpoint = 0,
                         limits = c(-1, 1), name = "Spearman's\nrho") +
    labs(title = title, x = NULL, y = NULL,
         subtitle = "Summed pathway abundance; BH within pathway") +
    theme_bw(base_size = 12, base_family = "Helvetica") +
    theme(panel.grid = element_blank(), panel.border = element_rect(linewidth = 0.75),
          plot.title = element_text(face = "bold", hjust = 0.5, size = 14),
          plot.subtitle = element_text(hjust = 0.5, size = 10.2),
          axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, colour = "black", size = 10.2),
          axis.text.y = element_text(colour = "black", size = 9.8),
          legend.title = element_text(face = "bold", size = 10), legend.text = element_text(size = 9),
          plot.margin = margin(4, 8, 3, 5))
}

p_r <- make_heatmap("Rhizomatousphere enriched") + theme(legend.position = "none")
p_e <- make_heatmap("Endophytic enriched")
legend <- get_legend(p_e + theme(legend.position = "right"))
p_e <- p_e + theme(legend.position = "none")
combined <- plot_grid(p_r, p_e, legend, nrow = 1, rel_widths = c(1, 1, 0.20), labels = c("C", "D"),
                      label_fontfamily = "Helvetica", label_fontface = "bold", label_size = 15)

ggsave(file.path(out_dir, "Fig_S4c_d_Spearman_summed_pathway_soil_BH.pdf"), combined, width = 15.6, height = 5.35, units = "in")
ggsave(file.path(out_dir, "Fig_S4c_d_Spearman_summed_pathway_soil_BH.png"), combined, width = 15.6, height = 5.35, units = "in", dpi = 600)
write.table(results %>% select(-Pathway_label), file.path(out_dir, "Spearman_summed_pathway_soil_BH_results.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
writeLines(c(
  "Spearman correlations between summed KO abundance within R/E-enriched pathways and soil physicochemical properties.",
  "R-enriched pathways use only R samples; E-enriched pathways use only E samples, matching the habitat-specific comparison-enrichment heatmap.",
  "No elevation adjustment was applied.",
  "BH correction was applied within each pathway across its 12 soil variables."
), file.path(out_dir, "README_Spearman_summed_pathway_soil.md"))
message("Completed ", nrow(results), " pathway-by-soil Spearman tests.")
