options(stringsAsFactors = FALSE, warn = 1)

script_arg <- grep("^--file=", commandArgs(), value = TRUE)
script_dir <- if (length(script_arg)) dirname(normalizePath(sub("^--file=", "", script_arg[1]))) else getwd()
project_dir <- normalizePath(file.path(script_dir, ".."), mustWork = FALSE)
work_root <- file.path(project_dir, "work")

n_source <- file.path(dirname(project_dir), "N_altitude_expanded5_2026-07-23", "work")
s_source <- file.path(dirname(project_dir), "S_altitude_cysEK_2026-07-23", "work")
base_feature_path <- "/Volumes/Haleo_DK01/WD_Takakia_F/metagenome-data/Figure7_CNPS/CNPS_EggNOG_split_results/prepared/CNPS_feature_renormalized_TPM_by_sample.tsv"
stopifnot(dir.exists(n_source), dir.exists(s_source), file.exists(base_feature_path))

sample_cols_for <- function(d) {
  out <- grep("^[PS](38|39|40|41|42)_[1-6]$", names(d), value = TRUE)
  stopifnot(length(out) == 60L)
  out
}

safe_wilcox <- function(x, y) {
  ok <- is.finite(x) & is.finite(y)
  x <- x[ok]
  y <- y[ok]
  if (length(x) < 2L || all(x == y)) return(1)
  out <- tryCatch(
    suppressWarnings(wilcox.test(x, y, paired = TRUE, exact = FALSE)$p.value),
    error = function(e) NA_real_
  )
  ifelse(is.finite(out), out, 1)
}

safe_kw <- function(x, g) {
  ok <- is.finite(x) & !is.na(g)
  x <- x[ok]
  g <- g[ok]
  if (length(unique(g)) < 2L || length(unique(x)) < 2L) return(1)
  out <- tryCatch(kruskal.test(x ~ factor(g))$p.value, error = function(e) NA_real_)
  ifelse(is.finite(out), out, 1)
}

star <- function(p, include_ns = FALSE) {
  out <- ifelse(p < 0.001, "***",
                ifelse(p < 0.01, "**",
                       ifelse(p < 0.05, "*", ifelse(include_ns, "ns", ""))))
  out[is.na(p)] <- if (include_ns) "ns" else ""
  out
}

build_stats <- function(catalog, value_mat, sample_meta, id_col, label_col, section_col) {
  mat <- as.matrix(value_mat)
  storage.mode(mat) <- "double"
  rownames(mat) <- catalog[[id_col]]
  alts <- sort(unique(sample_meta$Altitude))
  habitats <- c("R", "E")

  heat <- do.call(rbind, lapply(seq_len(nrow(catalog)), function(i) {
    pieces <- list()
    k <- 1L
    for (a in alts) {
      for (h in habitats) {
        idx <- sample_meta$Altitude == a & sample_meta$Habitat == h
        mu <- mean(mat[i, idx])
        pieces[[k]] <- data.frame(
          Feature_ID = catalog[[id_col]][i],
          Display_label = catalog[[label_col]][i],
          Section = catalog[[section_col]][i],
          Feature_order = catalog$Feature_order[i],
          Scope = "Altitude",
          Altitude = a,
          Habitat = h,
          Mean_TPM = mu,
          log2_mean_TPM_plus1 = log2(mu + 1),
          stringsAsFactors = FALSE
        )
        k <- k + 1L
      }
    }
    for (h in habitats) {
      idx <- sample_meta$Habitat == h
      mu <- mean(mat[i, idx])
      pieces[[k]] <- data.frame(
        Feature_ID = catalog[[id_col]][i],
        Display_label = catalog[[label_col]][i],
        Section = catalog[[section_col]][i],
        Feature_order = catalog$Feature_order[i],
        Scope = "Total",
        Altitude = NA_integer_,
        Habitat = h,
        Mean_TPM = mu,
        log2_mean_TPM_plus1 = log2(mu + 1),
        stringsAsFactors = FALSE
      )
      k <- k + 1L
    }
    do.call(rbind, pieces)
  }))

  altitude_tests <- do.call(rbind, lapply(seq_len(nrow(catalog)), function(i) {
    do.call(rbind, lapply(alts, function(a) {
      r_idx <- which(sample_meta$Altitude == a & sample_meta$Habitat == "R")
      e_idx <- which(sample_meta$Altitude == a & sample_meta$Habitat == "E")
      r_idx <- r_idx[order(sample_meta$Replicate[r_idx])]
      e_idx <- e_idx[order(sample_meta$Replicate[e_idx])]
      data.frame(
        Feature_ID = catalog[[id_col]][i],
        Display_label = catalog[[label_col]][i],
        Section = catalog[[section_col]][i],
        Altitude = a,
        n_pairs = length(r_idx),
        Mean_R = mean(mat[i, r_idx]),
        Mean_E = mean(mat[i, e_idx]),
        Higher = ifelse(mean(mat[i, r_idx]) >= mean(mat[i, e_idx]), "R", "E"),
        p_raw = safe_wilcox(mat[i, r_idx], mat[i, e_idx]),
        stringsAsFactors = FALSE
      )
    }))
  }))
  altitude_tests$p_BH_within_feature_5_altitudes <- ave(
    altitude_tests$p_raw, altitude_tests$Feature_ID,
    FUN = function(x) p.adjust(x, method = "BH")
  )
  altitude_tests$Significance <- star(altitude_tests$p_BH_within_feature_5_altitudes)

  total_tests <- do.call(rbind, lapply(seq_len(nrow(catalog)), function(i) {
    r_idx <- which(sample_meta$Habitat == "R")
    e_idx <- which(sample_meta$Habitat == "E")
    r_idx <- r_idx[order(sample_meta$Altitude[r_idx], sample_meta$Replicate[r_idx])]
    e_idx <- e_idx[order(sample_meta$Altitude[e_idx], sample_meta$Replicate[e_idx])]
    data.frame(
      Feature_ID = catalog[[id_col]][i],
      Display_label = catalog[[label_col]][i],
      Section = catalog[[section_col]][i],
      n_pairs = length(r_idx),
      Mean_R = mean(mat[i, r_idx]),
      Mean_E = mean(mat[i, e_idx]),
      Higher = ifelse(mean(mat[i, r_idx]) >= mean(mat[i, e_idx]), "R", "E"),
      p_raw = safe_wilcox(mat[i, r_idx], mat[i, e_idx]),
      stringsAsFactors = FALSE
    )
  }))
  total_tests$p_BH_all_displayed_rows <- p.adjust(total_tests$p_raw, method = "BH")
  total_tests$Significance <- star(total_tests$p_BH_all_displayed_rows)

  kw <- do.call(rbind, lapply(seq_len(nrow(catalog)), function(i) {
    do.call(rbind, lapply(habitats, function(h) {
      idx <- sample_meta$Habitat == h
      data.frame(
        Feature_ID = catalog[[id_col]][i],
        Display_label = catalog[[label_col]][i],
        Section = catalog[[section_col]][i],
        Habitat = h,
        p_raw = safe_kw(mat[i, idx], sample_meta$Altitude[idx]),
        stringsAsFactors = FALSE
      )
    }))
  }))
  kw$p_BH_all_displayed_rows_within_habitat <- ave(
    kw$p_raw, kw$Habitat, FUN = function(x) p.adjust(x, method = "BH")
  )
  kw$Significance <- star(kw$p_BH_all_displayed_rows_within_habitat, include_ns = TRUE)

  list(heat = heat, altitude = altitude_tests, total = total_tests, kw = kw)
}

write_element <- function(code, catalog, feature_tpm, sum_catalog, membership, note) {
  outdir <- file.path(work_root, code)
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  sample_cols <- sample_cols_for(feature_tpm)

  sample_meta <- data.frame(
    Sample = sample_cols,
    Habitat = ifelse(substr(sample_cols, 1, 1) == "S", "R", "E"),
    Altitude = as.integer(substr(sample_cols, 2, 3)) * 100,
    Replicate = as.integer(sub(".*_", "", sample_cols)),
    Pair_ID = paste(as.integer(substr(sample_cols, 2, 3)) * 100,
                    as.integer(sub(".*_", "", sample_cols)), sep = "_"),
    stringsAsFactors = FALSE
  )

  catalog <- catalog[order(catalog$Feature_order), , drop = FALSE]
  feature_tpm <- feature_tpm[match(catalog$Feature_ID, feature_tpm$Feature_ID), , drop = FALSE]
  stopifnot(identical(catalog$Feature_ID, feature_tpm$Feature_ID))
  gene_values <- feature_tpm[, sample_cols, drop = FALSE]
  gene_stats <- build_stats(catalog, gene_values, sample_meta,
                            "Feature_ID", "Display_label", "Section")

  stopifnot(all(membership$Feature_ID %in% catalog$Feature_ID))
  duplicate_within <- duplicated(membership[, c("Summary_ID", "Feature_ID")])
  stopifnot(!any(duplicate_within))
  sum_catalog <- sum_catalog[order(sum_catalog$Feature_order), , drop = FALSE]
  sum_values <- do.call(rbind, lapply(sum_catalog$Feature_ID, function(sid) {
    ids <- membership$Feature_ID[membership$Summary_ID == sid]
    stopifnot(length(ids) > 0L)
    colSums(gene_values[match(ids, catalog$Feature_ID), , drop = FALSE])
  }))
  rownames(sum_values) <- NULL
  colnames(sum_values) <- sample_cols
  sum_tpm <- cbind(sum_catalog, as.data.frame(sum_values, check.names = FALSE))
  sum_stats <- build_stats(sum_catalog, sum_values, sample_meta,
                           "Feature_ID", "Display_label", "Section")

  write.table(sample_meta, file.path(outdir, "sample_metadata.tsv"), sep = "\t",
              quote = FALSE, row.names = FALSE)
  write.table(catalog, file.path(outdir, paste0(code, "_gene_catalog.tsv")), sep = "\t",
              quote = FALSE, row.names = FALSE)
  write.table(feature_tpm, file.path(outdir, paste0(code, "_gene_TPM_by_sample.tsv")), sep = "\t",
              quote = FALSE, row.names = FALSE)
  write.table(gene_stats$heat, file.path(outdir, paste0(code, "_gene_heat_values.tsv")), sep = "\t",
              quote = FALSE, row.names = FALSE)
  write.table(gene_stats$altitude, file.path(outdir, paste0(code, "_gene_Wilcoxon_by_altitude.tsv")), sep = "\t",
              quote = FALSE, row.names = FALSE)
  write.table(gene_stats$total, file.path(outdir, paste0(code, "_gene_Wilcoxon_Total.tsv")), sep = "\t",
              quote = FALSE, row.names = FALSE)
  write.table(gene_stats$kw, file.path(outdir, paste0(code, "_gene_KW_altitude.tsv")), sep = "\t",
              quote = FALSE, row.names = FALSE)
  write.table(sum_catalog, file.path(outdir, paste0(code, "_sum_catalog.tsv")), sep = "\t",
              quote = FALSE, row.names = FALSE)
  write.table(membership, file.path(outdir, paste0(code, "_sum_membership.tsv")), sep = "\t",
              quote = FALSE, row.names = FALSE)
  write.table(sum_tpm, file.path(outdir, paste0(code, "_sum_TPM_by_sample.tsv")), sep = "\t",
              quote = FALSE, row.names = FALSE)
  write.table(sum_stats$heat, file.path(outdir, paste0(code, "_sum_heat_values.tsv")), sep = "\t",
              quote = FALSE, row.names = FALSE)
  write.table(sum_stats$altitude, file.path(outdir, paste0(code, "_sum_Wilcoxon_by_altitude.tsv")), sep = "\t",
              quote = FALSE, row.names = FALSE)
  write.table(sum_stats$total, file.path(outdir, paste0(code, "_sum_Wilcoxon_Total.tsv")), sep = "\t",
              quote = FALSE, row.names = FALSE)
  write.table(sum_stats$kw, file.path(outdir, paste0(code, "_sum_KW_altitude.tsv")), sep = "\t",
              quote = FALSE, row.names = FALSE)
  write.table(data.frame(
    Transform = "log2(mean raw TPM + 1)",
    Minimum = 0,
    Fixed_maximum = 15,
    Observed_gene_maximum = max(gene_stats$heat$log2_mean_TPM_plus1, na.rm = TRUE),
    Observed_sum_maximum = max(sum_stats$heat$log2_mean_TPM_plus1, na.rm = TRUE),
    Palette = "blue-white-red",
    Note = "Squares and circles share the established fixed CNPS scale"
  ), file.path(outdir, "shared_fill_scale.tsv"), sep = "\t",
  quote = FALSE, row.names = FALSE)
  writeLines(note, file.path(outdir, "membership_note.txt"))

  cat(sprintf("%s: %d gene rows, %d requested sum rows\n",
              code, nrow(catalog), nrow(sum_catalog)))
}

# Nitrogen: K00284 was already present as gltS. Rename only the display label to
# GLU; do not add a duplicate row or duplicate its abundance.
n_catalog <- read.delim(file.path(n_source, "N_feature_catalog.tsv"), check.names = FALSE)
n_tpm <- read.delim(file.path(n_source, "N_feature_TPM_by_sample.tsv"), check.names = FALSE)
n_catalog$Display_label[n_catalog$KO == "K00284"] <- "GLU"
n_tpm$Display_label[n_tpm$KO == "K00284"] <- "GLU"

n_sum_catalog <- data.frame(
  Feature_ID = sprintf("NS%02d", 1:5),
  Display_label = c(
    "NRT2 + nrtABCD",
    "nasAB + narB",
    "narGHI + napAB",
    "GDH family",
    "gltBD + GLU"
  ),
  Element = "Nitrogen",
  Section = "Requested N sums",
  Feature_order = 1:5,
  Source_type = "Within-sample raw TPM sum",
  stringsAsFactors = FALSE
)
n_members <- list(
  NS01 = c("N_NEW_K02575", "K15576", "K15577", "K15578", "K15579"),
  NS02 = c("K00360", "K00367", "K00372_nasA"),
  NS03 = c("K00370_narG", "K00371", "K00374", "K02567", "K02568"),
  NS04 = c("N_NEW_K00260", "K00261", "K00262", "K15371"),
  NS05 = c("K00265", "K00266", "K00284")
)
n_membership <- do.call(rbind, lapply(names(n_members), function(sid) {
  ids <- n_members[[sid]]
  data.frame(
    Summary_ID = sid,
    Summary_label = n_sum_catalog$Display_label[match(sid, n_sum_catalog$Feature_ID)],
    Feature_ID = ids,
    KO = n_catalog$KO[match(ids, n_catalog$Feature_ID)],
    Gene = n_catalog$Display_label[match(ids, n_catalog$Feature_ID)],
    Mapping_basis = c(
      NS01 = "User-defined NRT2 and nrtABCD transport group",
      NS02 = "User-defined nasAB and narB group; nasC and unresolved K00372 rows intentionally excluded",
      NS03 = "User-defined narGHI and napAB group; narZ and napC intentionally excluded",
      NS04 = "User-defined gudB/gdh family",
      NS05 = "User-defined gltBD and GLU group; GLT1/K00264 intentionally excluded"
    )[[sid]],
    stringsAsFactors = FALSE
  )
}))
write_element(
  "N", n_catalog, n_tpm, n_sum_catalog, n_membership,
  paste(
    "K00284 was already present as gltS and is displayed as GLU in this branch; it was not duplicated.",
    "NS02 follows the exact nasAB + narB definition: nasC and the unresolved K00372 row remain gene squares but are excluded from this circle.",
    "NS03 follows the requested narGHI + napAB definition: narZ and napC remain as gene squares but are excluded from this circle.",
    "NS05 follows gltBD + GLU: GLT1/K00264 remains a gene square but is excluded from this circle.",
    sep = "\n"
  )
)

# Sulfur: the user's "arsABC" is interpreted as the sulfur-cycle asrABC genes.
# arsABC would denote arsenic-resistance genes and is absent from the S catalog.
s_catalog <- read.delim(file.path(s_source, "S_feature_catalog.tsv"), check.names = FALSE)
s_tpm <- read.delim(file.path(s_source, "S_feature_TPM_by_sample.tsv"), check.names = FALSE)
s_sum_catalog <- data.frame(
  Feature_ID = sprintf("SS%02d", 1:6),
  Display_label = c(
    "PAPSS + sat + cysND",
    "asrABC",
    "sir + cysIJ",
    "cysEK",
    "fccAB",
    "doxD + tsdA"
  ),
  Element = "Sulfur",
  Section = "Requested S sums",
  Feature_order = 1:6,
  Source_type = "Within-sample raw TPM sum",
  stringsAsFactors = FALSE
)
s_members <- list(
  SS01 = c("S_NEW_K13811", "K00958", "K00956", "K00957"),
  SS02 = c("K16950", "K16951", "K00385"),
  SS03 = c("K00392", "K00380", "K00381"),
  SS04 = c("S_NEW_K00640", "S_NEW_K01738"),
  SS05 = c("K17230", "K17229"),
  SS06 = c("K16937", "K19713")
)
s_membership <- do.call(rbind, lapply(names(s_members), function(sid) {
  ids <- s_members[[sid]]
  data.frame(
    Summary_ID = sid,
    Summary_label = s_sum_catalog$Display_label[match(sid, s_sum_catalog$Feature_ID)],
    Feature_ID = ids,
    KO = s_catalog$KO[match(ids, s_catalog$Feature_ID)],
    Gene = s_catalog$Display_label[match(ids, s_catalog$Feature_ID)],
    Mapping_basis = if (sid == "SS02") {
      "User wrote arsABC; interpreted as sulfur-cycle asrABC (anaerobic sulfite reductase)"
    } else {
      "User-defined sulfur functional group"
    },
    stringsAsFactors = FALSE
  )
}))
write_element(
  "S", s_catalog, s_tpm, s_sum_catalog, s_membership,
  paste(
    "The requested 'arsABC' was interpreted as asrABC because asrABC is the sulfur-cycle anaerobic sulfite reductase family present in the list.",
    "The arsenic-resistance arsABC family is not present in the sulfur catalog and was not substituted.",
    "SS01 follows the exact requested four rows and does not include cysC or the bifunctional cysNC row.",
    sep = "\n"
  )
)

# Phosphorus: start from the current EggNOG-aware CNPS feature table. K01078 is
# the previously corrected PHO row and is relabelled accordingly.
base_all <- read.delim(base_feature_path, check.names = FALSE)
p_base <- base_all[base_all$Element == "Phosphorous", , drop = FALSE]
p_base <- p_base[order(p_base$List_order, p_base$Split_rank), , drop = FALSE]
stopifnot(nrow(p_base) == 45L)
# Keep the current main-text P scope: appA and phyC (phytase-related) and opd
# (the previously removed putative phosphotriesterase-related row) are excluded.
p_base <- p_base[!p_base$K_number %in% c("K01093", "K07048", "K01083"), , drop = FALSE]
stopifnot(nrow(p_base) == 42L)
p_catalog <- data.frame(
  Feature_ID = p_base$Feature_key,
  Display_label = p_base$Gene_name,
  Element = "Phosphorus",
  Section = p_base$Process,
  Feature_order = seq_len(nrow(p_base)),
  KO = p_base$K_number,
  Function_description = p_base$Function_description,
  Module_map = p_base$Module_or_Map_number,
  Source_type = p_base$Abundance_source,
  Source_url = "https://www.kegg.jp/pathway/map00440",
  stringsAsFactors = FALSE
)
p_catalog$Display_label[p_catalog$KO == "K01078"] <- "PHO"
p_values <- p_base[, sample_cols_for(p_base), drop = FALSE]
p_tpm <- cbind(p_catalog, p_values)

p_sum_catalog <- data.frame(
  Feature_ID = sprintf("PS%02d", 1:3),
  Display_label = c(
    "phoADN + aphA + PHO",
    "ugpABCEQ",
    "phoADN + aphA + PHO + ugpABCEQ"
  ),
  Element = "Phosphorus",
  Section = "Requested P sums",
  Feature_order = 1:3,
  Source_type = "Within-sample raw TPM sum",
  stringsAsFactors = FALSE
)
p_group1 <- c("K01077", "K01113", "K09474", "K03788", "K01078")
p_group2 <- c("K05814", "K05813", "K05816", "K05815", "K01126")
p_members <- list(
  PS01 = p_group1,
  PS02 = p_group2,
  PS03 = unique(c(p_group1, p_group2))
)
p_membership <- do.call(rbind, lapply(names(p_members), function(sid) {
  ids <- p_members[[sid]]
  data.frame(
    Summary_ID = sid,
    Summary_label = p_sum_catalog$Display_label[match(sid, p_sum_catalog$Feature_ID)],
    Feature_ID = ids,
    KO = p_catalog$KO[match(ids, p_catalog$Feature_ID)],
    Gene = p_catalog$Display_label[match(ids, p_catalog$Feature_ID)],
    Mapping_basis = if (sid == "PS03") {
      "Union of PS01 and PS02 member genes, summed once from sample-level TPM"
    } else {
      "User-defined phosphorus functional group"
    },
    stringsAsFactors = FALSE
  )
}))
write_element(
  "P", p_catalog, p_tpm, p_sum_catalog, p_membership,
  paste(
    "K01078 is displayed as PHO, following the user's corrected naming.",
    "appA/K01093, opd/K07048 and phyC/K01083 are excluded to match the current main-text P scope.",
    "PS03 is calculated from the unique union of PS01 and PS02 bottom-level rows; no pre-summed values are re-summed.",
    sep = "\n"
  )
)
