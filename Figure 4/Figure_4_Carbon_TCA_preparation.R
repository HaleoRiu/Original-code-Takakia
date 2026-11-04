options(stringsAsFactors = FALSE, warn = 1)

suppressPackageStartupMessages(library(readxl))

script_arg <- grep("^--file=", commandArgs(), value = TRUE)
script_dir <- if (length(script_arg)) dirname(normalizePath(sub("^--file=", "", script_arg[1]))) else getwd()
project_dir <- normalizePath(file.path(script_dir, ".."), mustWork = FALSE)
variant <- Sys.getenv("C_MAIN_VARIANT", unset = "V1_rTCA")
if (!variant %in% c("V1_rTCA", "V2_TCA")) stop("C_MAIN_VARIANT must be V1_rTCA or V2_TCA")
work <- file.path(project_dir, "work", variant)
dir.create(work, recursive = TRUE, showWarnings = FALSE)

catalog_path <- Sys.getenv(
  "C_CORE_CATALOG",
  unset = file.path(project_dir, "input", "Takakia_Carbon_maintext_rTCA_TCA_candidates_2026-07-22.xlsx")
)
tpm_path <- Sys.getenv(
  "C_CORE_TPM",
  unset = "/Volumes/Haleo_DK01/WD_Takakia_F/metagenome-data/tpm_renormalized.txt"
)

stopifnot(file.exists(catalog_path), file.exists(tpm_path))

label_sheet <- if (variant == "V1_rTCA") "Main_labels_V1_rTCA" else "Main_labels_V2_TCA"
ko_sheet <- if (variant == "V1_rTCA") "Main_KO_V1_rTCA" else "Main_KO_V2_TCA"
labels <- as.data.frame(read_excel(catalog_path, sheet = label_sheet), check.names = FALSE)
labels <- labels[!is.na(labels$Order), , drop = FALSE]
labels$Order <- as.integer(labels$Order)
labels$`KO count` <- as.integer(labels$`KO count`)
labels$`Detected KO count` <- as.integer(labels$`Detected KO count`)
labels <- labels[order(labels$Order), , drop = FALSE]

curated <- as.data.frame(read_excel(catalog_path, sheet = ko_sheet), check.names = FALSE)
curated <- curated[!is.na(curated$Order), , drop = FALSE]

expected_labels <- if (variant == "V1_rTCA") 39L else 44L
expected_kos <- if (variant == "V1_rTCA") 77L else 89L
stopifnot(nrow(labels) == expected_labels)
stopifnot(nrow(curated) == expected_kos)
stopifnot(identical(labels$Order, seq_len(nrow(labels))))
stopifnot(all(!is.na(labels$`Main group`)))

tpm <- read.delim(tpm_path, check.names = FALSE)
stopifnot("feature-id" %in% names(tpm))
sample_cols <- grep("^[PS](38|39|40|41|42)_[1-6]$", names(tpm), value = TRUE)
stopifnot(length(sample_cols) == 60L)

sample_meta <- data.frame(
  Sample = sample_cols,
  Habitat = ifelse(substr(sample_cols, 1, 1) == "S", "R", "E"),
  Altitude = as.integer(substr(sample_cols, 2, 3)) * 100,
  Replicate = as.integer(sub(".*_", "", sample_cols)),
  stringsAsFactors = FALSE
)
sample_meta$Pair_ID <- paste(sample_meta$Altitude, sample_meta$Replicate, sep = "_")
write.table(sample_meta, file.path(work, "sample_metadata.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

split_kos <- function(x) {
  z <- trimws(unlist(strsplit(as.character(x), ";", fixed = TRUE)))
  unique(z[nzchar(z)])
}

label_kos <- lapply(labels$`Included KOs`, split_kos)
stopifnot(all(lengths(label_kos) == labels$`KO count`))

membership <- do.call(rbind, lapply(seq_len(nrow(labels)), function(i) {
  data.frame(
    Feature_ID = sprintf("C%02d", labels$Order[i]),
    Feature_order = labels$Order[i],
    Carbon_domain = labels$`Carbon domain`[i],
    Section = labels$`Main group`[i],
    Display_label = labels$`Functional label`[i],
    Priority = labels$Priority[i],
    Figure_tier = labels$`Figure tier`[i],
    KO = label_kos[[i]],
    stringsAsFactors = FALSE
  )
}))
membership$Detected_in_TPM <- membership$KO %in% tpm$`feature-id`

detected_by_label <- aggregate(Detected_in_TPM ~ Feature_ID, membership, sum)
detected_check <- merge(
  data.frame(Feature_ID = sprintf("C%02d", labels$Order), Expected = labels$`Detected KO count`),
  transform(detected_by_label, Calculated = as.integer(Detected_in_TPM))[, c("Feature_ID", "Calculated")],
  by = "Feature_ID", all = TRUE
)
detected_check$Match <- detected_check$Expected == detected_check$Calculated
write.table(detected_check, file.path(work, "detected_KO_count_check.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
if (!all(detected_check$Match)) stop("Detected-KO counts disagree with the workbook; see detected_KO_count_check.tsv")

dup_across_labels <- aggregate(Feature_ID ~ KO, membership, function(x) paste(unique(x), collapse = ";"))
dup_across_labels$N_labels <- vapply(strsplit(dup_across_labels$Feature_ID, ";", fixed = TRUE), length, integer(1))
write.table(dup_across_labels[dup_across_labels$N_labels > 1, , drop = FALSE],
            file.path(work, "KO_membership_duplicates_across_labels.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

tpm_index <- match(membership$KO, tpm$`feature-id`)
ko_mat <- matrix(0, nrow = nrow(membership), ncol = length(sample_cols),
                 dimnames = list(membership$KO, sample_cols))
present <- !is.na(tpm_index)
ko_mat[present, ] <- as.matrix(tpm[tpm_index[present], sample_cols, drop = FALSE])
storage.mode(ko_mat) <- "double"

sum_member_rows <- function(row_index) {
  colSums(ko_mat[row_index, , drop = FALSE])
}

label_rows <- lapply(seq_len(nrow(labels)), function(i) {
  fid <- sprintf("C%02d", labels$Order[i])
  idx <- which(membership$Feature_ID == fid)
  data.frame(
    Feature_ID = fid,
    Display_label = labels$`Functional label`[i],
    Element = "Carbon",
    Carbon_domain = labels$`Carbon domain`[i],
    Section = labels$`Main group`[i],
    Feature_order = labels$Order[i],
    Priority = labels$Priority[i],
    Figure_tier = labels$`Figure tier`[i],
    KO_count = length(idx),
    Detected_KO_count = sum(membership$Detected_in_TPM[idx]),
    as.list(sum_member_rows(idx)),
    check.names = FALSE
  )
})
label_tpm <- do.call(rbind, label_rows)

section_order <- unique(labels$`Main group`)

section_rows <- lapply(seq_along(section_order), function(i) {
  sec <- section_order[i]
  idx <- which(membership$Section == sec)
  unique_ko <- unique(membership$KO[idx])
  unique_idx <- match(unique_ko, membership$KO)
  domains <- unique(membership$Carbon_domain[idx])
  data.frame(
    Feature_ID = sprintf("CS%02d", i),
    Display_label = sec,
    Element = "Carbon",
    Carbon_domain = domains[1],
    Section = sec,
    Feature_order = i,
    Priority = "Section sum (all priorities)",
    Figure_tier = "Trend draft",
    KO_count = length(unique_ko),
    Detected_KO_count = sum(unique_ko %in% tpm$`feature-id`),
    as.list(colSums(ko_mat[unique_idx, , drop = FALSE])),
    check.names = FALSE
  )
})
section_tpm <- do.call(rbind, section_rows)

section_membership <- do.call(rbind, lapply(seq_along(section_order), function(i) {
  sec <- section_order[i]
  idx <- which(membership$Section == sec)
  data.frame(
    Summary_ID = sprintf("CS%02d", i),
    Summary_label = sec,
    Summary_order = i,
    Section = sec,
    KO = unique(membership$KO[idx]),
    stringsAsFactors = FALSE
  )
}))

write.table(membership, file.path(work, "C_label_KO_membership.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(section_membership, file.path(work, "C_section_KO_membership.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(label_tpm, file.path(work, "C_label_summed_TPM_by_sample.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(section_tpm, file.path(work, "C_section_summed_TPM_by_sample.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(label_tpm[, 1:10], file.path(work, "C_label_catalog.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(section_tpm[, 1:10], file.path(work, "C_section_catalog.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

safe_wilcox <- function(x, y) {
  ok <- is.finite(x) & is.finite(y)
  x <- x[ok]; y <- y[ok]
  if (length(x) < 2L || all(x == y)) return(1)
  out <- tryCatch(suppressWarnings(wilcox.test(x, y, paired = TRUE, exact = FALSE)$p.value),
                  error = function(e) NA_real_)
  ifelse(is.finite(out), out, 1)
}

safe_kw <- function(x, g) {
  ok <- is.finite(x) & !is.na(g)
  x <- x[ok]; g <- g[ok]
  if (length(unique(g)) < 2L || length(unique(x)) < 2L) return(1)
  out <- tryCatch(kruskal.test(x ~ factor(g))$p.value, error = function(e) NA_real_)
  ifelse(is.finite(out), out, 1)
}

star <- function(p, include_ns = FALSE) {
  ans <- ifelse(p < 0.001, "***",
                ifelse(p < 0.01, "**",
                       ifelse(p < 0.05, "*", ifelse(include_ns, "ns", ""))))
  ans[is.na(p)] <- if (include_ns) "ns" else ""
  ans
}

compute_stats <- function(dat, block_name) {
  meta <- dat[, 1:10]
  mat <- as.matrix(dat[, sample_cols, drop = FALSE])
  storage.mode(mat) <- "double"
  alts <- sort(unique(sample_meta$Altitude))
  habitats <- c("R", "E")

  heat <- do.call(rbind, lapply(seq_len(nrow(dat)), function(i) {
    pieces <- list(); k <- 1L
    for (a in alts) for (h in habitats) {
      idx <- sample_meta$Altitude == a & sample_meta$Habitat == h
      mu <- mean(mat[i, idx])
      pieces[[k]] <- data.frame(
        Feature_ID = meta$Feature_ID[i], Display_label = meta$Display_label[i],
        Element = "Carbon", Carbon_domain = meta$Carbon_domain[i], Section = meta$Section[i],
        Feature_order = meta$Feature_order[i], Block = block_name,
        Scope = "Altitude", Altitude = a, Habitat = h,
        Mean_TPM = mu, log2_mean_TPM_plus1 = log2(mu + 1)
      ); k <- k + 1L
    }
    for (h in habitats) {
      idx <- sample_meta$Habitat == h
      mu <- mean(mat[i, idx])
      pieces[[k]] <- data.frame(
        Feature_ID = meta$Feature_ID[i], Display_label = meta$Display_label[i],
        Element = "Carbon", Carbon_domain = meta$Carbon_domain[i], Section = meta$Section[i],
        Feature_order = meta$Feature_order[i], Block = block_name,
        Scope = "Total", Altitude = NA_integer_, Habitat = h,
        Mean_TPM = mu, log2_mean_TPM_plus1 = log2(mu + 1)
      ); k <- k + 1L
    }
    do.call(rbind, pieces)
  }))

  altitude_tests <- do.call(rbind, lapply(seq_len(nrow(dat)), function(i) {
    do.call(rbind, lapply(alts, function(a) {
      r_idx <- which(sample_meta$Altitude == a & sample_meta$Habitat == "R")
      e_idx <- which(sample_meta$Altitude == a & sample_meta$Habitat == "E")
      r_idx <- r_idx[order(sample_meta$Replicate[r_idx])]
      e_idx <- e_idx[order(sample_meta$Replicate[e_idx])]
      data.frame(
        Feature_ID = meta$Feature_ID[i], Display_label = meta$Display_label[i],
        Carbon_domain = meta$Carbon_domain[i], Section = meta$Section[i],
        Block = block_name, Altitude = a, n_pairs = length(r_idx),
        Mean_R = mean(mat[i, r_idx]), Mean_E = mean(mat[i, e_idx]),
        Higher = ifelse(mean(mat[i, r_idx]) >= mean(mat[i, e_idx]), "R", "E"),
        p_raw = safe_wilcox(mat[i, r_idx], mat[i, e_idx])
      )
    }))
  }))
  altitude_tests$p_BH_within_feature_5_altitudes <- ave(
    altitude_tests$p_raw, altitude_tests$Feature_ID,
    FUN = function(x) p.adjust(x, method = "BH")
  )
  altitude_tests$p_BH_audit_within_altitude_all_features <- ave(
    altitude_tests$p_raw, altitude_tests$Altitude,
    FUN = function(x) p.adjust(x, method = "BH")
  )
  altitude_tests$Significance <- star(altitude_tests$p_BH_within_feature_5_altitudes)

  total_tests <- do.call(rbind, lapply(seq_len(nrow(dat)), function(i) {
    r_idx <- which(sample_meta$Habitat == "R")
    e_idx <- which(sample_meta$Habitat == "E")
    r_idx <- r_idx[order(sample_meta$Altitude[r_idx], sample_meta$Replicate[r_idx])]
    e_idx <- e_idx[order(sample_meta$Altitude[e_idx], sample_meta$Replicate[e_idx])]
    data.frame(
      Feature_ID = meta$Feature_ID[i], Display_label = meta$Display_label[i],
      Carbon_domain = meta$Carbon_domain[i], Section = meta$Section[i],
      Block = block_name, n_pairs = length(r_idx),
      Mean_R = mean(mat[i, r_idx]), Mean_E = mean(mat[i, e_idx]),
      Higher = ifelse(mean(mat[i, r_idx]) >= mean(mat[i, e_idx]), "R", "E"),
      p_raw = safe_wilcox(mat[i, r_idx], mat[i, e_idx])
    )
  }))
  total_tests$p_BH_all_features_in_block <- p.adjust(total_tests$p_raw, method = "BH")
  total_tests$Significance <- star(total_tests$p_BH_all_features_in_block)

  kw <- do.call(rbind, lapply(seq_len(nrow(dat)), function(i) {
    do.call(rbind, lapply(habitats, function(h) {
      idx <- sample_meta$Habitat == h
      data.frame(
        Feature_ID = meta$Feature_ID[i], Display_label = meta$Display_label[i],
        Carbon_domain = meta$Carbon_domain[i], Section = meta$Section[i],
        Block = block_name, Habitat = h,
        p_raw = safe_kw(mat[i, idx], sample_meta$Altitude[idx])
      )
    }))
  }))
  kw$p_BH_all_features_within_habitat <- ave(
    kw$p_raw, kw$Habitat, FUN = function(x) p.adjust(x, method = "BH")
  )
  kw$Significance <- star(kw$p_BH_all_features_within_habitat, include_ns = TRUE)

  list(heat = heat, altitude = altitude_tests, total = total_tests, kw = kw)
}

label_stats <- compute_stats(label_tpm, "Functional label")
section_stats <- compute_stats(section_tpm, "Section sum")

write.table(label_stats$heat, file.path(work, "C_label_heat_values.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(label_stats$altitude, file.path(work, "C_label_paired_Wilcoxon_by_altitude.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(label_stats$total, file.path(work, "C_label_paired_Wilcoxon_Total.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(label_stats$kw, file.path(work, "C_label_KW_altitude.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(section_stats$heat, file.path(work, "C_section_heat_values.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(section_stats$altitude, file.path(work, "C_section_paired_Wilcoxon_by_altitude.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(section_stats$total, file.path(work, "C_section_paired_Wilcoxon_Total.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(section_stats$kw, file.path(work, "C_section_KW_altitude.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

all_heat <- rbind(label_stats$heat, section_stats$heat)
observed_max <- max(all_heat$log2_mean_TPM_plus1, na.rm = TRUE)
write.table(data.frame(
  Transform = "log2(mean raw summed TPM + 1)",
  Minimum = 0,
  Fixed_maximum = 15,
  Observed_maximum = observed_max,
  Palette = "blue-white-red",
  Note = "Fixed to the previous complete CNPS altitude figure for direct visual comparability"
), file.path(work, "shared_fill_scale.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

priority_counts <- as.data.frame(table(labels$Priority), stringsAsFactors = FALSE)
names(priority_counts) <- c("Priority", "N_labels")
write.table(priority_counts, file.path(work, "included_priority_counts.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

cat(sprintf("Variant %s; included labels: %d; member KO rows: %d; main groups: %d\n",
            variant, nrow(label_tpm), nrow(membership), nrow(section_tpm)))
cat("Included labels by priority:\n")
print(priority_counts, row.names = FALSE)
cat(sprintf("Detected unique member KOs: %d of %d\n",
            sum(unique(membership[, c("KO", "Detected_in_TPM")])$Detected_in_TPM),
            length(unique(membership$KO))))
cat(sprintf("Observed log2(mean TPM + 1) maximum: %.4f; plotted scale fixed at 0-15\n", observed_max))
cat(sprintf("Label altitude significant cells: %d; Total significant labels: %d; KW significant label-habitat tests: %d\n",
            sum(nzchar(label_stats$altitude$Significance)),
            sum(nzchar(label_stats$total$Significance)),
            sum(label_stats$kw$Significance != "ns")))
