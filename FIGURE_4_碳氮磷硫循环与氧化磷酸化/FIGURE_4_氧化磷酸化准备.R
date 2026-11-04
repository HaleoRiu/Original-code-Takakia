options(stringsAsFactors = FALSE, warn = 1)

work <- "/Users/haleoriu/Documents/Takakia绘图和文章修改/outputs/CNPS_altitude_sum_revision_2026-07-21/V2_3_categories/work"
feature_tpm_path <- "/Volumes/Haleo_DK01/WD_Takakia_F/metagenome-data/Figure7_CNPS/CNPS_EggNOG_split_results/prepared/CNPS_feature_renormalized_TPM_by_sample.tsv"
c_snapshot_path <- "/Volumes/Haleo_DK01/WD_Takakia_F/metagenome-data/Figure7_CNPS/CNPS_altitude_pathway_2026-07-20/tables/latest_CNPS_feature_list_snapshot.tsv"
s_transport_path <- "/Volumes/Haleo_DK01/WD_Takakia_F/metagenome-data/Figure7_CNPS/S_transporter_screen_2026-07-20/retained_KO_renormalized_TPM_by_sample.tsv"
c_rule_path <- "/Volumes/Haleo_DK01/WD_Takakia_F/metagenome-data/Figure7_CNPS/CNPS_main_text_C_NC_grouped_results/raw_tables/Carbon_enzyme_grouping_rules.tsv"
c_tpm_path <- "/Volumes/Haleo_DK01/WD_Takakia_F/metagenome-data/Figure7_CNPS/CNPS_main_text_C_NC_grouped_results/raw_tables/Carbon_NC_grouped_TPM_by_sample.tsv"
op_path_tpm_path <- file.path(work, "oxidative_pathway_summed_TPM_by_sample.tsv")
op_map_path <- file.path(work, "oxidative_pathway_KO_mapping.tsv")

feature_tpm <- read.delim(feature_tpm_path, check.names = FALSE)
c_rules <- read.delim(c_rule_path, check.names = FALSE)
c_tpm <- read.delim(c_tpm_path, check.names = FALSE)
op_tpm <- read.delim(op_path_tpm_path, check.names = FALSE)
op_map <- read.delim(op_map_path, check.names = FALSE)
s_transport_tpm <- read.delim(s_transport_path, check.names = FALSE)

sample_cols <- grep("^[PS](38|39|40|41|42)_[1-6]$", names(feature_tpm), value = TRUE)
stopifnot(length(sample_cols) == 60L)
stopifnot(identical(sample_cols, grep("^[PS](38|39|40|41|42)_[1-6]$", names(c_tpm), value = TRUE)))

sample_meta <- data.frame(
  Sample = sample_cols,
  Habitat = ifelse(substr(sample_cols, 1, 1) == "S", "R", "E"),
  Altitude = as.integer(substr(sample_cols, 2, 3)) * 100,
  Replicate = as.integer(sub(".*_", "", sample_cols)),
  stringsAsFactors = FALSE
)
sample_meta$Pair_ID <- paste(sample_meta$Altitude, sample_meta$Replicate, sep = "_")
write.table(sample_meta, file.path(work, "sample_metadata.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

make_catalog <- function(id, label, element, section, order, tpm, source_type) {
  out <- data.frame(
    Feature_ID = id,
    Display_label = label,
    Element = element,
    Section = section,
    Feature_order = order,
    Source_type = source_type,
    stringsAsFactors = FALSE
  )
  cbind(out, tpm[, sample_cols, drop = FALSE])
}

# Upper block: Carbon is the already-reviewed 39-enzyme aggregation.
c_rules <- c_rules[order(c_rules$Panel_order), ]
c_tpm <- c_tpm[match(c_rules$Feature_set_ID, c_tpm$Feature_set_ID), ]
c_upper <- make_catalog(c_rules$Feature_set_ID, c_rules$Enzyme_label, "Carbon",
                        c_rules$Section, c_rules$Panel_order, c_tpm,
                        "Integrated enzyme group")

# Upper block: N, P and S retain the latest curated feature rows.
make_element_upper <- function(element_name, display_element = element_name) {
  d <- feature_tpm[feature_tpm$Element == element_name, ]
  d <- d[order(d$List_order, d$Split_rank), ]
  make_catalog(d$Feature_key, d$Gene_name, display_element, d$Process,
               seq_len(nrow(d)), d, "Curated CNPS feature")
}
n_upper <- make_element_upper("Nitrogen")
p_upper <- make_element_upper("Phosphorous", "Phosphorus")
s_upper <- make_element_upper("Sulfur")

# Append the 12 detected literature-backed sulfur-transporter KOs to the end of
# the sulfur upper block. Their abundance values are direct rows from
# tpm_renormalized.txt, as recorded in the transporter screening audit.
s_transport_meta <- data.frame(
  KO = c("K02045", "K02046", "K02047", "K02048", "K06203", "K03321",
         "K15553", "K15554", "K15555", "K15551", "K15552", "K10831"),
  Gene = c("cysA", "cysU/cysT", "cysW", "cysP/sbp", "cysZ", "sulP",
           "ssuA", "ssuC", "ssuB", "tauA", "tauC", "tauB"),
  Transport_family = c(rep("Inorganic sulfate/thiosulfate transport", 6),
                       rep("Alkanesulfonate ABC transporter", 3),
                       rep("Taurine ABC transporter", 3)),
  Function_description = c(
    "sulfate/thiosulfate ABC transporter ATP-binding protein",
    "sulfate/thiosulfate ABC transporter permease",
    "sulfate/thiosulfate ABC transporter permease",
    "thiosulfate/sulfate substrate-binding protein",
    "CysZ sulfate transporter", "SulP sulfate permease",
    "alkanesulfonate-binding protein", "alkanesulfonate transporter permease",
    "alkanesulfonate transporter ATP-binding protein", "taurine-binding protein",
    "taurine transporter permease", "taurine transporter ATP-binding protein"
  ),
  stringsAsFactors = FALSE
)
stopifnot(identical(s_transport_tpm$KO, s_transport_meta$KO))
stopifnot(identical(s_transport_tpm$Gene, s_transport_meta$Gene))
stopifnot(identical(sample_cols, names(s_transport_tpm)[-(1:2)]))

s_transport_upper <- make_catalog(
  paste0("S_TRANS_", s_transport_meta$KO), s_transport_meta$Gene, "Sulfur",
  s_transport_meta$Transport_family, nrow(s_upper) + seq_len(nrow(s_transport_meta)),
  s_transport_tpm, "Direct retained KO row from tpm_renormalized.txt"
)
s_upper <- rbind(s_upper, s_transport_upper)

s_transport_catalog <- cbind(
  Feature_ID = s_transport_upper$Feature_ID,
  s_transport_meta,
  Feature_order = s_transport_upper$Feature_order,
  Abundance_source = "Direct row from tpm_renormalized.txt",
  Literature_basis = "SCycDB transporter families; Liu et al. 2016 Table S8",
  stringsAsFactors = FALSE
)
write.table(s_transport_catalog, file.path(work, "S_transport_KO_catalog.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

# Preserve the previous curated snapshot and append the new transporter rows for
# a complete reviewable list without modifying the source workbook.
base_snapshot <- read.delim(c_snapshot_path, check.names = FALSE)
new_snapshot <- data.frame(
  Element = "Sulfur", Sheet_order = 4L,
  List_order = nrow(base_snapshot[base_snapshot$Element == "Sulfur", ]) + seq_len(nrow(s_transport_meta)),
  Process = s_transport_meta$Transport_family, Gene_name = s_transport_meta$Gene,
  K_number = s_transport_meta$KO, Function_description = s_transport_meta$Function_description,
  Module_or_Map_number = "map00920; map02010",
  Reference = "SCycDB; Liu et al. 2016 Supplementary Table S8",
  Original_List_order = NA_integer_, K_primary = s_transport_meta$KO,
  Feature_key = s_transport_upper$Feature_ID, Split_rank = 1L,
  Abundance_source = "Direct row from tpm_renormalized.txt",
  Feature_ID = max(base_snapshot$Feature_ID, na.rm = TRUE) + seq_len(nrow(s_transport_meta)),
  stringsAsFactors = FALSE
)
new_snapshot <- new_snapshot[, names(base_snapshot), drop = FALSE]
write.table(rbind(base_snapshot, new_snapshot), file.path(work, "latest_CNPS_feature_list_snapshot.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

# Upper block: every supplied oxidative-phosphorylation/associated Pathway is retained,
# but abundance is recalculated as the direct sum of detected unique member-KO TPMs.
op_class <- ifelse(op_tpm$Pathway_order <= 10,
                   "Core oxidative phosphorylation",
                   "Other supplied pathways")
op_upper <- make_catalog(sprintf("OP%02d", op_tpm$Pathway_order), op_tpm$Pathway,
                         "Oxidative phosphorylation", op_class,
                         op_tpm$Pathway_order, op_tpm,
                         "Detected-unique-KO direct TPM sum")

upper <- rbind(c_upper, n_upper, p_upper, s_upper, op_upper)
write.table(upper, file.path(work, "upper_feature_TPM_by_sample.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(upper[, 1:6], file.path(work, "upper_feature_catalog.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

# Pathway-summary membership for N/P/S. Entries not explicitly grouped in list-cnps remain singleton pathways.
init_membership <- function(d, prefix) {
  data.frame(
    Element = d$Element,
    Feature_ID = d$Feature_ID,
    Display_label = d$Display_label,
    Feature_order = d$Feature_order,
    Summary_ID = paste0(prefix, "_single_", d$Feature_ID),
    Summary_label = d$Display_label,
    Summary_order = 1000 + d$Feature_order,
    Mapping_basis = "Unlisted singleton retained",
    stringsAsFactors = FALSE
  )
}

assign_group <- function(m, ids, summary_id, summary_label, summary_order, basis) {
  hit <- m$Feature_ID %in% ids
  m$Summary_ID[hit] <- summary_id
  m$Summary_label[hit] <- summary_label
  m$Summary_order[hit] <- summary_order
  m$Mapping_basis[hit] <- basis
  m
}

# Nitrogen groups. Correct list-cnps ntrABCD -> the curated nrtABCD names.
n_map <- init_membership(n_upper, "N")
n_ids <- setNames(n_upper$Feature_ID, n_upper$Display_label)
ids_by_label <- function(d, labels) d$Feature_ID[d$Display_label %in% labels]
n_map <- assign_group(n_map, ids_by_label(n_upper, c("nifA","nifS","nifD","nifH","nifK")), "N01", "Nitrogen fixation", 1, "list-cnps: nifASDHK")
n_map <- assign_group(n_map, ids_by_label(n_upper, "hao"), "N02", "Ammonia oxidation", 2, "Unlisted singleton retained")
n_map <- assign_group(n_map, ids_by_label(n_upper, c("nasA","nasB","nasC","nasA/nasC unresolved","narB")), "N03", "Assimilatory nitrate reduction", 3, "nasABC (including unresolved nasA/nasC) plus narB")
n_map <- assign_group(n_map, ids_by_label(n_upper, c("narG","narZ","narH","narI","napA","napB","napC")), "N04", "Dissimilatory nitrate reduction", 4, "narGZHI plus napABC")
n_map <- assign_group(n_map, ids_by_label(n_upper, "nirK"), "N05", "Nitrite reduction (nirK)", 5, "Unlisted singleton retained")
n_map <- assign_group(n_map, ids_by_label(n_upper, c("norB","norC")), "N06", "Nitric oxide reduction", 6, "list-cnps: norBC")
n_map <- assign_group(n_map, ids_by_label(n_upper, "nosZ"), "N07", "Nitrous oxide reduction", 7, "Unlisted singleton retained")
n_map <- assign_group(n_map, ids_by_label(n_upper, c("nrfA","nrfH","nirB","nirD")), "N08", "DNRA nitrite reduction", 8, "list-cnps: nrfA, nrfH, nirB, nirD")
n_map <- assign_group(n_map, ids_by_label(n_upper, "nirA"), "N09", "Assimilatory nitrite reduction", 9, "Unlisted singleton retained")
n_map <- assign_group(n_map, ids_by_label(n_upper, c("nrtA","nrtB","nrtC","nrtD")), "N10", "Nitrate/nitrite transport", 10, "list-cnps ntrABCD interpreted as curated nrtABCD")
n_map <- assign_group(n_map, ids_by_label(n_upper, c("urtA","urtB","urtC","urtD","urtE")), "N11", "Urea transport", 11, "list-cnps: urtABCDE")
n_map <- assign_group(n_map, ids_by_label(n_upper, c("ureA","ureB","ureC")),
                      "N12", "Urea hydrolysis (ureABC)", 12,
                      "User revision 2026-07-21: ureABC summed separately")
n_map <- assign_group(n_map, ids_by_label(n_upper, c("gdhA","gdh","gltB","gltD","gltS","glnA","glsA")),
                      "N13", "N assimilation (gdh/glt/gln/gls)", 13,
                      "User revision 2026-07-21: 8 gdh/glt/gln/gls feature rows summed together")

# Phosphorus groups. Correct list-cnps piy -> curated pit. The 2026-07-21
# revision separates pst+pit, ugp, phnAWX, the remaining phn operon, and an
# overlapping all-phn family total.
p_map <- init_membership(p_upper, "P")
p_map <- assign_group(p_map, ids_by_label(p_upper, c("pstA","pstB","pstC","pstS","pit")),
                      "P01", "Inorganic phosphate uptake (pst/pit)", 1,
                      "User revision 2026-07-21: pstABCS plus pit")
p_map <- assign_group(p_map, ids_by_label(p_upper, c("ugpA","ugpB","ugpC","ugpE")),
                      "P02", "Glycerol-3-phosphate uptake (ugp)", 2,
                      "User revision 2026-07-21: ugpABCE")
p_map <- assign_group(p_map, ids_by_label(p_upper, c("phoB","phoP","phoR","phoU")),
                      "P03", "Phosphate-starvation regulation", 3, "list-cnps: phoBPRU")
p_map <- assign_group(p_map, ids_by_label(p_upper, c("phoA","phoD","phoN","appA","opd","ugpQ","aphA","olpA","phyC")),
                      "P04", "Other organic-P mineralisation", 4,
                      "User revision 2026-07-21: previous other-organic-P group after removing phnAWX")
p_map <- assign_group(p_map, ids_by_label(p_upper, c("phnA","phnW","phnX")),
                      "P05", "phnAWX", 5,
                      "User revision 2026-07-21: phnAWX summed separately")
p_phn_operon_labels <- c("phnC","phnD","phnE","phnF","phnN","phnL","phnK","phnO","phnM","phnJ","phnI","phnH","phnG","phnP")
p_map <- assign_group(p_map, ids_by_label(p_upper, p_phn_operon_labels),
                      "P06", "phn operon", 6,
                      "User revision 2026-07-21: all remaining phn-operon feature rows")
p_map <- assign_group(p_map, ids_by_label(p_upper, c("gcd","pqqB","pqqC")),
                      "P08", "Inorganic-P solubilisation", 8, "list-cnps: gcd, pqqBC")
p_map <- assign_group(p_map, ids_by_label(p_upper, "ppk1"),
                      "P09", "Polyphosphate synthesis", 9, "Unlisted singleton retained")
p_map <- assign_group(p_map, ids_by_label(p_upper, c("ppa","ppx")),
                      "P10", "Polyphosphate degradation", 10, "list-cnps: ppa, ppx")

# Overlapping total: phnAWX plus every remaining phn-operon row.
p_phn_family_membership <- p_map[p_map$Summary_ID %in% c("P05", "P06"), ]
p_phn_family_membership$Summary_ID <- "P07"
p_phn_family_membership$Summary_label <- "phn family"
p_phn_family_membership$Summary_order <- 7
p_phn_family_membership$Mapping_basis <- "User revision 2026-07-21: phnAWX plus all remaining phn-operon rows"
p_membership <- rbind(p_map, p_phn_family_membership)

# Sulfur groups. Correct list-cnps arsABC -> curated asrABC.
s_map <- init_membership(s_upper, "S")
s_map <- assign_group(s_map, ids_by_label(s_upper, c("sorA","soeA","soeB","soeC")), "S01", "Sulfite oxidation", 1, "list-cnps: sorA, soeABC")
s_map <- assign_group(s_map, ids_by_label(s_upper, c("sat","cysC","cysD","cysNC","cysN","cysH","cysJ","cysI")), "S02", "Assimilatory sulfate reduction", 2, "list-cnps: sat, cysCDHIJN, cysNC")
s_map <- assign_group(s_map, ids_by_label(s_upper, c("aprA","aprB")), "S03", "Dissimilatory sulfate activation", 3, "list-cnps: aprAB")
s_map <- assign_group(s_map, ids_by_label(s_upper, c("asrA","asrB","asrC","sir")), "S04", "Sulfite reduction", 4, "list-cnps arsABC interpreted as curated asrABC; plus sir")
s_map <- assign_group(s_map, ids_by_label(s_upper, "phsC"), "S05", "Thiosulfate reduction", 5, "Unlisted singleton retained")
s_map <- assign_group(s_map, ids_by_label(s_upper, c("sqr","fccA","fccB")), "S06", "Sulfide oxidation", 6, "list-cnps: fccAB, sqr")
s_map <- assign_group(s_map, ids_by_label(s_upper, c("soxA","soxX","soxB","soxC","soxY","soxZ")), "S07", "SOX thiosulfate oxidation", 7, "list-cnps: soxABCXYZ")
s_map <- assign_group(s_map, ids_by_label(s_upper, c("doxD","tsdA","ttrA")), "S08", "Tetrathionate transformation", 8, "list-cnps: ttrA, doxD, tsdA")
s_map <- assign_group(s_map, ids_by_label(s_upper, c("dmsA","dmsB","dmsC")), "S09", "DMSO reduction", 9, "list-cnps: dmsABC; both curated dmsBC KO pairs retained")
s_map <- assign_group(s_map, ids_by_label(s_upper, "mddA"), "S10", "Methanethiol formation", 10, "Unlisted singleton retained")
s_map <- assign_group(s_map, ids_by_label(s_upper, "ddhC"), "S11", "DMS dehydrogenation", 11, "Unlisted singleton retained")
s_map <- assign_group(s_map, ids_by_label(s_upper, "tmm"), "S12", "Trimethylamine monooxygenation", 12, "Unlisted singleton retained")
s_map <- assign_group(s_map, paste0("S_TRANS_", c("K02045","K02046","K02047","K02048","K06203","K03321")),
                      "S13", "Inorganic sulfate/thiosulfate transport", 13,
                      "Combined category: cysAPUW-type ABC + CysZ + SulP")
s_map <- assign_group(s_map, paste0("S_TRANS_", c("K15553","K15554","K15555")),
                      "S14", "Alkanesulfonate transport", 14,
                      "Combined category: ssuA + ssuC + ssuB")
s_map <- assign_group(s_map, paste0("S_TRANS_", c("K15551","K15552","K10831")),
                      "S15", "Taurine transport", 15,
                      "Combined category: tauA + tauC + tauB")

# Carbon: the 39 enzyme rows collapse to their existing 12 reviewed Sections.
c_map <- data.frame(
  Element = c_upper$Element,
  Feature_ID = c_upper$Feature_ID,
  Display_label = c_upper$Display_label,
  Feature_order = c_upper$Feature_order,
  Summary_ID = sprintf("CP%02d", c_rules$Section_order),
  Summary_label = c_upper$Section,
  Summary_order = c_rules$Section_order,
  Mapping_basis = "Existing Carbon 39-enzyme Section",
  stringsAsFactors = FALSE
)

# Oxidative block is already Pathway-level; reproduce the direct KO-sum Pathways as circle summaries.
op_members <- aggregate(KO ~ Pathway_order + Pathway, data = op_map[!op_map$Is_duplicate_within_pathway, ],
                        FUN = function(x) paste(x, collapse = ";"))
op_map_summary <- data.frame(
  Element = op_upper$Element,
  Feature_ID = op_upper$Feature_ID,
  Display_label = op_upper$Display_label,
  Feature_order = op_upper$Feature_order,
  Summary_ID = sprintf("OPP%02d", op_upper$Feature_order),
  Summary_label = op_upper$Display_label,
  Summary_order = op_upper$Feature_order,
  Mapping_basis = "Detected-unique-KO direct sum for supplied oxidative/associated Pathway",
  stringsAsFactors = FALSE
)

# Additional many-to-one audit membership for the requested core oxidative total.
op_core_membership <- transform(
  op_map_summary[op_map_summary$Feature_order <= 10, ],
  Summary_ID = "OP_CORE_TOTAL",
  Summary_label = "Core oxidative phosphorylation total",
  Summary_order = 22,
  Mapping_basis = "Sum of the first 10 Pathways classified as core oxidative phosphorylation"
)

membership <- rbind(c_map, n_map, p_membership, s_map, op_map_summary, op_core_membership)
membership <- membership[order(factor(membership$Element, levels = c("Carbon","Nitrogen","Phosphorus","Sulfur","Oxidative phosphorylation")), membership$Summary_order, membership$Feature_order), ]
write.table(membership, file.path(work, "CNPS_pathway_summary_membership.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

make_summary_tpm <- function(upper_element, member_element) {
  ids <- unique(member_element$Summary_ID[order(member_element$Summary_order)])
  rows <- lapply(ids, function(sid) {
    mm <- member_element[member_element$Summary_ID == sid, ]
    uu <- upper_element[match(mm$Feature_ID, upper_element$Feature_ID), ]
    data.frame(
      Feature_ID = sid,
      Display_label = mm$Summary_label[1],
      Element = mm$Element[1],
      Section = "Pathway summary",
      Feature_order = mm$Summary_order[1],
      Source_type = "Pathway sum",
      as.list(colSums(uu[, sample_cols, drop = FALSE])),
      check.names = FALSE
    )
  })
  do.call(rbind, rows)
}

summary_tpm <- rbind(
  make_summary_tpm(c_upper, c_map),
  make_summary_tpm(n_upper, n_map),
  make_summary_tpm(p_upper, p_membership),
  make_summary_tpm(s_upper, s_map),
  make_summary_tpm(op_upper, op_map_summary)
)

op_core_total <- data.frame(
  Feature_ID = "OP_CORE_TOTAL",
  Display_label = "Core oxidative phosphorylation total",
  Element = "Oxidative phosphorylation",
  Section = "Pathway summary",
  Feature_order = 22,
  Source_type = "Sum of 10 core oxidative-phosphorylation Pathways",
  as.list(colSums(op_upper[op_upper$Feature_order <= 10, sample_cols, drop = FALSE])),
  check.names = FALSE
)
summary_tpm <- rbind(summary_tpm, op_core_total)
write.table(summary_tpm, file.path(work, "pathway_summary_TPM_by_sample.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(summary_tpm[, 1:6], file.path(work, "pathway_summary_catalog.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

safe_wilcox <- function(x, y) {
  ok <- is.finite(x) & is.finite(y)
  x <- x[ok]; y <- y[ok]
  if (length(x) < 2L || all(x == y)) return(1)
  out <- tryCatch(suppressWarnings(wilcox.test(x, y, paired = TRUE, exact = FALSE)$p.value), error = function(e) NA_real_)
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
  ans <- ifelse(p < 0.001, "***", ifelse(p < 0.01, "**", ifelse(p < 0.05, "*", ifelse(include_ns, "ns", ""))))
  ans[is.na(p)] <- if (include_ns) "ns" else ""
  ans
}

compute_stats <- function(dat, block_name) {
  meta <- dat[, 1:6]
  mat <- as.matrix(dat[, sample_cols, drop = FALSE])
  storage.mode(mat) <- "double"
  alts <- sort(unique(sample_meta$Altitude))
  habitats <- c("R", "E")

  heat <- do.call(rbind, lapply(seq_len(nrow(dat)), function(i) {
    pieces <- list()
    k <- 1L
    for (a in alts) for (h in habitats) {
      idx <- sample_meta$Altitude == a & sample_meta$Habitat == h
      mu <- mean(mat[i, idx])
      pieces[[k]] <- data.frame(Feature_ID = meta$Feature_ID[i], Display_label = meta$Display_label[i],
                                Element = meta$Element[i], Feature_order = meta$Feature_order[i],
                                Block = block_name, Scope = "Altitude", Altitude = a, Habitat = h,
                                Mean_TPM = mu, log2_mean_TPM_plus1 = log2(mu + 1))
      k <- k + 1L
    }
    for (h in habitats) {
      idx <- sample_meta$Habitat == h
      mu <- mean(mat[i, idx])
      pieces[[k]] <- data.frame(Feature_ID = meta$Feature_ID[i], Display_label = meta$Display_label[i],
                                Element = meta$Element[i], Feature_order = meta$Feature_order[i],
                                Block = block_name, Scope = "Total", Altitude = NA_integer_, Habitat = h,
                                Mean_TPM = mu, log2_mean_TPM_plus1 = log2(mu + 1))
      k <- k + 1L
    }
    do.call(rbind, pieces)
  }))

  alt_tests <- do.call(rbind, lapply(seq_len(nrow(dat)), function(i) {
    do.call(rbind, lapply(alts, function(a) {
      r_idx <- which(sample_meta$Altitude == a & sample_meta$Habitat == "R")
      e_idx <- which(sample_meta$Altitude == a & sample_meta$Habitat == "E")
      r_idx <- r_idx[order(sample_meta$Replicate[r_idx])]
      e_idx <- e_idx[order(sample_meta$Replicate[e_idx])]
      data.frame(Feature_ID = meta$Feature_ID[i], Display_label = meta$Display_label[i],
                 Element = meta$Element[i], Block = block_name, Altitude = a,
                 n_pairs = length(r_idx), Mean_R = mean(mat[i, r_idx]), Mean_E = mean(mat[i, e_idx]),
                 Higher = ifelse(mean(mat[i, r_idx]) >= mean(mat[i, e_idx]), "R", "E"),
                 p_raw = safe_wilcox(mat[i, r_idx], mat[i, e_idx]))
    }))
  }))
  alt_tests$p_BH_within_feature_5_altitudes <- ave(alt_tests$p_raw, alt_tests$Feature_ID, FUN = function(x) p.adjust(x, "BH"))
  alt_tests$p_BH_audit_within_altitude_all_features <- ave(alt_tests$p_raw, alt_tests$Altitude, FUN = function(x) p.adjust(x, "BH"))
  alt_tests$Significance <- star(alt_tests$p_BH_within_feature_5_altitudes, FALSE)

  total_tests <- do.call(rbind, lapply(seq_len(nrow(dat)), function(i) {
    r_idx <- which(sample_meta$Habitat == "R")
    e_idx <- which(sample_meta$Habitat == "E")
    r_idx <- r_idx[order(sample_meta$Altitude[r_idx], sample_meta$Replicate[r_idx])]
    e_idx <- e_idx[order(sample_meta$Altitude[e_idx], sample_meta$Replicate[e_idx])]
    data.frame(Feature_ID = meta$Feature_ID[i], Display_label = meta$Display_label[i],
               Element = meta$Element[i], Block = block_name, n_pairs = length(r_idx),
               Mean_R = mean(mat[i, r_idx]), Mean_E = mean(mat[i, e_idx]),
               Higher = ifelse(mean(mat[i, r_idx]) >= mean(mat[i, e_idx]), "R", "E"),
               p_raw = safe_wilcox(mat[i, r_idx], mat[i, e_idx]))
  }))
  total_tests$p_BH_all_features_in_block <- p.adjust(total_tests$p_raw, "BH")
  total_tests$p_BH_within_element <- ave(total_tests$p_raw, total_tests$Element, FUN = function(x) p.adjust(x, "BH"))
  total_tests$Significance <- star(total_tests$p_BH_all_features_in_block, FALSE)

  kw <- do.call(rbind, lapply(seq_len(nrow(dat)), function(i) {
    do.call(rbind, lapply(habitats, function(h) {
      idx <- sample_meta$Habitat == h
      data.frame(Feature_ID = meta$Feature_ID[i], Display_label = meta$Display_label[i],
                 Element = meta$Element[i], Block = block_name, Habitat = h,
                 p_raw = safe_kw(mat[i, idx], sample_meta$Altitude[idx]))
    }))
  }))
  kw$p_BH_all_features_within_habitat <- ave(kw$p_raw, kw$Habitat, FUN = function(x) p.adjust(x, "BH"))
  kw$p_BH_within_element_and_habitat <- ave(kw$p_raw, interaction(kw$Element, kw$Habitat), FUN = function(x) p.adjust(x, "BH"))
  kw$Significance <- star(kw$p_BH_all_features_within_habitat, TRUE)

  list(heat = heat, altitude = alt_tests, total = total_tests, kw = kw)
}

upper_stats <- compute_stats(upper, "Feature heatmap")
summary_stats <- compute_stats(summary_tpm, "Pathway summary")

write.table(upper_stats$heat, file.path(work, "upper_heat_values.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(upper_stats$altitude, file.path(work, "upper_paired_Wilcoxon_by_altitude.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(upper_stats$total, file.path(work, "upper_paired_Wilcoxon_Total.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(upper_stats$kw, file.path(work, "upper_KW_altitude.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(summary_stats$heat, file.path(work, "summary_heat_values.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(summary_stats$altitude, file.path(work, "summary_paired_Wilcoxon_by_altitude.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(summary_stats$total, file.path(work, "summary_paired_Wilcoxon_Total.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(summary_stats$kw, file.path(work, "summary_KW_altitude.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

all_heat <- rbind(upper_stats$heat, summary_stats$heat)
common_fill_max <- ceiling(max(all_heat$log2_mean_TPM_plus1, na.rm = TRUE) * 2) / 2
write.table(data.frame(Transform = "log2(mean raw summed TPM + 1)", Minimum = 0,
                       Maximum = common_fill_max, Palette = "blue-white-red"),
            file.path(work, "shared_fill_scale.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

cat("Upper rows by element:\n")
print(table(upper$Element))
cat("Summary rows by element:\n")
print(table(summary_tpm$Element))
cat(sprintf("Common fill scale: 0 to %.1f\n", common_fill_max))
cat(sprintf("Upper altitude significant (within-feature BH): %d\n", sum(upper_stats$altitude$Significance != "")))
cat(sprintf("Upper Total significant (block-wide BH): %d\n", sum(upper_stats$total$Significance != "")))
cat(sprintf("Summary altitude significant (within-summary BH): %d\n", sum(summary_stats$altitude$Significance != "")))
cat(sprintf("Summary Total significant (block-wide BH): %d\n", sum(summary_stats$total$Significance != "")))
