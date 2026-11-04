options(stringsAsFactors = FALSE, scipen = 999)

root <- "/Volumes/Haleo_DK01/WD_Takakia_F/metagenome-data/Figure7_CNPS"
out <- file.path(Sys.getenv("TAKAKIA_PROJECT_DIR", unset = "~/Documents/Takakia"), "Figure4_final_tables_2026-07-26/work")
dir.create(out, recursive = TRUE, showWarnings = FALSE)

read_tsv <- function(path) {
  read.delim(path, check.names = FALSE, quote = "", comment.char = "",
             na.strings = c("", "NA"))
}

sample_cols <- function(x) grep("^[PS][0-9]{2}_[1-6]$", names(x), value = TRUE)

ko_tpm <- read_tsv("/Volumes/Haleo_DK01/WD_Takakia_F/metagenome-data/tpm_renormalized.txt")
names(ko_tpm)[1] <- "KO"
ko_samples <- sample_cols(ko_tpm)
ko_tpm[ko_samples] <- lapply(ko_tpm[ko_samples], as.numeric)

final_cat <- read_tsv(file.path(root, "Figure7_final_CNPS_KO_catalog_all_candidates.tsv"))
official_defs <- read_tsv(file.path(root, "Figure7_KEGG_official_KO_definitions.tsv"))

ref_files <- file.path(root, paste0("Figure7_CNPS_literature_recurated_", c("C", "N", "P", "S"), ".tsv"))
refs <- do.call(rbind, lapply(ref_files, read_tsv))
names(refs)[names(refs) == "K-number"] <- "KO"

collapse_unique <- function(x) {
  x <- unique(trimws(as.character(x[!is.na(x) & nzchar(trimws(as.character(x)))])))
  if (!length(x)) "" else paste(x, collapse = " | ")
}

ref_by_ko <- do.call(rbind, lapply(split(refs, refs$KO), function(z) {
  data.frame(
    KO = z$KO[1],
    Literature_reference = collapse_unique(z$reference),
    Literature_pathway = collapse_unique(z$`Process / pathway`),
    Literature_function_description = collapse_unique(z$`function description`),
    Literature_note = collapse_unique(z$`Curation note`),
    stringsAsFactors = FALSE
  )
}))
row.names(ref_by_ko) <- NULL

def_by_ko <- do.call(rbind, lapply(split(final_cat, final_cat$KO), function(z) {
  data.frame(
    KO = z$KO[1],
    KEGG_definition_master = collapse_unique(z$KEGG_definition),
    KEGG_module_or_pathway_master = collapse_unique(z$KEGG_module_or_pathway),
    KEGG_URL_master = collapse_unique(z$KEGG_URL),
    stringsAsFactors = FALSE
  )
}))
row.names(def_by_ko) <- NULL
def_by_ko <- merge(
  official_defs[, c("KO", "KEGG_definition", "KEGG_URL")],
  def_by_ko, by = "KO", all = TRUE, sort = FALSE
)
def_by_ko$KEGG_definition_master <- ifelse(
  !is.na(def_by_ko$KEGG_definition) &
    nzchar(def_by_ko$KEGG_definition),
  def_by_ko$KEGG_definition,
  def_by_ko$KEGG_definition_master
)
def_by_ko$KEGG_URL_master <- ifelse(
  !is.na(def_by_ko$KEGG_URL) & nzchar(def_by_ko$KEGG_URL),
  def_by_ko$KEGG_URL,
  def_by_ko$KEGG_URL_master
)
def_by_ko <- def_by_ko[, c(
  "KO", "KEGG_definition_master", "KEGG_module_or_pathway_master",
  "KEGG_URL_master"
)]

features <- data.frame(
  Feature_ID = character(), Panel = character(), Element = character(),
  Feature_level = character(), Section = character(), Display_label = character(),
  Figure_order = integer(), Source_system = character(), stringsAsFactors = FALSE
)
membership <- data.frame(
  Feature_ID = character(), Member_ID = character(), KO = character(),
  Gene_label = character(), Function_description = character(),
  Module_map = character(), Source_URL = character(),
  Abundance_source = character(), stringsAsFactors = FALSE
)

add_feature <- function(id, panel, element, level, section, label, order, system) {
  features <<- rbind(features, data.frame(
    Feature_ID = id, Panel = panel, Element = element,
    Feature_level = level, Section = section, Display_label = label,
    Figure_order = order, Source_system = system, stringsAsFactors = FALSE
  ))
}

add_member <- function(fid, member_id, ko, gene, fun = "", module = "",
                       url = "", source = "") {
  membership <<- rbind(membership, data.frame(
    Feature_ID = fid, Member_ID = member_id, KO = ko, Gene_label = gene,
    Function_description = fun, Module_map = module, Source_URL = url,
    Abundance_source = source, stringsAsFactors = FALSE
  ))
}

## Carbon: select only the tiles actually retained in final Figure 4.
c_v2_dir <- file.path(root, "C_maintext_dual_TCA_rTCA_2026-07-22/work/V2_TCA")
c_v2_cat <- read_tsv(file.path(c_v2_dir, "C_label_catalog.tsv"))
c_v2_mem <- read_tsv(file.path(c_v2_dir, "C_label_KO_membership.tsv"))
c_rt_dir <- file.path(root, "C_maintext_rTCA_3HP_3labels_2026-07-22/work/V1_rTCA_3HP3")
c_rt_cat <- read_tsv(file.path(c_rt_dir, "C_label_catalog.tsv"))
c_rt_mem <- read_tsv(file.path(c_rt_dir, "C_label_KO_membership.tsv"))

c_base_labels <- c(
  "Chitin depolymerisation", "Cellulose depolymerisation",
  "Arabinoxylan debranching", "Mannan depolymerisation",
  "Xylan depolymerisation", "Pectate depolymerisation",
  "Pectin de-esterification", "Pectin depolymerisation",
  "Alpha-amylase", "Pullulanase", "Glucoamylase", "Isopullulanase",
  "Protocatechuate 3,4-ring cleavage", "Protocatechuate 4,5-ring cleavage",
  "Vanillate O-demethylation", "Vanillin dehydrogenase",
  "Phenolic compound oxidation", "Glyoxylate shunt",
  "Citrate synthase", "Aconitate hydratase", "Isocitrate dehydrogenase",
  "2-Oxoglutarate dehydrogenase complex", "Succinyl-CoA synthetase",
  "Succinate dehydrogenase", "Fumarate hydratase", "Malate dehydrogenase",
  "sMMO — mmoX", "sMMO — mmoC", "Aerobic CO oxidation",
  "CFeSP methyltransferase", "CODH/ACS complex gamma subunit",
  "W–L methyl branch", "Pyruvate-to-acetate fermentation",
  "Ethanol fermentation", "Butyrate fermentation",
  "L-Lactate fermentation", "D-Lactate fermentation"
)

c_sections <- c(
  "Chitin", "Cellulose and hemicellulose", "Cellulose and hemicellulose",
  "Cellulose and hemicellulose", "Cellulose and hemicellulose",
  "Pectin", "Pectin", "Pectin", "Starch", "Starch", "Starch", "Starch",
  "Aromatics", "Aromatics", "Phenolic compounds", "Phenolic compounds",
  "Phenolic compounds", "Glyoxylate shunt",
  rep("Tricarboxylic acid (TCA) cycle", 8),
  "Methane cycling", "Methane cycling", "Aerobic CO oxidation",
  rep("Wood–Ljungdahl pathway", 3), rep("Fermentation", 5)
)

c_order <- 1L
for (i in seq_along(c_base_labels)) {
  lab <- c_base_labels[i]
  src <- c_v2_cat[c_v2_cat$Display_label == lab, , drop = FALSE]
  if (nrow(src) != 1) stop("Carbon label not unique/missing: ", lab)
  fid <- paste0("C_TILE_", sprintf("%02d", c_order))
  display_lab <- gsub("—", "–", lab, fixed = TRUE)
  add_feature(fid, "a", "Carbon", "Figure tile", c_sections[i], display_lab,
              c_order, "KO TPM")
  mm <- c_v2_mem[c_v2_mem$Feature_ID == src$Feature_ID, , drop = FALSE]
  for (j in seq_len(nrow(mm))) {
    add_member(fid, mm$KO[j], mm$KO[j], "", "", "",
               paste0("https://www.kegg.jp/entry/", mm$KO[j]),
               "tpm_renormalized.txt")
  }
  c_order <- c_order + 1L
}

## Split the CBB tile exactly as displayed in the final figure.
for (spec in list(
  list("Phosphoribulokinase", "K00855"),
  list("RuBisCO", c("K01601", "K01602"))
)) {
  fid <- paste0("C_TILE_", sprintf("%02d", c_order))
  add_feature(fid, "a", "Carbon", "Figure tile",
              "Calvin–Benson–Bassham cycle", spec[[1]], c_order, "KO TPM")
  for (ko in spec[[2]]) add_member(fid, ko, ko, "", "", "",
                                   paste0("https://www.kegg.jp/entry/", ko),
                                   "tpm_renormalized.txt")
  c_order <- c_order + 1L
}

rt_labels <- c(
  "Pyruvate:ferredoxin oxidoreductase (PFOR)",
  "2-Oxoglutarate:ferredoxin oxidoreductase",
  "Citryl-CoA route (incomplete)",
  "3-Hydroxypropionate formation and propionyl-CoA carboxylation",
  "First-cycle closure and glyoxylate formation",
  "Second-cycle C5 intermediate conversion"
)
rt_display <- c(
  "Pyruvate:ferredoxin oxidoreductase",
  "2-Oxoglutarate:ferredoxin oxidoreductase",
  "Citryl-CoA route", "3-HP formation", "First-cycle", "Second-cycle"
)
rt_sections <- c(rep("Reductive TCA cycle", 3), rep("3-Hydroxypropionate (3-HP)", 3))
for (i in seq_along(rt_labels)) {
  src <- c_rt_cat[c_rt_cat$Display_label == rt_labels[i], , drop = FALSE]
  if (nrow(src) != 1) stop("rTCA/3HP label not unique/missing: ", rt_labels[i])
  fid <- paste0("C_TILE_", sprintf("%02d", c_order))
  add_feature(fid, "a", "Carbon", "Figure tile", rt_sections[i],
              rt_display[i], c_order, "KO TPM")
  mm <- c_rt_mem[c_rt_mem$Feature_ID == src$Feature_ID, , drop = FALSE]
  for (j in seq_len(nrow(mm))) {
    add_member(fid, mm$KO[j], mm$KO[j], "", "", "",
               paste0("https://www.kegg.jp/entry/", mm$KO[j]),
               "tpm_renormalized.txt")
  }
  c_order <- c_order + 1L
}

## Nitrogen tiles, including EggNOG-resolved narG/narZ and the displayed ureDEFHJ group.
n_dir <- file.path(root, "N_ure_sums_2026-07-24/work/N")
n_cat <- read_tsv(file.path(n_dir, "N_gene_catalog.tsv"))
n_tpm <- read_tsv(file.path(n_dir, "N_gene_TPM_by_sample.tsv"))
n_samples <- sample_cols(n_tpm)
n_single_ids <- c(
  "N_NEW_K02575", "K15576", "K15577", "K15578", "K15579",
  "K00372_nasA", "K00360", "K00367",
  "K00370_narG", "K00370_narZ", "K00371", "K00374",
  "K02567", "K02568", "K02569", "K10535", "K00368",
  "K04561", "K02305", "K00376", "K02584", "K04487", "K02586",
  "K02588", "K02591", "K03385", "K15876", "K00362", "K00363",
  "K00366", "K01427", "K01428", "K01429", "K01430", "K14048",
  "N_NEW_K00260", "K00261", "K00262", "K15371", "K00284",
  "K00265", "K00266", "K01915",
  "K11959", "K11960", "K11961", "K11962", "K11963"
)
n_overrides <- c(
  N_NEW_K02575 = "NRT2", N_NEW_K00260 = "gudB",
  K00262 = "gdh_K00262", K15371 = "gdh_K15371"
)
n_order <- 1L
for (sid in n_single_ids) {
  z <- n_cat[n_cat$Feature_ID == sid, , drop = FALSE]
  if (nrow(z) != 1) stop("Nitrogen source ID not unique/missing: ", sid)
  fid <- paste0("N_TILE_", sprintf("%02d", n_order))
  lab <- if (sid %in% names(n_overrides)) unname(n_overrides[sid]) else z$Display_label
  add_feature(fid, "b", "Nitrogen", "Figure tile", z$Section, lab,
              n_order, "N resolved-feature TPM")
  add_member(fid, sid, z$KO, lab, z$Function_description, z$Module_map,
             z$Source_url, z$Source_type)
  n_order <- n_order + 1L
}
n_acc_ids <- c("K03187", "K03188", "K03189", "K03190", "K03192")
fid <- paste0("N_TILE_", sprintf("%02d", n_order))
add_feature(fid, "b", "Nitrogen", "Figure tile",
            "Urea hydrolysis and urease maturation", "ureDEFHJ",
            n_order, "N resolved-feature TPM")
for (sid in n_acc_ids) {
  z <- n_cat[n_cat$Feature_ID == sid, , drop = FALSE]
  add_member(fid, sid, z$KO, z$Display_label, z$Function_description,
             z$Module_map, z$Source_url, z$Source_type)
}

## Phosphorus tiles.
p_dir <- file.path(root, "N_S_P_altitude_sums_2026-07-23/work/P")
p_cat <- read_tsv(file.path(p_dir, "P_gene_catalog.tsv"))
p_tpm <- read_tsv(file.path(p_dir, "P_gene_TPM_by_sample.tsv"))
p_samples <- sample_cols(p_tpm)
p_single_labels <- c(
  "phnC", "phnD", "phnE", "phnF", "phnA", "phnW", "phnX",
  "gcd", "pqqB", "pqqC", "pstA", "pstB", "pstC", "pstS", "pit",
  "phoB", "phoP", "phoR", "phoU",
  "phoA", "phoD", "phoN", "aphA", "PHO",
  "ugpQ", "ugpA", "ugpB", "ugpC", "ugpE",
  "ppk1", "ppa", "ppx"
)
p_order <- 1L
for (lab in p_single_labels) {
  z <- p_cat[p_cat$Display_label == lab, , drop = FALSE]
  if (nrow(z) != 1) stop("Phosphorus label not unique/missing: ", lab)
  fid <- paste0("P_TILE_", sprintf("%02d", p_order))
  add_feature(fid, "c", "Phosphorus", "Figure tile", z$Section, lab,
              p_order, "P KO TPM")
  add_member(fid, z$Feature_ID, z$KO, lab, z$Function_description,
             z$Module_map, z$Source_url, z$Source_type)
  p_order <- p_order + 1L
}
p_phn_group <- c("phnG", "phnH", "phnI", "phnJ", "phnK",
                 "phnL", "phnM", "phnN", "phnO", "phnP")
fid <- paste0("P_TILE_", sprintf("%02d", p_order))
add_feature(fid, "c", "Phosphorus", "Figure tile",
            "Phosphonate uptake and C-P cleavage", "phnGHIJKLMNOP",
            p_order, "P KO TPM")
for (lab in p_phn_group) {
  z <- p_cat[p_cat$Display_label == lab, , drop = FALSE]
  add_member(fid, z$Feature_ID, z$KO, lab, z$Function_description,
             z$Module_map, z$Source_url, z$Source_type)
}

## Sulfur tiles.
s_dir <- file.path(root, "N_S_P_altitude_sums_2026-07-23/work/S")
s_cat <- read_tsv(file.path(s_dir, "S_gene_catalog.tsv"))
s_tpm <- read_tsv(file.path(s_dir, "S_gene_TPM_by_sample.tsv"))
s_samples <- sample_cols(s_tpm)
s_single_ids <- c(
  "S_NEW_K13811", "K00958", "K00957", "K00956", "K00860", "K00390",
  "K05301", "K21307", "K21308", "K21309",
  "S_TRANS_K02045", "S_TRANS_K02046", "S_TRANS_K02047",
  "S_TRANS_K02048", "S_TRANS_K06203", "S_TRANS_K03321",
  "K17222", "K17223", "K17224", "K17225", "K17226", "K17227",
  "K08357", "K16937", "K19713", "K08354",
  "K16950", "K16951", "K00385", "K00392", "K00381", "K00380",
  "K17230", "K17229", "K17218",
  "S_NEW_K00640", "S_NEW_K01738",
  "K07306", "K00184", "K00185", "K16966",
  "S_TRANS_K15553", "S_TRANS_K15554", "S_TRANS_K15555",
  "S_TRANS_K15551", "S_TRANS_K15552", "S_TRANS_K10831"
)
s_overrides <- c(S_TRANS_K02046 = "cysU", S_TRANS_K02048 = "cysP")
s_order <- 1L
for (sid in s_single_ids) {
  z <- s_cat[s_cat$Feature_ID == sid, , drop = FALSE]
  if (nrow(z) != 1) stop("Sulfur source ID not unique/missing: ", sid)
  fid <- paste0("S_TILE_", sprintf("%02d", s_order))
  lab <- if (sid %in% names(s_overrides)) unname(s_overrides[sid]) else z$Display_label
  add_feature(fid, "d", "Sulfur", "Figure tile", z$Section, lab,
              s_order, "S KO TPM")
  add_member(fid, sid, z$KO, lab, z$Function_description, z$Module_map,
             z$Source_url, z$Source_type)
  s_order <- s_order + 1L
}
fid <- paste0("S_TILE_", sprintf("%02d", s_order))
add_feature(fid, "d", "Sulfur", "Figure tile", "Dissimilatory sulfate reduction",
            "aprAB", s_order, "S KO TPM")
for (sid in c("K00394", "K00395")) {
  z <- s_cat[s_cat$Feature_ID == sid, , drop = FALSE]
  add_member(fid, sid, z$KO, z$Display_label, z$Function_description,
             z$Module_map, z$Source_url, z$Source_type)
}

## Oxidative phosphorylation: exactly the ten pathway rows shown beneath panel c.
ox_dir <- file.path(root, "CNPS_altitude_sum_revision_2026-07-21/V2_3_categories/work")
ox_map <- read_tsv(file.path(ox_dir, "oxidative_pathway_KO_mapping.tsv"))
ox_tpm <- read_tsv(file.path(ox_dir, "oxidative_pathway_summed_TPM_by_sample.tsv"))
ox_names <- c(
  "F-type ATPase", "V/A-type ATPase", "NADH-quinone oxidoreductase",
  "NAD(P)H-quinone oxidoreductase", "Succinate dehydrogenase (ubiquinone)",
  "Cytochrome c oxidase, cbb3-type", "Cytochrome bd ubiquinol oxidase",
  "Cytochrome o ubiquinol oxidase",
  "Cytochrome c oxidase, prokaryotes, aa3-type", "Cytochrome bc1 complex"
)
for (i in seq_along(ox_names)) {
  fid <- paste0("OP_TILE_", sprintf("%02d", i))
  add_feature(fid, "c", "Oxidative phosphorylation", "Pathway tile",
              "Core oxidative phosphorylation", ox_names[i], i,
              "Direct KO sum")
  z <- ox_map[ox_map$Pathway == ox_names[i] & !duplicated(ox_map$KO), , drop = FALSE]
  for (j in seq_len(nrow(z))) {
    add_member(fid, z$KO[j], z$KO[j], "", z$KEGG_definition[j],
               "map00190", paste0("https://www.kegg.jp/entry/", z$KO[j]),
               "tpm_renormalized.txt; pathway definition from P_pathway/DiTing review")
  }
}

## Attach authoritative KO definitions and literature references.
membership <- merge(membership, def_by_ko, by = "KO", all.x = TRUE, sort = FALSE)
membership <- merge(membership, ref_by_ko, by = "KO", all.x = TRUE, sort = FALSE)
membership$Function_description[
  is.na(membership$Function_description) | !nzchar(membership$Function_description)
] <- membership$KEGG_definition_master[
  is.na(membership$Function_description) | !nzchar(membership$Function_description)
]
membership$Function_description[
  is.na(membership$Function_description) | !nzchar(membership$Function_description)
] <- membership$Literature_function_description[
  is.na(membership$Function_description) | !nzchar(membership$Function_description)
]
membership$Module_map[
  is.na(membership$Module_map) | !nzchar(membership$Module_map)
] <- membership$KEGG_module_or_pathway_master[
  is.na(membership$Module_map) | !nzchar(membership$Module_map)
]
membership$Source_URL[
  is.na(membership$Source_URL) | !nzchar(membership$Source_URL)
] <- membership$KEGG_URL_master[
  is.na(membership$Source_URL) | !nzchar(membership$Source_URL)
]

## Build sample-level TPM for every displayed tile.
samples <- ko_samples
tile_tpm <- features
for (s in samples) tile_tpm[[s]] <- 0

sum_ko_rows <- function(kos) {
  kos <- unique(kos)
  z <- ko_tpm[ko_tpm$KO %in% kos, samples, drop = FALSE]
  if (!nrow(z)) return(setNames(rep(0, length(samples)), samples))
  colSums(z, na.rm = TRUE)
}

sum_source_rows <- function(tab, ids) {
  sc <- sample_cols(tab)
  z <- tab[tab$Feature_ID %in% unique(ids), sc, drop = FALSE]
  if (!nrow(z)) return(setNames(rep(0, length(sc)), sc))
  colSums(z, na.rm = TRUE)
}

for (i in seq_len(nrow(features))) {
  fid <- features$Feature_ID[i]
  mm <- membership[membership$Feature_ID == fid, , drop = FALSE]
  el <- features$Element[i]
  vals <- if (el %in% c("Carbon", "Oxidative phosphorylation")) {
    sum_ko_rows(mm$KO)
  } else if (el == "Nitrogen") {
    sum_source_rows(n_tpm, mm$Member_ID)
  } else if (el == "Phosphorus") {
    sum_source_rows(p_tpm, mm$Member_ID)
  } else {
    sum_source_rows(s_tpm, mm$Member_ID)
  }
  tile_tpm[i, names(vals)] <- as.list(vals)
}

## Aggregate rows correspond to the circle/pathway sums used for interpretation.
aggregates <- data.frame(
  Feature_ID = character(), Panel = character(), Element = character(),
  Feature_level = character(), Section = character(), Display_label = character(),
  Figure_order = integer(), Source_system = character(), stringsAsFactors = FALSE
)
agg_members <- data.frame(Aggregate_ID = character(), Tile_ID = character(),
                          stringsAsFactors = FALSE)

add_aggregate <- function(id, panel, element, section, label, tile_ids, order) {
  aggregates <<- rbind(aggregates, data.frame(
    Feature_ID = id, Panel = panel, Element = element,
    Feature_level = "Aggregate (circle/pathway sum)", Section = section,
    Display_label = label, Figure_order = order,
    Source_system = "Within-sample sum of displayed tiles",
    stringsAsFactors = FALSE
  ))
  agg_members <<- rbind(agg_members, data.frame(
    Aggregate_ID = id, Tile_ID = tile_ids, stringsAsFactors = FALSE
  ))
}

tile_ids_by <- function(element, labels) {
  features$Feature_ID[features$Element == element & features$Display_label %in% labels]
}

## Carbon circles.
c_agg_specs <- list(
  Chitin = c("Chitin depolymerisation"),
  `Cellulose and hemicellulose` = c("Cellulose depolymerisation", "Arabinoxylan debranching",
                                    "Mannan depolymerisation", "Xylan depolymerisation"),
  Pectin = c("Pectate depolymerisation", "Pectin de-esterification", "Pectin depolymerisation"),
  Starch = c("Alpha-amylase", "Pullulanase", "Glucoamylase", "Isopullulanase"),
  Aromatics = c("Protocatechuate 3,4-ring cleavage", "Protocatechuate 4,5-ring cleavage"),
  `Phenolic compounds` = c("Vanillate O-demethylation", "Vanillin dehydrogenase",
                           "Phenolic compound oxidation"),
  `Glyoxylate shunt` = c("Glyoxylate shunt"),
  `Tricarboxylic acid (TCA) cycle` = c("Citrate synthase", "Aconitate hydratase",
    "Isocitrate dehydrogenase", "2-Oxoglutarate dehydrogenase complex",
    "Succinyl-CoA synthetase", "Succinate dehydrogenase", "Fumarate hydratase",
    "Malate dehydrogenase"),
  `Reductive TCA cycle` = c("Pyruvate:ferredoxin oxidoreductase",
    "2-Oxoglutarate:ferredoxin oxidoreductase", "Citryl-CoA route"),
  `Methane cycling` = c("sMMO – mmoX", "sMMO – mmoC"),
  `Aerobic CO oxidation` = c("Aerobic CO oxidation"),
  `Calvin–Benson–Bassham cycle` = c("Phosphoribulokinase", "RuBisCO"),
  `Wood–Ljungdahl pathway` = c("CFeSP methyltransferase",
    "CODH/ACS complex gamma subunit", "W–L methyl branch"),
  `3-Hydroxypropionate (3-HP)` = c("3-HP formation", "First-cycle", "Second-cycle"),
  Fermentation = c("Pyruvate-to-acetate fermentation", "Ethanol fermentation",
                   "Butyrate fermentation", "L-Lactate fermentation",
                   "D-Lactate fermentation")
)
ao <- 1L
for (nm in names(c_agg_specs)) {
  add_aggregate(paste0("C_SUM_", sprintf("%02d", ao)), "a", "Carbon",
                nm, nm, tile_ids_by("Carbon", c_agg_specs[[nm]]), ao)
  ao <- ao + 1L
}

## N, P, S and oxidative-phosphorylation sums.
agg_specs <- list(
  list("N_SUM_01", "b", "Nitrogen", "Nitrate/nitrite transport",
       "NRT2 + nrtABCD", c("NRT2", "nrtA", "nrtB", "nrtC", "nrtD")),
  list("N_SUM_02", "b", "Nitrogen", "Assimilatory nitrate reduction",
       "nasAB + narB", c("nasA", "nasB", "narB")),
  list("N_SUM_03", "b", "Nitrogen", "Membrane/periplasmic nitrate reduction",
       "narGHI + napAB", c("narG", "narH", "narI", "napA", "napB")),
  list("N_SUM_04", "b", "Nitrogen", "Nitrogen fixation",
       "Nitrogen fixation", c("nifA", "nifS", "nifD", "nifH", "nifK")),
  list("N_SUM_05", "b", "Nitrogen", "DNRA",
       "DNRA nitrite reduction", c("nrfA", "nrfH", "nirB", "nirD")),
  list("N_SUM_06", "b", "Nitrogen", "Urea transport",
       "Urea transport", c("urtA", "urtB", "urtC", "urtD", "urtE")),
  list("N_SUM_07", "b", "Nitrogen", "Urea hydrolysis",
       "URE + ureAB + ureCBA", c("URE", "ureAB", "ureC", "ureB", "ureA")),
  list("N_SUM_08", "b", "Nitrogen", "Urease maturation",
       "Urease accessory system", c("ureDEFHJ")),
  list("N_SUM_09", "b", "Nitrogen", "Urea hydrolysis",
       "Urease total", c("URE", "ureAB", "ureC", "ureB", "ureA", "ureDEFHJ")),
  list("N_SUM_10", "b", "Nitrogen", "Ammonium assimilation",
       "GDH + GS-GOGAT", c("gudB", "gdhA", "gdh_K00262", "gdh_K15371",
                           "GLU", "gltB", "gltD", "glnA")),
  list("P_SUM_01", "c", "Phosphorus", "Organic phosphorus mineralisation",
       "phoADN + aphA + PHO", c("phoA", "phoD", "phoN", "aphA", "PHO")),
  list("P_SUM_02", "c", "Phosphorus", "Glycerol-3-phosphate uptake/mineralisation",
       "ugpABCEQ", c("ugpA", "ugpB", "ugpC", "ugpE", "ugpQ")),
  list("P_SUM_03", "c", "Phosphorus", "Organic phosphorus accessibility",
       "phoADN + aphA + PHO + ugpABCEQ",
       c("phoA", "phoD", "phoN", "aphA", "PHO", "ugpA", "ugpB", "ugpC", "ugpE", "ugpQ")),
  list("P_SUM_04", "c", "Phosphorus", "Inorganic phosphate uptake",
       "pstABCS + pit", c("pstA", "pstB", "pstC", "pstS", "pit")),
  list("P_SUM_05", "c", "Phosphorus", "Phosphonate uptake",
       "phnCDEF", c("phnC", "phnD", "phnE", "phnF")),
  list("P_SUM_06", "c", "Phosphorus", "C-P lyase",
       "phnGHIJKLMNOP", c("phnGHIJKLMNOP")),
  list("P_SUM_07", "c", "Phosphorus", "Phosphonoacetate utilisation",
       "phnAWX", c("phnA", "phnW", "phnX")),
  list("S_SUM_01", "d", "Sulfur", "Sulfate activation",
       "PAPSS + sat + cysND", c("PAPSS", "sat", "cysN", "cysD")),
  list("S_SUM_02", "d", "Sulfur", "Anaerobic sulfite reduction",
       "asrABC", c("asrA", "asrB", "asrC")),
  list("S_SUM_03", "d", "Sulfur", "Assimilatory sulfite reduction",
       "sir + cysIJ", c("sir", "cysI", "cysJ")),
  list("S_SUM_04", "d", "Sulfur", "Cysteine biosynthesis",
       "cysEK", c("cysE", "cysK")),
  list("S_SUM_05", "d", "Sulfur", "Sulfide oxidation",
       "fccAB", c("fccA", "fccB")),
  list("S_SUM_06", "d", "Sulfur", "Tetrathionate intermediate oxidation",
       "doxD + tsdA", c("doxD", "tsdA")),
  list("S_SUM_07", "d", "Sulfur", "Sulfate/thiosulfate transport",
       "cysAUWP + cysZ + sulP", c("cysA", "cysU", "cysW", "cysP", "cysZ", "sulP")),
  list("S_SUM_08", "d", "Sulfur", "Alkanesulfonate transport",
       "ssuABC", c("ssuA", "ssuB", "ssuC")),
  list("S_SUM_09", "d", "Sulfur", "Taurine transport",
       "tauABC", c("tauA", "tauB", "tauC")),
  list("S_SUM_10", "d", "Sulfur", "SOX thiosulfate oxidation",
       "soxAXBCYZ", c("soxA", "soxX", "soxB", "soxC", "soxY", "soxZ")),
  list("S_SUM_11", "d", "Sulfur", "Sulfite oxidation",
       "sorA + soeABC", c("sorA", "soeA", "soeB", "soeC")),
  list("S_SUM_12", "d", "Sulfur", "DMSO reduction",
       "dmsABC", c("dmsA", "dmsB", "dmsC")),
  list("OP_SUM_01", "c", "Oxidative phosphorylation", "Core oxidative phosphorylation",
       "Core oxidative phosphorylation total", ox_names)
)
for (i in seq_along(agg_specs)) {
  z <- agg_specs[[i]]
  add_aggregate(z[[1]], z[[2]], z[[3]], z[[4]], z[[5]],
                tile_ids_by(z[[3]], z[[6]]), i)
}

## Fail loudly if any aggregate label failed to resolve.
agg_counts <- table(agg_members$Aggregate_ID)
empty_aggs <- aggregates$Feature_ID[!aggregates$Feature_ID %in% names(agg_counts)]
if (length(empty_aggs)) stop("Aggregates with no members: ", paste(empty_aggs, collapse = ", "))

agg_tpm <- aggregates
for (s in samples) agg_tpm[[s]] <- 0
for (i in seq_len(nrow(aggregates))) {
  ids <- agg_members$Tile_ID[agg_members$Aggregate_ID == aggregates$Feature_ID[i]]
  z <- tile_tpm[tile_tpm$Feature_ID %in% ids, samples, drop = FALSE]
  agg_tpm[i, samples] <- as.list(colSums(z, na.rm = TRUE))
}

all_feature_tpm <- rbind(tile_tpm, agg_tpm)

## Statistical tests on raw TPM. Altitude BH correction is applied within each
## displayed feature across the five altitudes, as requested by the user.
star <- function(q) {
  ifelse(is.na(q), "ns",
    ifelse(q < 0.001, "***", ifelse(q < 0.01, "**", ifelse(q < 0.05, "*", "ns"))))
}
safe_wilcox <- function(x, y) {
  if (all(is.na(x)) || all(is.na(y))) return(NA_real_)
  d <- x - y
  if (all(abs(d[is.finite(d)]) < .Machine$double.eps^0.5)) return(1)
  suppressWarnings(tryCatch(wilcox.test(x, y, paired = TRUE, exact = FALSE)$p.value,
                            error = function(e) NA_real_))
}
safe_kw <- function(v, g) {
  if (length(unique(g[is.finite(v)])) < 2) return(NA_real_)
  suppressWarnings(tryCatch(kruskal.test(v ~ factor(g))$p.value,
                            error = function(e) NA_real_))
}

alts <- c("38", "39", "40", "41", "42")
alt_long <- data.frame()
total_stats <- data.frame()
kw_stats <- data.frame()
for (i in seq_len(nrow(all_feature_tpm))) {
  z <- all_feature_tpm[i, ]
  tmp <- data.frame()
  for (a in alts) {
    ecols <- paste0("P", a, "_", 1:6)
    rcols <- paste0("S", a, "_", 1:6)
    e <- as.numeric(z[ecols])
    r <- as.numeric(z[rcols])
    tmp <- rbind(tmp, data.frame(
      Feature_ID = z$Feature_ID, Panel = z$Panel, Element = z$Element,
      Feature_level = z$Feature_level, Section = z$Section,
      Display_label = z$Display_label, Altitude_m = as.integer(a) * 100,
      Mean_R_TPM = mean(r, na.rm = TRUE), Mean_E_TPM = mean(e, na.rm = TRUE),
      Higher = ifelse(mean(r, na.rm = TRUE) > mean(e, na.rm = TRUE), "R",
                      ifelse(mean(e, na.rm = TRUE) > mean(r, na.rm = TRUE), "E", "Equal")),
      Wilcoxon_p = safe_wilcox(r, e), stringsAsFactors = FALSE
    ))
  }
  tmp$BH_q_within_feature <- p.adjust(tmp$Wilcoxon_p, method = "BH")
  tmp$Significance <- star(tmp$BH_q_within_feature)
  alt_long <- rbind(alt_long, tmp)

  e_all <- as.numeric(z[paste0("P", rep(alts, each = 6), "_", rep(1:6, times = 5))])
  r_all <- as.numeric(z[paste0("S", rep(alts, each = 6), "_", rep(1:6, times = 5))])
  total_stats <- rbind(total_stats, data.frame(
    Feature_ID = z$Feature_ID, Mean_R_TPM = mean(r_all, na.rm = TRUE),
    Mean_E_TPM = mean(e_all, na.rm = TRUE),
    Higher = ifelse(mean(r_all, na.rm = TRUE) > mean(e_all, na.rm = TRUE), "R",
                    ifelse(mean(e_all, na.rm = TRUE) > mean(r_all, na.rm = TRUE), "E", "Equal")),
    Wilcoxon_p = safe_wilcox(r_all, e_all), stringsAsFactors = FALSE
  ))
  kw_stats <- rbind(kw_stats, data.frame(
    Feature_ID = z$Feature_ID,
    KW_R_p = safe_kw(r_all, rep(as.integer(alts) * 100, each = 6)),
    KW_E_p = safe_kw(e_all, rep(as.integer(alts) * 100, each = 6)),
    stringsAsFactors = FALSE
  ))
}
total_stats$BH_q_across_features <- p.adjust(total_stats$Wilcoxon_p, method = "BH")
total_stats$Significance <- star(total_stats$BH_q_across_features)
kw_stats$KW_R_BH_q <- p.adjust(kw_stats$KW_R_p, method = "BH")
kw_stats$KW_E_BH_q <- p.adjust(kw_stats$KW_E_p, method = "BH")
kw_stats$KW_R_Significance <- star(kw_stats$KW_R_BH_q)
kw_stats$KW_E_Significance <- star(kw_stats$KW_E_BH_q)

## Reader-friendly compact matrix: direction plus significance.
summary_tab <- all_feature_tpm[, c("Feature_ID", "Panel", "Element",
                                   "Feature_level", "Section", "Display_label")]
for (a in as.integer(alts) * 100) {
  z <- alt_long[alt_long$Altitude_m == a, ]
  z <- z[match(summary_tab$Feature_ID, z$Feature_ID), ]
  summary_tab[[paste0(a, " m")]] <- ifelse(
    z$Significance == "ns", "ns",
    paste0(z$Higher, z$Significance)
  )
}
tot <- total_stats[match(summary_tab$Feature_ID, total_stats$Feature_ID), ]
kw <- kw_stats[match(summary_tab$Feature_ID, kw_stats$Feature_ID), ]
summary_tab$Total <- ifelse(tot$Significance == "ns", "ns",
                            paste0(tot$Higher, tot$Significance))
summary_tab$Mean_R_TPM_total <- tot$Mean_R_TPM
summary_tab$Mean_E_TPM_total <- tot$Mean_E_TPM
summary_tab$KW_R <- kw$KW_R_Significance
summary_tab$KW_E <- kw$KW_E_Significance

## Final KO list: one row per KO membership under a Figure 4 tile.
final_list <- merge(features, membership, by = "Feature_ID", all.x = TRUE, sort = FALSE)
final_list <- final_list[order(match(final_list$Feature_ID, features$Feature_ID),
                               final_list$KO), ]
missing_fun <- is.na(final_list$Function_description) |
  !nzchar(final_list$Function_description)
final_list$Function_description[missing_fun] <- paste0(
  "KO member of Figure 4 feature: ", final_list$Display_label[missing_fun]
)
missing_gene <- is.na(final_list$Gene_label) | !nzchar(final_list$Gene_label)
derived_gene <- sub(";.*$", "", final_list$Function_description)
derived_gene[grepl("^KO member of Figure 4 feature:", derived_gene)] <-
  final_list$KO[grepl("^KO member of Figure 4 feature:", derived_gene)]
final_list$Gene_label[missing_gene] <- derived_gene[missing_gene]
final_list$Reference_basis <- ifelse(
  !is.na(final_list$Literature_reference) &
    nzchar(final_list$Literature_reference),
  paste0("KEGG Orthology; ", final_list$Literature_reference),
  "KEGG Orthology"
)
final_list$Detected_in_TPM <- vapply(final_list$KO, function(ko) {
  z <- ko_tpm[ko_tpm$KO == ko, samples, drop = FALSE]
  nrow(z) > 0 && sum(as.matrix(z), na.rm = TRUE) > 0
}, logical(1))
final_list$KEGG_URL <- ifelse(
  is.na(final_list$Source_URL) | !nzchar(final_list$Source_URL),
  paste0("https://www.kegg.jp/entry/", final_list$KO),
  final_list$Source_URL
)
final_list$Figure_scope_note <- "Included because the tile/pathway is present in final Figure 4"
final_list <- final_list[, c(
  "Panel", "Element", "Feature_level", "Section", "Figure_order",
  "Display_label", "KO", "Gene_label", "Function_description", "Module_map",
  "Detected_in_TPM", "Literature_reference", "Literature_pathway",
  "Literature_note", "Reference_basis", "KEGG_URL", "Abundance_source",
  "Figure_scope_note",
  "Feature_ID", "Member_ID"
)]

## Aggregate membership as an auditable companion table.
agg_members_out <- merge(agg_members,
                         aggregates[, c("Feature_ID", "Panel", "Element", "Section",
                                        "Display_label")],
                         by.x = "Aggregate_ID", by.y = "Feature_ID",
                         all.x = TRUE, sort = FALSE)
names(agg_members_out)[names(agg_members_out) == "Display_label"] <- "Aggregate_label"
agg_members_out <- merge(agg_members_out,
                         features[, c("Feature_ID", "Display_label")],
                         by.x = "Tile_ID", by.y = "Feature_ID",
                         all.x = TRUE, sort = FALSE)
names(agg_members_out)[names(agg_members_out) == "Display_label"] <- "Member_tile"
agg_members_out <- agg_members_out[, c("Panel", "Element", "Section",
                                       "Aggregate_ID", "Aggregate_label",
                                       "Tile_ID", "Member_tile")]

## Add total and KW columns to the detailed altitude table.
alt_long <- merge(alt_long, total_stats, by = "Feature_ID", all.x = TRUE,
                  suffixes = c("", "_Total"), sort = FALSE)
alt_long <- merge(alt_long, kw_stats, by = "Feature_ID", all.x = TRUE, sort = FALSE)
alt_long <- alt_long[order(match(alt_long$Feature_ID, all_feature_tpm$Feature_ID),
                           alt_long$Altitude_m), ]

write.table(final_list, file.path(out, "Final_Figure4_KO_list.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE, na = "")
write.table(features, file.path(out, "Figure4_tile_catalog.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE, na = "")
write.table(membership, file.path(out, "Figure4_tile_KO_membership.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE, na = "")
write.table(aggregates, file.path(out, "Figure4_aggregate_catalog.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE, na = "")
write.table(agg_members_out, file.path(out, "Figure4_aggregate_membership.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE, na = "")
write.table(all_feature_tpm, file.path(out, "Figure4_feature_TPM_by_sample.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE, na = "")
write.table(summary_tab, file.path(out, "Altitude_comparison_summary.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE, na = "")
write.table(alt_long, file.path(out, "Altitude_comparison_detailed.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE, na = "")

audit <- data.frame(
  Check = c(
    "Final Figure 4 PDF", "Figure tiles/pathways", "Aggregate rows",
    "Final list membership rows", "Unique KO count",
    "Undetected intended KOs", "Altitude correction",
    "Total correction", "Habitat mapping", "phnF identity"
  ),
  Result = c(
    "/Volumes/Haleo_DK01/WD_Takakia_F/Final_version/FIGURE 4(16).pdf",
    nrow(features), nrow(aggregates), nrow(final_list),
    length(unique(final_list$KO)),
    paste(final_list$KO[!final_list$Detected_in_TPM], collapse = ", "),
    "BH within each Figure feature across five altitudes",
    "BH across all submitted features for Total R/E",
    "S-prefix = Rhizomatousphere (R); P-prefix = Endophytic (E)",
    ifelse(any(final_list$KO == "K02043" & final_list$Display_label == "phnF"),
           "Confirmed: phnF = K02043", "FAILED")
  ),
  stringsAsFactors = FALSE
)
write.table(audit, file.path(out, "audit_summary.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE, na = "")

message("Prepared ", nrow(features), " Figure tiles/pathways, ",
        nrow(aggregates), " aggregate rows, ", nrow(final_list),
        " tile-KO membership rows and ", nrow(alt_long),
        " altitude-comparison rows.")
