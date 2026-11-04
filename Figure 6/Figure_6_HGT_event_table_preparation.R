options(stringsAsFactors = FALSE)

root <- Sys.getenv("TAKAKIA_PROJECT_DIR", unset = "~/Documents/Takakia")
proj <- file.path(root, "HGT_full_event_sheet_2026-07-27")
dir.create(file.path(proj, "work"), recursive = TRUE, showWarnings = FALSE)

raw_path <- "/Volumes/Haleo_DK01/WD_Takakia_F/all_g_detected_HGTs.txt"
joined_path <- "/Volumes/Haleo_DK01/WD_Takakia_F/metagenome-data/HGT_direction_COG_summary/HGT_events_with_MAG_group_and_COG.tsv"
donor_path <- "/Volumes/Haleo_DK01/WD_Takakia_F/metagenome-data/donor.tsv"
recipient_path <- "/Volumes/Haleo_DK01/WD_Takakia_F/metagenome-data/recipient.tsv"

raw <- read.delim(raw_path, check.names = FALSE, quote = "", comment.char = "")
hgt <- read.delim(joined_path, check.names = FALSE, quote = "", comment.char = "")
donor_anno <- read.delim(donor_path, check.names = FALSE, quote = "", comment.char = "")
recipient_anno <- read.delim(recipient_path, check.names = FALSE, quote = "", comment.char = "")
rownames(donor_anno) <- donor_anno$`#query`
rownames(recipient_anno) <- recipient_anno$`#query`

stopifnot(nrow(raw) == 452L, nrow(hgt) == 452L)
stopifnot(identical(raw$Gene_1, hgt$Gene_1))
stopifnot(identical(raw$Gene_2, hgt$Gene_2))
stopifnot(identical(raw$direction, hgt$direction))

clean_value <- function(x) {
  x[is.na(x) | !nzchar(trimws(x)) | x == "-"] <- "—"
  x
}

paired_cell <- function(donor, recipient) {
  donor <- clean_value(donor)
  recipient <- clean_value(recipient)
  ifelse(donor == recipient, donor, paste0("Donor: ", donor, "␤Recipient: ", recipient))
}

paired_by_gene <- function(donor, recipient, same_gene) {
  donor <- clean_value(donor)
  recipient <- clean_value(recipient)
  ifelse(same_gene, donor, paired_cell(donor, recipient))
}

strip_ko <- function(x) sub("^ko:", "", clean_value(x))

gene_annotation <- function(preferred, ko) {
  preferred <- clean_value(preferred)
  ko <- strip_ko(ko)
  if (preferred != "—" && ko != "—") return(paste0(preferred, " (", ko, ")"))
  if (preferred != "—") return(preferred)
  if (ko != "—") return(ko)
  "Unannotated"
}

cog_map <- c(
  A="RNA processing and modification", B="Chromatin structure and dynamics",
  C="Energy production and conversion", D="Cell cycle control, cell division, chromosome partitioning",
  E="Amino acid transport and metabolism", F="Nucleotide transport and metabolism",
  G="Carbohydrate transport and metabolism", H="Coenzyme transport and metabolism",
  I="Lipid transport and metabolism", J="Translation, ribosomal structure and biogenesis",
  K="Transcription", L="Replication, recombination and repair",
  M="Cell wall/membrane/envelope biogenesis", N="Cell motility",
  O="Posttranslational modification, protein turnover, and chaperones",
  P="Inorganic ion transport and metabolism", Q="Secondary metabolites biosynthesis, transport and catabolism",
  R="General function prediction only", S="Function unknown", T="Signal transduction mechanisms",
  U="Intracellular trafficking, secretion, and vesicular transport", V="Defense mechanisms",
  W="Extracellular structures", X="Mobilome: prophages and transposons",
  Y="Nuclear structure", Z="Cytoskeleton"
)

expand_cog <- function(x) {
  x <- clean_value(x)
  if (x == "—") return(x)
  letters <- strsplit(gsub("[^A-Z]", "", x), "")[[1]]
  paste(vapply(letters, function(letter) paste0(letter, " — ", ifelse(letter %in% names(cog_map), cog_map[[letter]], "Unclassified")), ""), collapse="; ")
}

cog_letters <- function(x) {
  x <- clean_value(x)
  if (x == "—" || !grepl("^[A-Z]+$", x)) return("—")
  x
}

anno_field <- function(tab, gene, field) {
  if (!gene %in% rownames(tab)) return("—")
  clean_value(tab[gene, field])
}

database_annotation <- function(tab, gene, description, ko) {
  description <- clean_value(description)
  ko <- strip_ko(ko)
  eggnog_ogs <- anno_field(tab, gene, "eggNOG_OGs")
  ec <- anno_field(tab, gene, "EC")
  cazy <- anno_field(tab, gene, "CAZy")
  pfam <- anno_field(tab, gene, "PFAMs")
  parts <- character(0)
  if (eggnog_ogs != "—") parts <- c(parts, paste0("EggNOG OGs: ", eggnog_ogs))
  if (description != "—") parts <- c(parts, paste0("Description: ", description))
  if (ko != "—") parts <- c(parts, paste0("KEGG KO: ", ko))
  if (ec != "—") parts <- c(parts, paste0("EC: ", ec))
  if (cazy != "—") parts <- c(parts, paste0("CAZy: ", cazy))
  if (!length(parts) && pfam != "—") parts <- c(parts, paste0("PFAM: ", pfam))
  if (!length(parts)) "Unannotated" else paste(parts, collapse="; ")
}

resolve_direction_gene <- function(gene_1, gene_2, bin_id, fallback) {
  prefix <- paste0(bin_id, "_")
  if (startsWith(gene_1, prefix)) return(gene_1)
  if (startsWith(gene_2, prefix)) return(gene_2)
  fallback <- clean_value(fallback)
  if (fallback != "—") return(fallback)
  "—"
}

donor_gene_full <- mapply(resolve_direction_gene, hgt$Gene_1, hgt$Gene_2, hgt$donor_bin, hgt$donor_gene, USE.NAMES = FALSE)
recipient_gene_full <- mapply(resolve_direction_gene, hgt$Gene_1, hgt$Gene_2, hgt$recipient_bin, hgt$recipient_gene, USE.NAMES = FALSE)
stopifnot(all(donor_gene_full != "—"), all(recipient_gene_full != "—"))

prefer_direct <- function(direct, fallback) {
  direct <- clean_value(direct)
  fallback <- clean_value(fallback)
  ifelse(direct != "—", direct, fallback)
}

donor_preferred <- mapply(function(gene, fallback) prefer_direct(anno_field(donor_anno, gene, "Preferred_name"), fallback), donor_gene_full, hgt$donor_Preferred_name, USE.NAMES = FALSE)
recipient_preferred <- mapply(function(gene, fallback) prefer_direct(anno_field(recipient_anno, gene, "Preferred_name"), fallback), recipient_gene_full, hgt$recipient_Preferred_name, USE.NAMES = FALSE)
donor_ko <- mapply(function(gene, fallback) prefer_direct(anno_field(donor_anno, gene, "KEGG_ko"), fallback), donor_gene_full, hgt$donor_KEGG_ko, USE.NAMES = FALSE)
recipient_ko <- mapply(function(gene, fallback) prefer_direct(anno_field(recipient_anno, gene, "KEGG_ko"), fallback), recipient_gene_full, hgt$recipient_KEGG_ko, USE.NAMES = FALSE)
donor_description <- mapply(function(gene, fallback) prefer_direct(anno_field(donor_anno, gene, "Description"), fallback), donor_gene_full, hgt$donor_Description, USE.NAMES = FALSE)
recipient_description <- mapply(function(gene, fallback) prefer_direct(anno_field(recipient_anno, gene, "Description"), fallback), recipient_gene_full, hgt$recipient_Description, USE.NAMES = FALSE)
donor_cog <- mapply(function(gene, fallback) prefer_direct(anno_field(donor_anno, gene, "COG_category"), fallback), donor_gene_full, hgt$donor_COG_category, USE.NAMES = FALSE)
recipient_cog <- mapply(function(gene, fallback) prefer_direct(anno_field(recipient_anno, gene, "COG_category"), fallback), recipient_gene_full, hgt$recipient_COG_category, USE.NAMES = FALSE)

donor_gene_annotation <- mapply(gene_annotation, donor_preferred, donor_ko, USE.NAMES=FALSE)
recipient_gene_annotation <- mapply(gene_annotation, recipient_preferred, recipient_ko, USE.NAMES=FALSE)
donor_ko_clean <- strip_ko(donor_ko)
recipient_ko_clean <- strip_ko(recipient_ko)
same_gene <- (donor_ko_clean != "—" & donor_ko_clean == recipient_ko_clean) |
  (clean_value(donor_preferred) != "—" & clean_value(donor_preferred) == clean_value(recipient_preferred))
same_gene_display <- ifelse(same_gene & clean_value(recipient_preferred) != "—", recipient_gene_annotation, donor_gene_annotation)
donor_db <- mapply(function(gene, desc, ko) database_annotation(donor_anno, gene, desc, ko), donor_gene_full, donor_description, donor_ko, USE.NAMES=FALSE)
recipient_db <- mapply(function(gene, desc, ko) database_annotation(recipient_anno, gene, desc, ko), recipient_gene_full, recipient_description, recipient_ko, USE.NAMES=FALSE)

out <- data.frame(
  "HGT event" = paste0(donor_gene_full, " → ", recipient_gene_full),
  "Identity (%)" = round(hgt$Identity, 3),
  "Transfer direction" = hgt$direction_group,
  "Transferred gene" = ifelse(same_gene, same_gene_display, paired_cell(donor_gene_annotation, recipient_gene_annotation)),
  "COG category" = paired_by_gene(vapply(donor_cog, cog_letters, ""), vapply(recipient_cog, cog_letters, ""), same_gene),
  "Database description" = paired_by_gene(donor_db, recipient_db, same_gene),
  check.names = FALSE
)

direction_order <- c("R-R", "E-E", "R-E", "E-R")
event_number <- as.integer(sub("HGT_", "", hgt$event_id))
ord <- order(match(out[["Transfer direction"]], direction_order), event_number)
out <- out[ord, , drop = FALSE]

counts <- table(factor(out[["Transfer direction"]], levels = direction_order))
stopifnot(identical(as.integer(counts), c(284L, 109L, 37L, 22L)))

write.table(out, file.path(proj, "work/HGT_all_452_events_display.tsv"), sep = "\t", quote = FALSE, row.names = FALSE, na = "")

audit <- data.frame(
  Item = c("Total HGT events", direction_order, "Different donor/recipient gene annotations", "Different donor/recipient COG categories", "Different donor/recipient database descriptions"),
  Value = c(nrow(out), as.integer(counts),
            sum(!same_gene),
            sum(vapply(donor_cog, expand_cog, "") != vapply(recipient_cog, expand_cog, "")),
            sum(donor_db != recipient_db))
)
write.table(audit, file.path(proj, "work/HGT_all_452_events_audit.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
print(audit)
