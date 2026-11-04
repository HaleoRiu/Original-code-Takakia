#!/usr/bin/env Rscript
# Within-group PERMANOVA (adonis2): does ALTITUDE affect community structure
# separately within the Rhizomatousphere (R) and Endophytic (E) groups?
# Datasets = the four ordination (PCoA) inputs:
#   Amplicon (Figure 1): bacteria ASV, fungi ASV
#   Metagenome (Figure 3): MAG abundance, functional-gene (KO) TPM
# Distance = Bray-Curtis (matches the PCoA plots). Altitude treated as a
# 5-level categorical factor (3800/3900/4000/4100/4200 m); 6 reps per level.

suppressPackageStartupMessages({ library(vegan) })

amp <- "/Volumes/Haleo_DK01/WD_Takakia_F/amplicon-data"
mg  <- "/Volumes/Haleo_DK01/WD_Takakia_F/metagenome-data"
out_csv <- "/Volumes/Haleo_DK01/WD_Takakia_F/adonis_within_group_altitude_results.csv"
out_txt <- "/Volumes/Haleo_DK01/WD_Takakia_F/adonis_within_group_altitude_summary.txt"

read_comm <- function(path) {
  raw <- read.delim(path, header = TRUE, row.names = 1, check.names = FALSE,
                    stringsAsFactors = FALSE, comment.char = "")
  scols <- grep("^[SP][0-9]{2}_[0-9]+$", colnames(raw), value = TRUE)
  comm <- t(as.matrix(raw[, scols, drop = FALSE]))
  storage.mode(comm) <- "numeric"
  comm[is.na(comm)] <- 0
  comm
}

meta_of <- function(samp) {
  data.frame(
    Sample   = samp,
    Group    = ifelse(grepl("^S", samp), "R",
                ifelse(grepl("^P", samp), "E", NA_character_)),
    Altitude = factor(as.numeric(substr(sub("^[SP]", "", samp), 1, 2)) * 100),
    stringsAsFactors = FALSE
  )
}

# pick the functional-gene KO/TPM table: first candidate with exactly 60 S/P cols
pick_func_table <- function() {
  cands <- c("tpm_renormalized.txt", "filtered_ko_1.txt", "div.txt")
  for (f in cands) {
    p <- file.path(mg, f)
    if (!file.exists(p)) next
    hd <- read.delim(p, header = TRUE, nrows = 1, check.names = FALSE, comment.char = "")
    n  <- length(grep("^[SP][0-9]{2}_[0-9]+$", colnames(hd)))
    if (n == 60) return(p)
  }
  stop("No functional-gene table with 60 S/P sample columns found among: ",
       paste(cands, collapse = ", "))
}

datasets <- list(
  list(name = "Amplicon_Bacteria(Fig1)", path = file.path(amp, "asv_bacteria.txt")),
  list(name = "Amplicon_Fungi(Fig1)",    path = file.path(amp, "asv_fungi.txt")),
  list(name = "Metagenome_MAG(Fig3)",    path = file.path(mg,  "MAGs_with_taxonomy_70.txt")),
  list(name = "Metagenome_FuncGene_KO(Fig3)", path = pick_func_table())
)

run_one <- function(name, path) {
  comm <- read_comm(path)
  meta <- meta_of(rownames(comm))
  res <- list()
  for (g in c("R", "E")) {
    idx <- meta$Group == g
    cm  <- comm[meta$Sample[idx], , drop = FALSE]
    md  <- meta[idx, , drop = FALSE]
    md$Altitude <- droplevels(md$Altitude)
    cm  <- cm[, colSums(cm) > 0, drop = FALSE]        # drop features absent in subset
    d   <- vegdist(cm, method = "bray")
    set.seed(20260712)
    ad  <- adonis2(d ~ Altitude, data = md, permutations = 999)
    res[[g]] <- data.frame(
      dataset   = name,
      group     = ifelse(g == "R", "R (Rhizomatousphere)", "E (Endophytic)"),
      n_samples = nrow(cm),
      n_features= ncol(cm),
      altitudes = nlevels(md$Altitude),
      Df        = ad$Df[1],
      SumOfSqs  = round(ad$SumOfSqs[1], 4),
      R2        = round(ad$R2[1], 4),
      F         = round(ad$F[1], 4),
      p_value   = ad$`Pr(>F)`[1],
      stringsAsFactors = FALSE
    )
  }
  do.call(rbind, res)
}

all_res <- do.call(rbind, lapply(datasets, function(d) run_one(d$name, d$path)))
# BH-FDR across the 8 tests
all_res$p_adj_BH <- round(p.adjust(all_res$p_value, method = "BH"), 4)
rownames(all_res) <- NULL

write.csv(all_res, out_csv, row.names = FALSE)

sink(out_txt)
cat("Within-group PERMANOVA (adonis2, Bray-Curtis, 999 perm)\n")
cat("H0: altitude has no effect on community structure within the group\n")
cat("Altitude = 5-level factor (3800-4200 m), 6 reps each\n")
cat("Functional-gene table used: ", basename(datasets[[4]]$path), "\n\n", sep = "")
print(all_res, row.names = FALSE)
sink()

cat("=== DONE ===\n")
cat("Functional-gene table:", basename(datasets[[4]]$path), "\n\n")
print(all_res, row.names = FALSE)
