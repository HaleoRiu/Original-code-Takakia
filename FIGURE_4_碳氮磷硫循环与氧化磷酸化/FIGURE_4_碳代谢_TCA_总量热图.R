options(stringsAsFactors = FALSE, warn = 1)
suppressPackageStartupMessages(library(grid))

script_arg <- grep("^--file=", commandArgs(), value = TRUE)
script_dir <- if (length(script_arg)) dirname(normalizePath(sub("^--file=", "", script_arg[1]))) else getwd()
project_dir <- normalizePath(file.path(script_dir, ".."), mustWork = FALSE)
variant <- Sys.getenv("C_MAIN_VARIANT", unset = "V1_rTCA")
if (!variant %in% c("V1_rTCA", "V2_TCA")) stop("C_MAIN_VARIANT must be V1_rTCA or V2_TCA")
work <- file.path(project_dir, "work", variant)

labels <- read.delim(file.path(work, "C_label_catalog.tsv"), check.names = FALSE)
sections <- read.delim(file.path(work, "C_section_catalog.tsv"), check.names = FALSE)
label_heat <- read.delim(file.path(work, "C_label_heat_values.tsv"), check.names = FALSE)
label_total <- read.delim(file.path(work, "C_label_paired_Wilcoxon_Total.tsv"), check.names = FALSE)
section_heat <- read.delim(file.path(work, "C_section_heat_values.tsv"), check.names = FALSE)
section_total <- read.delim(file.path(work, "C_section_paired_Wilcoxon_Total.tsv"), check.names = FALSE)
scale_info <- read.delim(file.path(work, "shared_fill_scale.tsv"), check.names = FALSE)

label_heat <- label_heat[label_heat$Scope == "Total", , drop = FALSE]
section_heat <- section_heat[section_heat$Scope == "Total", , drop = FALSE]
stopifnot(nrow(label_heat) == 2L * nrow(labels), nrow(section_heat) == 2L * nrow(sections))

cell_pt <- 32.33
font_pt <- 28
section_gap <- 40
left_margin <- 105
right_margin <- 105
bottom_margin <- 55
top_margin <- 80
y_line <- 190
y_section_label <- 178
y_circle <- 270
y_circle_key <- 235
y_e <- 365
y_r <- y_e + cell_pt
y_label_anchor <- y_r + cell_pt + 16
page_height_pt <- 960
fill_max <- scale_info$Fixed_maximum[1]
palette <- colorRampPalette(c("#2C7BB6", "#F7F7F7", "#D7191C"))(512)

fill_for <- function(z) {
  if (!is.finite(z)) z <- 0
  idx <- pmax(1, pmin(length(palette), floor((z / fill_max) * (length(palette) - 1)) + 1))
  palette[idx]
}

label_sections <- unique(labels$Section)
section_gap_before <- setNames((seq_along(label_sections) - 1) * section_gap, label_sections)
labels$x_left <- left_margin + (seq_len(nrow(labels)) - 1) * cell_pt + section_gap_before[labels$Section]
page_width_pt <- max(labels$x_left) + cell_pt + right_margin

section_span <- do.call(rbind, lapply(seq_len(nrow(sections)), function(i) {
  z <- labels[labels$Section == sections$Section[i], , drop = FALSE]
  if (!nrow(z)) stop("Section has no labels: ", sections$Section[i])
  data.frame(
    Feature_ID = sections$Feature_ID[i],
    Section = sections$Section[i],
    x_start = min(z$x_left),
    x_end = max(z$x_left) + cell_pt,
    x_mid = (min(z$x_left) + max(z$x_left) + cell_pt) / 2,
    stringsAsFactors = FALSE
  )
}))

draw_text <- function(txt, x, y, just = "center", rot = 0, face = "plain", size = font_pt) {
  txt <- gsub("—", "-", txt, fixed = TRUE)
  txt <- gsub("–", "-", txt, fixed = TRUE)
  grid.text(txt, x = unit(x, "pt"), y = unit(y, "pt"), just = just, rot = rot,
            gp = gpar(fontfamily = "Helvetica", fontface = face, fontsize = size, col = "black"))
}

heat_value <- function(df, fid, habitat) {
  z <- df[df$Feature_ID == fid & df$Habitat == habitat, , drop = FALSE]
  if (nrow(z) != 1L) stop("Expected one heat row for ", fid, " / ", habitat)
  z[1, ]
}

stat_value <- function(df, fid) {
  z <- df[df$Feature_ID == fid, , drop = FALSE]
  if (nrow(z) != 1L) stop("Expected one total-stat row for ", fid)
  z[1, ]
}

star_color <- function(fill) {
  rgb <- col2rgb(fill) / 255
  lum <- 0.2126 * rgb[1] + 0.7152 * rgb[2] + 0.0722 * rgb[3]
  if (lum < 0.48) "white" else "black"
}

draw_tile <- function(x, y, rec, shape = c("square", "circle"), star = "") {
  shape <- match.arg(shape)
  fill <- fill_for(rec$log2_mean_TPM_plus1)
  if (shape == "square") {
    grid.rect(x = unit(x, "pt"), y = unit(y, "pt"), width = unit(cell_pt, "pt"), height = unit(cell_pt, "pt"),
              just = c("left", "bottom"), gp = gpar(fill = fill, col = "white", lwd = 0.65))
    cx <- x + cell_pt / 2
    cy <- y + cell_pt / 2
  } else {
    cx <- x
    cy <- y
    grid.circle(x = unit(cx, "pt"), y = unit(cy, "pt"), r = unit(cell_pt / 2, "pt"),
                gp = gpar(fill = fill, col = "white", lwd = 0.65))
  }
  if (isTRUE(all.equal(rec$Mean_TPM, 0))) {
    draw_text("X", cx, cy, face = "bold")
  }
  if (nzchar(star)) {
    grid.text(star, x = unit(cx, "pt"), y = unit(cy, "pt"),
              gp = gpar(fontfamily = "Helvetica", fontface = "bold", fontsize = 20,
                        col = star_color(fill)))
  }
}

draw_legend <- function() {
  draw_text("R - Rhizomatousphere    E - Endophytic", left_margin, page_height_pt - 35,
            just = c("left", "center"))
  legend_w <- 360
  legend_h <- 22
  lx <- page_width_pt - right_margin - legend_w
  ly <- page_height_pt - 90
  for (i in seq_len(180)) {
    x <- lx + (i - 1) * legend_w / 180
    grid.rect(x = unit(x, "pt"), y = unit(ly, "pt"), width = unit(legend_w / 180 + 0.2, "pt"),
              height = unit(legend_h, "pt"), just = c("left", "bottom"),
              gp = gpar(fill = fill_for((i - 1) / 179 * fill_max), col = NA))
  }
  draw_text("log2(mean TPM + 1)", lx + legend_w / 2, ly + 54, face = "bold")
  for (tick in seq(0, fill_max, by = 3)) {
    tx <- lx + legend_w * tick / fill_max
    grid.segments(unit(tx, "pt"), unit(ly, "pt"), unit(tx, "pt"), unit(ly - 8, "pt"),
                  gp = gpar(col = "black", lwd = 0.8))
    draw_text(as.character(tick), tx, ly - 24)
  }
}

draw_page <- function() {
  grid.newpage()
  pushViewport(viewport(xscale = c(0, page_width_pt), yscale = c(0, page_height_pt)))
  grid.rect(gp = gpar(fill = "white", col = NA))
  draw_legend()

  draw_text("R", left_margin - 18, y_r + cell_pt / 2, just = c("right", "center"), face = "bold")
  draw_text("E", left_margin - 18, y_e + cell_pt / 2, just = c("right", "center"), face = "bold")

  for (i in seq_len(nrow(labels))) {
    fid <- labels$Feature_ID[i]
    st <- stat_value(label_total, fid)
    for (hab in c("R", "E")) {
      rec <- heat_value(label_heat, fid, hab)
      yy <- if (hab == "R") y_r else y_e
      star <- if (st$Higher == hab) st$Significance else ""
      draw_tile(labels$x_left[i], yy, rec, "square", star)
    }
    draw_text(labels$Display_label[i], labels$x_left[i] + cell_pt / 2, y_label_anchor,
              just = c("left", "center"), rot = 90, size = 24)
  }

  for (i in seq_len(nrow(section_span))) {
    ss <- section_span[i, ]
    fid <- ss$Feature_ID
    st <- stat_value(section_total, fid)
    cx_r <- ss$x_mid - cell_pt * 0.58
    cx_e <- ss$x_mid + cell_pt * 0.58
    rec_r <- heat_value(section_heat, fid, "R")
    rec_e <- heat_value(section_heat, fid, "E")
    draw_tile(cx_r, y_circle, rec_r, "circle", if (st$Higher == "R") st$Significance else "")
    draw_tile(cx_e, y_circle, rec_e, "circle", if (st$Higher == "E") st$Significance else "")
    draw_text("R", cx_r, y_circle_key, face = "bold")
    draw_text("E", cx_e, y_circle_key, face = "bold")
    grid.segments(unit(ss$x_start, "pt"), unit(y_line, "pt"), unit(ss$x_end, "pt"), unit(y_line, "pt"),
                  gp = gpar(col = "#4D4D4D", lwd = 1.25))
    section_label <- if (ss$Section == "rTCA-related machinery (incomplete)")
      "rTCA-related (incomplete)" else ss$Section
    draw_text(section_label, ss$x_mid, y_section_label, just = c("right", "center"), rot = 35, size = 20)
  }
  popViewport()
}

pdf_name <- if (variant == "V1_rTCA")
  "C_maintext_V1_rTCA_total_only.pdf" else "C_maintext_V2_TCA_total_only.pdf"
pdf_path <- file.path(work, pdf_name)
pdf(pdf_path, width = page_width_pt / 72, height = page_height_pt / 72,
    family = "Helvetica", useDingbats = FALSE, onefile = TRUE)
draw_page()
dev.off()

write.table(data.frame(
  Item = c("PDF_width_pt", "PDF_height_pt", "Square_width_pt", "Square_height_pt",
           "Circle_diameter_pt", "Font_size_pt", "Shared_fill_min", "Shared_fill_max",
           "Displayed_label_squares", "Displayed_section_circle_pairs", "Zero_mark", "Test"),
  Value = c(page_width_pt, page_height_pt, cell_pt, cell_pt, cell_pt, font_pt, 0, fill_max,
            nrow(labels), nrow(sections), "X", "Paired Wilcoxon on raw TPM; BH within square or circle block")
), file.path(work, "maintext_total_figure_geometry.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

cat(sprintf("Main-text Total-only C panel: %s\n", pdf_path))
cat(sprintf("Figure: %.2f x %.2f pt; %d squares; %d pathway circle pairs; cell %.2f pt; font %.1f pt\n",
            page_width_pt, page_height_pt, nrow(labels), nrow(sections), cell_pt, font_pt))
