options(stringsAsFactors = FALSE)
library(grid)

work <- file.path(Sys.getenv("TAKAKIA_PROJECT_DIR", unset = "~/Documents/Takakia"), "outputs/CNPS_altitude_sum_revision_2026-07-21/V2_3_categories/work")
upper <- read.delim(file.path(work, "upper_feature_catalog.tsv"), check.names = FALSE)
summary_cat <- read.delim(file.path(work, "pathway_summary_catalog.tsv"), check.names = FALSE)
upper_heat <- read.delim(file.path(work, "upper_heat_values.tsv"), check.names = FALSE)
upper_alt <- read.delim(file.path(work, "upper_paired_Wilcoxon_by_altitude.tsv"), check.names = FALSE)
upper_total <- read.delim(file.path(work, "upper_paired_Wilcoxon_Total.tsv"), check.names = FALSE)
upper_kw <- read.delim(file.path(work, "upper_KW_altitude.tsv"), check.names = FALSE)
sum_heat <- read.delim(file.path(work, "summary_heat_values.tsv"), check.names = FALSE)
sum_alt <- read.delim(file.path(work, "summary_paired_Wilcoxon_by_altitude.tsv"), check.names = FALSE)
sum_total <- read.delim(file.path(work, "summary_paired_Wilcoxon_Total.tsv"), check.names = FALSE)
sum_kw <- read.delim(file.path(work, "summary_KW_altitude.tsv"), check.names = FALSE)
scale_info <- read.delim(file.path(work, "shared_fill_scale.tsv"), check.names = FALSE)

cell_pt <- 32.33
font_pt <- 28
line_col <- "#4D4D4D"
alts <- c(3800, 3900, 4000, 4100, 4200)
habitats <- c("R", "E")
fill_max <- scale_info$Maximum[1]
palette <- colorRampPalette(c("#2C7BB6", "#F7F7F7", "#D7191C"))(512)
fill_for <- function(z) {
  idx <- pmax(1, pmin(length(palette), floor((z / fill_max) * (length(palette) - 1)) + 1))
  palette[idx]
}

label_widths <- c(col1 = 470, col2 = 540, col3 = 700)
kw_col_pt <- 48
chart_width <- 10 * cell_pt + 16 + 2 * kw_col_pt + 10 + 2 * cell_pt
col_gap <- 52
left_margin <- 34
right_margin <- 34
top_margin <- 170
bottom_margin <- 60

block_rows <- function(element) {
  1.45 + 2.55 + sum(upper$Element == element) + 0.75 + 1.15 + sum(summary_cat$Element == element) + 1.35
}
col1_rows <- block_rows("Carbon") + block_rows("Nitrogen")
col2_rows <- block_rows("Phosphorus") + block_rows("Sulfur")
col3_rows <- block_rows("Oxidative phosphorylation")
page_height_pt <- top_margin + max(col1_rows, col2_rows, col3_rows) * cell_pt + bottom_margin
page_width_pt <- left_margin + sum(label_widths + chart_width) + 2 * col_gap + right_margin

col_x <- c(
  col1 = left_margin,
  col2 = left_margin + unname(label_widths[1]) + chart_width + col_gap,
  col3 = left_margin + unname(label_widths[1]) + chart_width + col_gap + unname(label_widths[2]) + chart_width + col_gap
)

draw_text <- function(label, x, y, just = "centre", face = "plain", size = font_pt, rot = 0, col = "black") {
  grid.text(label, x = unit(x, "native"), y = unit(y, "native"), just = just, rot = rot,
            gp = gpar(fontfamily = "Helvetica", fontface = face, fontsize = size, col = col))
}

draw_header <- function(x0, y0, label_w) {
  hx <- x0 + label_w
  draw_text("Feature", hx - 14, y0 + 1.25 * cell_pt, just = "right", face = "bold")
  for (i in seq_along(alts)) {
    cx <- hx + ((i - 1) * 2 + 1) * cell_pt
    draw_text(as.character(alts[i]), cx, y0 + 0.75 * cell_pt, face = "bold", rot = 45)
    draw_text("R", cx - 0.5 * cell_pt, y0 + 2.00 * cell_pt, face = "bold")
    draw_text("E", cx + 0.5 * cell_pt, y0 + 2.00 * cell_pt, face = "bold")
  }
  kw0 <- hx + 10 * cell_pt + 16
  draw_text("KW", kw0 + kw_col_pt, y0 + 0.75 * cell_pt, face = "bold")
  draw_text("R", kw0 + 0.5 * kw_col_pt, y0 + 2.00 * cell_pt, face = "bold")
  draw_text("E", kw0 + 1.5 * kw_col_pt, y0 + 2.00 * cell_pt, face = "bold")
  total0 <- kw0 + 2 * kw_col_pt + 10
  draw_text("Total", total0 + cell_pt, y0 + 0.75 * cell_pt, face = "bold")
  draw_text("R", total0 + 0.5 * cell_pt, y0 + 2.00 * cell_pt, face = "bold")
  draw_text("E", total0 + 1.5 * cell_pt, y0 + 2.00 * cell_pt, face = "bold")
}

heat_value <- function(heat_df, fid, scope, altitude, habitat) {
  hit <- heat_df$Feature_ID == fid & heat_df$Scope == scope & heat_df$Habitat == habitat
  if (scope == "Altitude") hit <- hit & heat_df$Altitude == altitude
  heat_df$log2_mean_TPM_plus1[hit][1]
}

draw_row <- function(fid, label, element, y, x0, label_w, shape,
                     heat_df, alt_df, total_df, kw_df) {
  hx <- x0 + label_w
  is_summary_row <- fid %in% summary_cat$Feature_ID
  face <- if (shape == "square" && !is_summary_row && element %in% c("Nitrogen", "Phosphorus", "Sulfur")) "italic" else "plain"
  draw_text(label, hx - 14, y, just = "right", face = face)

  for (i in seq_along(alts)) for (j in seq_along(habitats)) {
    cx <- hx + ((i - 1) * 2 + j - 0.5) * cell_pt
    z <- heat_value(heat_df, fid, "Altitude", alts[i], habitats[j])
    if (shape == "square") {
      grid.rect(x = unit(cx, "native"), y = unit(y, "native"),
                width = unit(cell_pt, "pt"), height = unit(cell_pt, "pt"),
                gp = gpar(fill = fill_for(z), col = "white", lwd = 0.45))
    } else {
      grid.circle(x = unit(cx, "native"), y = unit(y, "native"),
                  r = unit(cell_pt / 2, "pt"),
                  gp = gpar(fill = fill_for(z), col = line_col, lwd = 0.65))
    }
  }

  kw0 <- hx + 10 * cell_pt + 16
  for (j in seq_along(habitats)) {
    sig <- kw_df$Significance[kw_df$Feature_ID == fid & kw_df$Habitat == habitats[j]][1]
    draw_text(sig, kw0 + (j - 0.5) * kw_col_pt, y, face = "plain")
  }

  total0 <- kw0 + 2 * kw_col_pt + 10
  for (j in seq_along(habitats)) {
    cx <- total0 + (j - 0.5) * cell_pt
    z <- heat_value(heat_df, fid, "Total", NA_integer_, habitats[j])
    if (shape == "square") {
      grid.rect(x = unit(cx, "native"), y = unit(y, "native"),
                width = unit(cell_pt, "pt"), height = unit(cell_pt, "pt"),
                gp = gpar(fill = fill_for(z), col = "white", lwd = 0.45))
    } else {
      grid.circle(x = unit(cx, "native"), y = unit(y, "native"),
                  r = unit(cell_pt / 2, "pt"),
                  gp = gpar(fill = fill_for(z), col = line_col, lwd = 0.65))
    }
  }

  # Paired-Wilcoxon stars are placed on the higher R/E cell.
  aa <- alt_df[alt_df$Feature_ID == fid, ]
  for (i in seq_along(alts)) {
    arow <- aa[aa$Altitude == alts[i], ]
    if (nrow(arow) && nzchar(arow$Significance[1])) {
      j <- match(arow$Higher[1], habitats)
      cx <- hx + ((i - 1) * 2 + j - 0.5) * cell_pt
      draw_text(arow$Significance[1], cx, y, face = "plain")
    }
  }
  tt <- total_df[total_df$Feature_ID == fid, ]
  if (nrow(tt) && nzchar(tt$Significance[1])) {
    j <- match(tt$Higher[1], habitats)
    draw_text(tt$Significance[1], total0 + (j - 0.5) * cell_pt, y, face = "plain")
  }
}

draw_separators <- function(d, x0, y_start, label_w) {
  if (nrow(d) < 2) return(invisible(NULL))
  breaks <- which(d$Section[-1] != d$Section[-nrow(d)])
  for (b in breaks) {
    yy <- y_start + b * cell_pt
    grid.lines(x = unit(c(x0 + 4, x0 + label_w + chart_width), "native"),
               y = unit(c(yy, yy), "native"), gp = gpar(col = line_col, lwd = 1.05))
  }
}

draw_element <- function(element, x0, y0, label_w) {
  u <- upper[upper$Element == element, ]
  u <- u[order(u$Feature_order), ]
  s <- summary_cat[summary_cat$Element == element, ]
  s <- s[order(s$Feature_order), ]

  title_y <- y0 + 0.58 * cell_pt
  draw_text(element, x0, title_y, just = "left", face = "bold")
  if (element == "Oxidative phosphorylation") {
    draw_text("all 21 Pathways retained; detected-KO direct sums", x0 + label_w + chart_width, title_y,
              just = "right", face = "plain", size = font_pt)
  }
  header_y <- y0 + 1.45 * cell_pt
  draw_header(x0, header_y, label_w)
  u_start <- header_y + 2.55 * cell_pt
  for (i in seq_len(nrow(u))) {
    draw_row(u$Feature_ID[i], u$Display_label[i], element,
             u_start + (i - 0.5) * cell_pt, x0, label_w, "square",
             upper_heat, upper_alt, upper_total, upper_kw)
  }
  draw_separators(u, x0, u_start, label_w)

  sum_label_y <- u_start + nrow(u) * cell_pt + 0.75 * cell_pt
  grid.lines(x = unit(c(x0 + 4, x0 + label_w + chart_width), "native"),
             y = unit(c(sum_label_y - 0.35 * cell_pt, sum_label_y - 0.35 * cell_pt), "native"),
             gp = gpar(col = line_col, lwd = 1.5))
  draw_text("Pathway sums", x0, sum_label_y, just = "left", face = "bold")
  s_start <- sum_label_y + 1.15 * cell_pt
  for (i in seq_len(nrow(s))) {
    summary_shape <- if (s$Feature_ID[i] == "P06") "square" else "circle"
    draw_row(s$Feature_ID[i], s$Display_label[i], element,
             s_start + (i - 0.5) * cell_pt, x0, label_w, summary_shape,
             sum_heat, sum_alt, sum_total, sum_kw)
  }
  s_start + nrow(s) * cell_pt + 1.35 * cell_pt
}

draw_legend <- function() {
  draw_text("Altitude-resolved CNPS, sulfur-transport and oxidative-phosphorylation functions", left_margin, 46,
            just = "left", face = "bold", size = font_pt)
  draw_text("Squares: feature/enzyme/Pathway rows and phn operon sum   Circles: other summed small Pathways", left_margin, 88,
            just = "left", size = font_pt)
  lx <- page_width_pt - 720
  ly <- 80
  legend_w <- 520
  n <- 130
  for (i in seq_len(n)) {
    x <- lx + (i - 0.5) * legend_w / n
    grid.rect(x = unit(x, "native"), y = unit(ly, "native"),
              width = unit(legend_w / n + 0.2, "pt"), height = unit(28, "pt"),
              gp = gpar(fill = palette[round(1 + (i - 1) * (length(palette) - 1) / (n - 1))], col = NA))
  }
  draw_text("log2(mean TPM + 1)", lx + legend_w / 2, ly - 35, face = "bold")
  tick_values <- seq(0, fill_max, by = 3)
  if (tail(tick_values, 1) != fill_max) tick_values <- c(tick_values, fill_max)
  for (tv in tick_values) {
    tx <- lx + (tv / fill_max) * legend_w
    grid.lines(x = unit(c(tx, tx), "native"), y = unit(c(ly + 14, ly + 23), "native"),
               gp = gpar(col = line_col, lwd = 1.0))
    draw_text(format(tv, trim = TRUE, scientific = FALSE), tx, ly + 48, just = "centre")
  }
}

draw_page <- function() {
  grid.newpage()
  pushViewport(viewport(xscale = c(0, page_width_pt), yscale = c(page_height_pt, 0)))
  grid.rect(gp = gpar(fill = "white", col = NA))
  draw_legend()
  y1 <- top_margin
  y1 <- draw_element("Carbon", col_x["col1"], y1, label_widths["col1"])
  draw_element("Nitrogen", col_x["col1"], y1, label_widths["col1"])
  y2 <- top_margin
  y2 <- draw_element("Phosphorus", col_x["col2"], y2, label_widths["col2"])
  draw_element("Sulfur", col_x["col2"], y2, label_widths["col2"])
  draw_element("Oxidative phosphorylation", col_x["col3"], top_margin, label_widths["col3"])
  popViewport()
}

pdf(file.path(work, "CNPS_altitude_sum_revision_V2_3categories_28pt.pdf"),
    width = page_width_pt / 72, height = page_height_pt / 72,
    family = "Helvetica", useDingbats = FALSE, onefile = TRUE)
draw_page()
dev.off()

png(file.path(work, "CNPS_altitude_sum_revision_V2_3categories_28pt.png"),
    width = page_width_pt / 72 * 120, height = page_height_pt / 72 * 120,
    res = 120, type = "cairo", bg = "white")
draw_page()
dev.off()

write.table(data.frame(
  Item = c("PDF_width_pt", "PDF_height_pt", "Square_width_pt", "Square_height_pt",
           "Circle_diameter_pt", "Font_size_pt", "Shared_fill_min", "Shared_fill_max"),
  Value = c(page_width_pt, page_height_pt, cell_pt, cell_pt, cell_pt, font_pt, 0, fill_max)
), file.path(work, "figure_geometry.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

cat(sprintf("Figure: %.2f x %.2f pt (%.2f x %.2f in)\n", page_width_pt, page_height_pt, page_width_pt/72, page_height_pt/72))
cat(sprintf("Cells: %.2f x %.2f pt; font: %.1f pt; fill: 0-%.1f\n", cell_pt, cell_pt, font_pt, fill_max))
