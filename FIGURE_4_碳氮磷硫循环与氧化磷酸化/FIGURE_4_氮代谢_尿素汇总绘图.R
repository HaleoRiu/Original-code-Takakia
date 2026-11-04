options(stringsAsFactors = FALSE, warn = 1)
suppressPackageStartupMessages(library(grid))

script_arg <- grep("^--file=", commandArgs(), value = TRUE)
script_dir <- if (length(script_arg)) dirname(normalizePath(sub("^--file=", "", script_arg[1]))) else getwd()
project_dir <- normalizePath(file.path(script_dir, ".."), mustWork = FALSE)
work_root <- file.path(project_dir, "work")

plot_element <- function(code, element_name, split_at, label_width,
                         panel_titles, panel_subtitles) {
  work <- file.path(work_root, code)
  catalog <- read.delim(file.path(work, paste0(code, "_gene_catalog.tsv")), check.names = FALSE)
  heat <- read.delim(file.path(work, paste0(code, "_gene_heat_values.tsv")), check.names = FALSE)
  alt_tests <- read.delim(file.path(work, paste0(code, "_gene_Wilcoxon_by_altitude.tsv")), check.names = FALSE)
  total_tests <- read.delim(file.path(work, paste0(code, "_gene_Wilcoxon_Total.tsv")), check.names = FALSE)
  kw <- read.delim(file.path(work, paste0(code, "_gene_KW_altitude.tsv")), check.names = FALSE)
  sum_catalog <- read.delim(file.path(work, paste0(code, "_sum_catalog.tsv")), check.names = FALSE)
  sum_heat <- read.delim(file.path(work, paste0(code, "_sum_heat_values.tsv")), check.names = FALSE)
  sum_alt_tests <- read.delim(file.path(work, paste0(code, "_sum_Wilcoxon_by_altitude.tsv")), check.names = FALSE)
  sum_total_tests <- read.delim(file.path(work, paste0(code, "_sum_Wilcoxon_Total.tsv")), check.names = FALSE)
  sum_kw <- read.delim(file.path(work, paste0(code, "_sum_KW_altitude.tsv")), check.names = FALSE)
  scale_info <- read.delim(file.path(work, "shared_fill_scale.tsv"), check.names = FALSE)

  catalog <- catalog[order(catalog$Feature_order), , drop = FALSE]
  sum_catalog <- sum_catalog[order(sum_catalog$Feature_order), , drop = FALSE]
  left <- catalog[catalog$Feature_order <= split_at, , drop = FALSE]
  right <- catalog[catalog$Feature_order > split_at, , drop = FALSE]
  sum_left_n <- ceiling(nrow(sum_catalog) / 2)
  sum_left <- sum_catalog[seq_len(sum_left_n), , drop = FALSE]
  sum_right <- sum_catalog[-seq_len(sum_left_n), , drop = FALSE]

  cell_pt <- 32.33
  font_pt <- 28
  line_col <- "#4D4D4D"
  alts <- c(3800, 3900, 4000, 4100, 4200)
  habitats <- c("R", "E")
  fill_max <- scale_info$Fixed_maximum[1]
  palette <- colorRampPalette(c("#2C7BB6", "#F7F7F7", "#D7191C"))(512)

  fill_for <- function(z) {
    if (!is.finite(z)) z <- 0
    idx <- pmax(1, pmin(length(palette), floor((z / fill_max) * (length(palette) - 1)) + 1))
    palette[idx]
  }

  kw_col_pt <- 48
  chart_width <- 10 * cell_pt + 16 + 2 * kw_col_pt + 10 + 2 * cell_pt
  col_gap <- 64
  left_margin <- 38
  right_margin <- 38
  top_margin <- 190
  bottom_margin <- 70
  panel_block_rows <- function(n) 1.45 + 2.55 + n + 1.35
  gene_block_rows <- max(panel_block_rows(nrow(left)), panel_block_rows(nrow(right)))
  sum_y0 <- top_margin + gene_block_rows * cell_pt + 34
  sum_block_rows <- max(panel_block_rows(nrow(sum_left)), panel_block_rows(nrow(sum_right)))
  page_height_pt <- sum_y0 + sum_block_rows * cell_pt + bottom_margin
  page_width_pt <- left_margin + 2 * (label_width + chart_width) + col_gap + right_margin
  col_x <- c(left = left_margin, right = left_margin + label_width + chart_width + col_gap)

  draw_text <- function(label, x, y, just = "centre", face = "plain",
                        size = font_pt, rot = 0, col = "black") {
    label <- gsub("[\u2013\u2014]", "-", label)
    grid.text(label, x = unit(x, "native"), y = unit(y, "native"),
              just = just, rot = rot,
              gp = gpar(fontfamily = "Helvetica", fontface = face,
                        fontsize = size, col = col))
  }

  draw_header <- function(x0, y0, first_label) {
    hx <- x0 + label_width
    draw_text(first_label, hx - 14, y0 + 1.25 * cell_pt, just = "right", face = "bold")
    for (i in seq_along(alts)) {
      cx <- hx + ((i - 1) * 2 + 1) * cell_pt
      draw_text(as.character(alts[i]), cx, y0 + 0.75 * cell_pt,
                face = "bold", rot = 45)
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

  heat_record <- function(table, fid, scope, altitude, habitat) {
    hit <- table$Feature_ID == fid & table$Scope == scope & table$Habitat == habitat
    if (scope == "Altitude") hit <- hit & table$Altitude == altitude
    out <- table[hit, , drop = FALSE]
    if (nrow(out) != 1L) stop("Expected one heat row for ", fid, " / ", scope, " / ", altitude, " / ", habitat)
    out[1, ]
  }

  draw_cell <- function(cx, cy, rec, sig = "", shape = "square") {
    gp <- gpar(fill = fill_for(rec$log2_mean_TPM_plus1), col = "white", lwd = 0.55)
    if (shape == "circle") {
      grid.circle(x = unit(cx, "native"), y = unit(cy, "native"),
                  r = unit(cell_pt / 2, "pt"), gp = gp)
    } else {
      grid.rect(x = unit(cx, "native"), y = unit(cy, "native"),
                width = unit(cell_pt, "pt"), height = unit(cell_pt, "pt"), gp = gp)
    }
    if (isTRUE(all.equal(rec$Mean_TPM, 0))) draw_text("X", cx, cy, face = "bold")
    if (nzchar(sig)) draw_text(sig, cx, cy)
  }

  draw_row <- function(rec, y, x0, tables, shape = "square", italic = TRUE) {
    fid <- rec$Feature_ID
    hx <- x0 + label_width
    draw_text(rec$Display_label, hx - 14, y, just = "right",
              face = if (italic) "italic" else "plain")
    aa <- tables$alt[tables$alt$Feature_ID == fid, , drop = FALSE]
    tt <- tables$total[tables$total$Feature_ID == fid, , drop = FALSE]
    kk <- tables$kw[tables$kw$Feature_ID == fid, , drop = FALSE]

    for (i in seq_along(alts)) {
      arow <- aa[aa$Altitude == alts[i], , drop = FALSE]
      for (j in seq_along(habitats)) {
        hab <- habitats[j]
        cx <- hx + ((i - 1) * 2 + j - 0.5) * cell_pt
        hrec <- heat_record(tables$heat, fid, "Altitude", alts[i], hab)
        sig <- if (nrow(arow) == 1L && arow$Higher[1] == hab) arow$Significance[1] else ""
        draw_cell(cx, y, hrec, sig, shape)
      }
    }

    kw0 <- hx + 10 * cell_pt + 16
    for (j in seq_along(habitats)) {
      sig <- kk$Significance[kk$Habitat == habitats[j]][1]
      draw_text(sig, kw0 + (j - 0.5) * kw_col_pt, y)
    }

    total0 <- kw0 + 2 * kw_col_pt + 10
    for (j in seq_along(habitats)) {
      hab <- habitats[j]
      cx <- total0 + (j - 0.5) * cell_pt
      hrec <- heat_record(tables$heat, fid, "Total", NA_integer_, hab)
      sig <- if (nrow(tt) == 1L && tt$Higher[1] == hab) tt$Significance[1] else ""
      draw_cell(cx, y, hrec, sig, shape)
    }
  }

  draw_separators <- function(d, x0, y_start) {
    if (nrow(d) < 2L) return(invisible(NULL))
    breaks <- which(d$Section[-1] != d$Section[-nrow(d)])
    for (b in breaks) {
      yy <- y_start + b * cell_pt
      grid.lines(
        x = unit(c(x0 + 4, x0 + label_width + chart_width), "native"),
        y = unit(c(yy, yy), "native"),
        gp = gpar(col = line_col, lwd = 1.05)
      )
    }
  }

  draw_panel <- function(panel_title, subtitle, d, x0, y0, tables,
                         shape = "square", italic = TRUE, first_label = "Gene") {
    draw_text(panel_title, x0, y0 + 0.58 * cell_pt, just = "left", face = "bold")
    draw_text(subtitle, x0 + label_width + chart_width, y0 + 0.58 * cell_pt,
              just = "right")
    header_y <- y0 + 1.45 * cell_pt
    draw_header(x0, header_y, first_label)
    row_start <- header_y + 2.55 * cell_pt
    if (nrow(d)) {
      for (i in seq_len(nrow(d))) {
        draw_row(d[i, ], row_start + (i - 0.5) * cell_pt, x0, tables, shape, italic)
      }
      if (shape == "square") draw_separators(d, x0, row_start)
    }
  }

  draw_legend <- function() {
    draw_text(paste0(element_name, "-cycle genes and requested sums across altitude"),
              left_margin, 46, just = "left", face = "bold")
    draw_text("Squares: individual KO/gene rows; circles: within-sample TPM sums",
              left_margin, 88, just = "left")
    lx <- page_width_pt - 590
    ly <- 78
    legend_w <- 420
    n <- 140
    for (i in seq_len(n)) {
      x <- lx + (i - 0.5) * legend_w / n
      z <- (i - 1) / (n - 1) * fill_max
      grid.rect(x = unit(x, "native"), y = unit(ly, "native"),
                width = unit(legend_w / n + 0.2, "pt"), height = unit(28, "pt"),
                gp = gpar(fill = fill_for(z), col = NA))
    }
    draw_text("log2(mean TPM + 1)", lx + legend_w / 2, ly - 35, face = "bold")
    for (tv in seq(0, fill_max, by = 3)) {
      tx <- lx + tv / fill_max * legend_w
      grid.lines(x = unit(c(tx, tx), "native"),
                 y = unit(c(ly + 14, ly + 23), "native"),
                 gp = gpar(col = line_col, lwd = 1.0))
      draw_text(as.character(tv), tx, ly + 48)
    }
  }

  gene_tables <- list(heat = heat, alt = alt_tests, total = total_tests, kw = kw)
  sum_tables <- list(heat = sum_heat, alt = sum_alt_tests, total = sum_total_tests, kw = sum_kw)

  draw_page <- function() {
    grid.newpage()
    pushViewport(viewport(xscale = c(0, page_width_pt), yscale = c(page_height_pt, 0)))
    grid.rect(gp = gpar(fill = "white", col = NA))
    draw_legend()
    draw_panel(panel_titles[1], panel_subtitles[1], left, col_x["left"], top_margin,
               gene_tables, "square", TRUE, "Gene")
    draw_panel(panel_titles[2], panel_subtitles[2], right, col_x["right"], top_margin,
               gene_tables, "square", TRUE, "Gene")
    draw_panel(paste0(element_name, " sums I"), "Requested functional sums",
               sum_left, col_x["left"], sum_y0, sum_tables, "circle", FALSE, "Sum")
    draw_panel(paste0(element_name, " sums II"), "Requested functional sums",
               sum_right, col_x["right"], sum_y0, sum_tables, "circle", FALSE, "Sum")
    popViewport()
  }

  pdf_path <- file.path(work, paste0(code, "_genes_and_requested_sums_altitude_28pt.pdf"))
  pdf(pdf_path, width = page_width_pt / 72, height = page_height_pt / 72,
      family = "Helvetica", useDingbats = FALSE, onefile = TRUE)
  draw_page()
  dev.off()

  write.table(data.frame(
    Item = c("PDF_width_pt", "PDF_height_pt", "Square_width_pt", "Square_height_pt",
             "Circle_diameter_pt", "Font_size_pt", "Shared_fill_min", "Shared_fill_max",
             "Displayed_gene_rows", "Displayed_sum_rows"),
    Value = c(page_width_pt, page_height_pt, cell_pt, cell_pt, cell_pt, font_pt,
              0, fill_max, nrow(catalog), nrow(sum_catalog))
  ), file.path(work, "figure_geometry.tsv"), sep = "\t",
  quote = FALSE, row.names = FALSE)

  cat(sprintf("%s PDF: %s\n", code, pdf_path))
}

plot_element(
  "N", "Nitrogen", 29, 610,
  c("Nitrogen I", "Nitrogen II"),
  c("Fixation, oxidation and nitrate/nitrite reduction",
    "N uptake, urea hydrolysis/maturation and ammonium assimilation")
)
