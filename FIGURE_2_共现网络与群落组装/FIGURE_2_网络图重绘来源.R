#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(igraph)
  library(ggplot2)
  library(dplyr)
  library(cowplot)
})

args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
script_path <- if (length(file_arg)) {
  normalizePath(sub("^--file=", "", file_arg[1]))
} else {
  normalizePath("redraw_networks.R")
}
network_dir <- dirname(script_path)

edge_colours <- c(
  "Positive" = "#F17C9A",
  "Negative" = "#57B894"
)

habitat_colours <- c(
  "Rhizomatousphere" = "#82B3D7",
  "Endophytic" = "#CF6C72"
)

module_colours <- c(
  "Module 1" = "#66A5AD",
  "Module 2" = "#E9BE63",
  "Module 3" = "#D98B86",
  "Module 4" = "#8FBF88",
  "Module 5" = "#B49AD6",
  "Module 6" = "#E7A3BF",
  "Others" = "#BDBDBD"
)

network_specs <- tibble::tribble(
  ~id, ~title, ~kind, ~file, ~node_file, ~edge_file,
  "bacteria_across", "Bacteria across-habitats", "across", NA_character_, "bacteria_node_n.csv", "bacteria_edge_n.csv",
  "fungi_across", "Fungi across-habitats", "across", NA_character_, "fungi_node.csv", "fungi_edge.csv",
  "bacteria_R", "Bacteria Rhizomatousphere", "module", "bac_rhizo.graphml", NA_character_, NA_character_,
  "fungi_R", "Fungi Rhizomatousphere", "module", "fungi_rhizo.graphml", NA_character_, NA_character_,
  "bacteria_E", "Bacteria Endophytic", "module", "bac_endo.graphml", NA_character_, NA_character_,
  "fungi_E", "Fungi Endophytic", "module", "fungi_endo.graphml", NA_character_, NA_character_
)

read_across_graph <- function(node_file, edge_file) {
  nodes <- read.csv(
    file.path(network_dir, node_file),
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  edges <- read.csv(
    file.path(network_dir, edge_file),
    check.names = FALSE,
    stringsAsFactors = FALSE
  )

  colnames(nodes)[colnames(nodes) == "ID"] <- "name"
  used_nodes <- unique(c(edges$source, edges$target))
  nodes <- nodes %>%
    filter(name %in% used_nodes) %>%
    mutate(
      Habitat = case_when(
        group == "fun" ~ "Rhizomatousphere",
        group == "bac" ~ "Endophytic",
        TRUE ~ "Others"
      ),
      NodeGroup = Habitat
    )

  edges <- edges %>%
    transmute(
      from = source,
      to = target,
      corr = as.numeric(correlation),
      EdgeSign = if_else(cor == "+", "Positive", "Negative"),
      weight = abs(corr)
    )

  graph_from_data_frame(edges, directed = FALSE, vertices = nodes)
}

read_module_graph <- function(graphml_file) {
  g <- read_graph(file.path(network_dir, graphml_file), format = "graphml")
  g <- as_undirected(g, mode = "collapse")
  E(g)$corr_num <- suppressWarnings(as.numeric(E(g)$corr))
  E(g)$corr_num[is.na(E(g)$corr_num)] <- suppressWarnings(as.numeric(E(g)$weight[is.na(E(g)$corr_num)]))
  E(g)$corr_num[is.na(E(g)$corr_num)] <- 1
  E(g)$EdgeSign <- ifelse(
    grepl("neg", E(g)$cor_type, ignore.case = TRUE) | E(g)$corr_num < 0,
    "Negative",
    "Positive"
  )
  E(g)$weight <- abs(E(g)$corr_num)

  community <- cluster_louvain(g, weights = E(g)$weight)
  module_id <- membership(community)
  module_sizes <- sort(table(module_id), decreasing = TRUE)
  top_modules <- names(module_sizes)[seq_len(min(6, length(module_sizes)))]
  module_labels <- setNames(
    paste0("Module ", seq_along(top_modules)),
    top_modules
  )
  V(g)$NodeGroup <- ifelse(
    as.character(module_id) %in% top_modules,
    module_labels[as.character(module_id)],
    "Others"
  )
  g
}

make_graph <- function(spec) {
  if (spec$kind == "across") {
    read_across_graph(spec$node_file, spec$edge_file)
  } else {
    read_module_graph(spec$file)
  }
}

normalize_layout <- function(layout_matrix) {
  x <- layout_matrix[, 1]
  y <- layout_matrix[, 2]
  x <- x - mean(range(x))
  y <- y - mean(range(y))
  scale_factor <- max(diff(range(x)), diff(range(y)))
  if (scale_factor == 0) scale_factor <- 1
  cbind(x / scale_factor, y / scale_factor)
}

circularize_layout <- function(layout_matrix, strength = 0.42) {
  layout_matrix <- normalize_layout(layout_matrix)
  x <- layout_matrix[, 1]
  y <- layout_matrix[, 2]
  r <- sqrt(x^2 + y^2)
  theta <- atan2(y, x)
  zero_r <- r == 0 | is.na(r)
  if (any(zero_r)) {
    theta[zero_r] <- seq(0, 2 * pi, length.out = sum(zero_r) + 1)[-1]
    r[zero_r] <- 0.01
  }
  r_rank <- rank(r, ties.method = "average") / length(r)
  target_r <- 0.05 + 0.47 * sqrt(r_rank)
  circular_x <- target_r * cos(theta)
  circular_y <- target_r * sin(theta)
  blended <- cbind(
    (1 - strength) * x + strength * circular_x,
    (1 - strength) * y + strength * circular_y
  )
  normalize_layout(blended)
}


prepare_plot_data <- function(g, title, seed) {
  set.seed(seed)
  weights <- E(g)$weight
  weights[is.na(weights) | weights <= 0] <- 1
  initial_layout <- layout_in_circle(g)
  layout_matrix <- layout_with_fr(
    g,
    coords = initial_layout,
    weights = weights,
    niter = 4000,
    start.temp = sqrt(vcount(g)) * 1.7,
    grid = "nogrid"
  )
  layout_matrix <- circularize_layout(layout_matrix, strength = 0.58)

  nodes <- igraph::as_data_frame(g, what = "vertices") %>%
    mutate(
      name = if ("name" %in% colnames(.)) name else as.character(seq_len(n())),
      x = layout_matrix[, 1],
      y = layout_matrix[, 2],
      degree_value = degree(g),
      degree_sqrt = sqrt(degree_value),
      NodeGroup = if_else(is.na(NodeGroup), "Others", NodeGroup),
      title = title
    )

  edges <- igraph::as_data_frame(g, what = "edges") %>%
    left_join(nodes %>% select(name, x, y), by = c("from" = "name")) %>%
    rename(x_from = x, y_from = y) %>%
    left_join(nodes %>% select(name, x, y), by = c("to" = "name")) %>%
    rename(x_to = x, y_to = y) %>%
    mutate(
      EdgeSign = if_else(is.na(EdgeSign), "Positive", EdgeSign),
      title = title
    )

  list(nodes = nodes, edges = edges)
}

rescale_size <- function(x, local_min, local_max, out_min = 0.85, out_max = 4.85) {
  if (local_max == local_min) {
    rep(mean(c(out_min, out_max)), length(x))
  } else {
    out_min + (x - local_min) / (local_max - local_min) * (out_max - out_min)
  }
}

graphs <- lapply(seq_len(nrow(network_specs)), function(i) make_graph(network_specs[i, ]))
plot_data <- lapply(seq_along(graphs), function(i) {
  prepare_plot_data(graphs[[i]], network_specs$title[i], seed = 20260630 + i)
})

all_degree_sqrt <- unlist(lapply(plot_data, function(x) x$nodes$degree_sqrt))
global_min_degree <- min(all_degree_sqrt, na.rm = TRUE)
global_max_degree <- max(all_degree_sqrt, na.rm = TRUE)

plot_data <- lapply(plot_data, function(x) {
  local_min_degree <- min(x$nodes$degree_sqrt, na.rm = TRUE)
  local_max_degree <- max(x$nodes$degree_sqrt, na.rm = TRUE)
  x$nodes <- x$nodes %>%
    mutate(
      NodeSize = rescale_size(degree_sqrt, local_min_degree, local_max_degree),
      local_degree_sqrt_min = local_min_degree,
      local_degree_sqrt_max = local_max_degree,
      global_degree_sqrt_min = global_min_degree,
      global_degree_sqrt_max = global_max_degree
    )
  x
})

plot_network <- function(x, title, kind, label = NULL) {
  node_palette <- if (kind == "across") habitat_colours else module_colours

  ggplot() +
    geom_segment(
      data = x$edges,
      aes(x = x_from, y = y_from, xend = x_to, yend = y_to, colour = EdgeSign),
      linewidth = 0.18,
      alpha = 0.58,
      lineend = "round"
    ) +
    geom_point(
      data = x$nodes,
      aes(x = x, y = y, fill = NodeGroup, size = NodeSize),
      shape = 21,
      colour = "#2F2F2F",
      stroke = 0.14,
      alpha = 0.96
    ) +
    scale_colour_manual(values = edge_colours, drop = FALSE) +
    scale_fill_manual(values = node_palette, drop = FALSE) +
    scale_size_identity() +
    coord_equal(clip = "off") +
    labs(title = title) +
    theme_void(base_family = "Helvetica") +
    theme(
      plot.title = element_text(face = "bold", size = 13, hjust = 0.5),
      plot.margin = margin(4, 4, 4, 4),
      legend.position = "none"
    ) +
    {
      if (!is.null(label)) {
        annotate(
          "text",
          x = -0.56,
          y = 0.56,
          label = label,
          hjust = 0,
          vjust = 1,
          family = "Helvetica",
          fontface = "bold",
          size = 7
        )
      }
    }
}

panel_labels <- paste0("(", letters[seq_len(nrow(network_specs))], ")")
network_plots <- lapply(seq_len(nrow(network_specs)), function(i) {
  plot_network(
    plot_data[[i]],
    network_specs$title[i],
    network_specs$kind[i],
    panel_labels[i]
  )
})

names(network_plots) <- network_specs$id

combined_plot <- plot_grid(
  network_plots[["bacteria_across"]],
  network_plots[["fungi_across"]],
  network_plots[["bacteria_R"]],
  network_plots[["fungi_R"]],
  network_plots[["bacteria_E"]],
  network_plots[["fungi_E"]],
  ncol = 2,
  align = "hv"
)

save_plot_pair <- function(filename_stem, plot, width, height) {
  ggsave(
    file.path(network_dir, paste0(filename_stem, ".pdf")),
    plot,
    width = width,
    height = height,
    units = "in",
    device = "pdf",
    useDingbats = FALSE,
    bg = "white"
  )
  ggsave(
    file.path(network_dir, paste0(filename_stem, ".png")),
    plot,
    width = width,
    height = height,
    units = "in",
    dpi = 600,
    bg = "white"
  )
}

for (i in seq_len(nrow(network_specs))) {
  save_plot_pair(
    paste0(network_specs$id[i], "_network_redraw"),
    network_plots[[i]],
    width = 4.0,
    height = 4.0
  )
}

save_plot_pair(
  "six_networks_redraw_combined",
  combined_plot,
  width = 8.2,
  height = 11.4
)

summary_table <- bind_rows(lapply(seq_len(nrow(network_specs)), function(i) {
  tibble(
    id = network_specs$id[i],
    title = network_specs$title[i],
    nodes = vcount(graphs[[i]]),
    edges = ecount(graphs[[i]]),
    raw_degree_min = min(plot_data[[i]]$nodes$degree_value),
    raw_degree_max = max(plot_data[[i]]$nodes$degree_value),
    global_raw_degree_min = min(unlist(lapply(plot_data, function(x) x$nodes$degree_value))),
    global_raw_degree_max = max(unlist(lapply(plot_data, function(x) x$nodes$degree_value))),
    size_rule = "panel degree min/max mapped to shared visual size range 0.85-4.85",
    node_size_min = min(plot_data[[i]]$nodes$NodeSize),
    node_size_max = max(plot_data[[i]]$nodes$NodeSize)
  )
}))

write.csv(
  summary_table,
  file.path(network_dir, "network_redraw_summary.csv"),
  row.names = FALSE
)

message("Created six redrawn network panels and combined layout.")
