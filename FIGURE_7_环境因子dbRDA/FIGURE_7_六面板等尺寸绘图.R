rm(list = ls())
setwd('/Users/haleoriu/Documents/Takakia绘图和文章修改/outputs/Figure7_environment_association_20261013/evidence_priority_no_AN_v1')
library(ggplot2)
library(ggrepel)
library(cowplot)
library(grid)
set.seed(20261013)

# 使用已经拟合并核验的排序坐标和模型统计量。
preset <- 'Elevation_evidencePriorityVIF10'
sites <- read.delim('tables/dbRDA_sites.tsv', check.names = FALSE)
arrows <- read.delim('tables/dbRDA_arrows.tsv', check.names = FALSE)
statistics <- read.delim('tables/dbRDA_statistics.tsv', check.names = FALSE)
profile_order <- c('R bacteria', 'E bacteria', 'R fungi', 'E fungi', 'R KO', 'E KO')

plots <- lapply(profile_order, function(profile) {
  ss <- sites[sites$Preset == preset & sites$Profile == profile, ]
  aa <- arrows[arrows$Preset == preset & arrows$Profile == profile, ]
  st <- statistics[statistics$Preset == preset & statistics$Profile == profile, ]
  stopifnot(nrow(ss) == 30, nrow(aa) == 5, nrow(st) == 1)
  ggplot(ss, aes(x, y)) +
    geom_hline(yintercept = 0, color = 'grey85', linewidth = .25) +
    geom_vline(xintercept = 0, color = 'grey85', linewidth = .25) +
    geom_segment(data = aa, aes(x = 0, y = 0, xend = x, yend = y),
                 color = 'grey45', linewidth = .45,
                 arrow = grid::arrow(length = grid::unit(1.8, 'mm'))) +
    geom_point(shape = 21,
               fill = if (substr(profile, 1, 1) == 'R') '#82B3D7' else '#CF6C72',
               color = 'grey25', size = 2.3, stroke = .3) +
    geom_text_repel(data = aa, aes(label = Label), size = 3,
                    family = 'Helvetica', seed = 20261013,
                    max.overlaps = Inf, min.segment.length = 0) +
    labs(title = profile,
         subtitle = sprintf('adj. R² = %.3f; p = %.3f', st$Adjusted_R2, st$p),
         x = sprintf('db-RDA1 (%.1f%%)', st$Axis1pct),
         y = sprintf('db-RDA2 (%.1f%%)', st$Axis2pct)) +
    theme_bw(base_size = 11, base_family = 'Helvetica') +
    theme(panel.grid = element_blank(),
          plot.title = element_text(face = 'bold', hjust = .5),
          plot.subtitle = element_text(size = 9, hjust = .5),
          plot.margin = margin(5.5, 5.5, 5.5, 5.5))
})

# 全部六图共用相同的外部排版尺寸，使绘图区边框宽高一致。
# 坐标轴刻度文字宽度不同，也不会改变单个面板的绘图区大小。
grobs <- lapply(plots, ggplotGrob)
stopifnot(length(unique(vapply(grobs, function(g) length(g$widths), integer(1)))) == 1,
          length(unique(vapply(grobs, function(g) length(g$heights), integer(1)))) == 1)
common_widths <- do.call(grid::unit.pmax, lapply(grobs, function(g) g$widths))
common_heights <- do.call(grid::unit.pmax, lapply(grobs, function(g) g$heights))
grobs <- lapply(grobs, function(g) {
  g$widths <- common_widths
  g$heights <- common_heights
  g
})
figure <- plot_grid(plotlist = grobs, ncol = 2, nrow = 3,
                    rel_widths = c(1, 1), rel_heights = c(1, 1, 1),
                    labels = letters[1:6], label_size = 15)
ggsave('figures/Figure7_dbRDA_Elevation_equal_panels.pdf', figure,
       width = 9, height = 11.8, device = 'pdf', family = 'Helvetica',
       useDingbats = FALSE)
ggsave('figures/Figure7_dbRDA_Elevation_equal_panels.png', figure,
       width = 9, height = 11.8, dpi = 300)
