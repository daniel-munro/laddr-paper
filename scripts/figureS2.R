## Effect on xQTLs of including more GTEx tissues and TCGA data for latent models

library(tidyverse)
library(patchwork)

qtls_gtex5_full <- read_tsv("data/processed/gtex5-full.qtls.tsv.gz", col_types = "cciccd")

qtls_gtex_full <- read_tsv("data/processed/gtex-full.qtls.tsv.gz", col_types = "cciccd")

qtls_gtextcga_full <- read_tsv("data/processed/gtextcga-full.qtls.tsv.gz", col_types = "cciccd")

tissues5 <- read_lines("data/info/tissues.gtex5.txt")

gtex_colors <- read_tsv(
  "data/pantry/gtex/tissueInfo.tsv",
  col_types = cols(tissueSiteDetailAbbr = "c", colorHex = "c", .default = "-")
) |>
  mutate(colorHex = str_c("#", colorHex)) |>
  deframe()

qtl_counts <- full_join(
  qtls_gtex5_full |>
    count(tissue, name = "n_gtex5"),
  qtls_gtex_full |>
    count(tissue, name = "n_gtex"),
  by = "tissue",
  relationship = "one-to-one"
) |>
  full_join(
    qtls_gtextcga_full |>
      count(tissue, name = "n_gtextcga"),
    by = "tissue",
    relationship = "one-to-one"
  )

xy_limit <- with(qtl_counts, max(n_gtex5, n_gtex, n_gtextcga)) / 1000 * 1.04

p1 <- qtl_counts |>
  mutate(in_5_tissues = if_else(tissue %in% tissues5, "True", "False")) |>
  ggplot(aes(x = n_gtex5 / 1000, y = n_gtex / 1000, fill = tissue, shape = in_5_tissues)) +
  geom_abline(slope = 1, intercept = 0, color = "black", linewidth = 0.3) +
  annotate("text", label = "y = x", x = 8, y = 7, hjust = 0) +
  geom_point(alpha = 0.75) +
  # scale_color_manual(values = c("black", "red")) +
  scale_fill_manual(values = gtex_colors, guide = "none") +
  scale_shape_manual(values = c(21, 24)) +
  expand_limits(x = c(0, xy_limit), y = c(0, xy_limit)) +
  coord_fixed(expand = 0) +
  theme_bw() +
  theme(
    axis.text = element_text(color = "black"),
    legend.position = "inside",
    legend.position.inside = c(0.75, 0.2),
    panel.grid = element_blank(),
  ) +
  xlab("xQTLs (×1000), 5-tissue models") +
  ylab("xQTLs (×1000), GTEx models") +
  labs(shape = "Tissue used for\n5-tissue models")

p2 <- qtl_counts |>
  ggplot(aes(x = n_gtex / 1000, y = n_gtextcga / 1000, fill = tissue)) +
  geom_abline(slope = 1, intercept = 0, color = "black", linewidth = 0.3) +
  annotate("text", label = "y = x", x = 8, y = 7, hjust = 0) +
  geom_point(shape = 21, color = "black", alpha = 0.75, show.legend = FALSE) +
  scale_fill_manual(values = gtex_colors) +
  expand_limits(x = c(0, xy_limit), y = c(0, xy_limit)) +
  coord_fixed(expand = 0) +
  theme_bw() +
  theme(
    axis.text = element_text(color = "black"),
    legend.position = "inside",
    legend.position.inside = c(0.9, 0.1),
    panel.grid = element_blank(),
  ) +
  xlab("xQTLs (×1000), GTEx models") +
  ylab("xQTLs (×1000), GTEx + TCGA models")

p1 + p2 + plot_annotation(tag_levels = "a") & theme(plot.tag = element_text(face = "bold"))

ggsave("figures/figureS2.png", width = 7.5, height = 3.75, device = png)
