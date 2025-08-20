## Factors influencing performance

library(tidyverse)

####################
## Panels a and b ## Effect on xQTLs of including more GTEx tissues and TCGA data for latent models
####################

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

qtl_counts |>
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

ggsave("figures/figure5/figure5a.png", width = 3.75, height = 3.75, device = png)

qtl_counts |>
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

ggsave("figures/figure5/figure5b.png", width = 3.75, height = 3.75, device = png)

#############
## Panel c ## Pruned annotation xQTLs
#############

modalities <- c(
  latent = "Latent",
  alt_TSS = "Alt. TSS",
  alt_polyA = "Alt. polyA",
  isoforms = "Isoform ratio",
  stability = "RNA stability",
  splicing = "Intron excision",
  expression = "Expression"
)

modality_colors <- c(
  `Latent` = "#13918d",
  `Alt. TSS` = "#896090",
  `Alt. polyA` = "#d97f26",
  `Isoform ratio` = "#6a90cd",
  `RNA stability` = "#ddb23c",
  `Intron excision` = "#59a257",
  Expression = "#bf4042"
)

map_groups <- c(
  latent = "Data-driven",
  pantry = "Knowledge-driven"
)

qtls_prune <- read_tsv("data/processed/prune-BRNCTXB.qtls.tsv.gz", col_types = "cicciccd") |>
  filter(map_group %in% names(map_groups)) |>
  mutate(modality = factor(modalities[modality], levels = names(modality_colors)),
         map_group = factor(map_groups[map_group], levels = map_groups),
         pruning = fct_reorder(as.character(pruning), pruning))

ylims <- qtls_prune |>
  count(map_group, pruning) |>
  mutate(ylim = max(n) * 1.03 / 1000, .by = map_group)

qtls_prune |>
  count(map_group, pruning, modality) |>
  filter(map_group %in% c("Data-driven", "Knowledge-driven")) |>
  ggplot(aes(x = pruning, y = n / 1000, fill = modality)) +
  facet_wrap(~map_group, scales = "free_y") +
  geom_col(width = 0.8) +
  geom_point(aes(y = ylim, fill = NULL), data = ylims, color = "white", show.legend = FALSE) +
  scale_y_continuous(breaks = c(0, 5, 10, 15, 20), expand = c(0, 0)) +
  scale_fill_manual(values = modality_colors) +
  theme_bw() +
  theme(
    axis.text = element_text(color = "black"),
    legend.key.size = unit(12, "pt"),
    panel.grid = element_blank(),
  ) +
  xlab("% of non-canonical isoforms kept") +
  ylab("xQTLs (×1000)") +
  labs(fill = "Modality")

ggsave("figures/figure5/figure5c.png", width = 5, height = 3, device = png)
