## Additional colocalization to xQTL ratios

library(tidyverse)
library(patchwork)

modalities <- c(
  expression = "Expression",
  isoforms = "Isoform ratio",
  splicing = "Intron excision",
  alt_TSS = "Alt. TSS",
  alt_polyA = "Alt. polyA",
  stability = "RNA stability",
  latent_residual = "Residual DD"
)

modality_colors <- c(
  Expression = "#bf4042",
  `Isoform ratio` = "#6a90cd",
  `Intron excision` = "#59a257",
  `Alt. TSS` = "#896090",
  `Alt. polyA` = "#d97f26",
  `RNA stability` = "#ddb23c",
  `Residual DD` = "#1ce6df"
)

qtls_hybrid <- read_tsv(
  "data/processed/gtex-residual-cross.qtls.tsv.gz",
  col_types = "ccicccd"
) |>
  mutate(modality = factor(modalities[modality], levels = rev(modalities)))

twas <- bind_rows(
  read_tsv("data/processed/gtex-residual.twas_hits.tsv.gz", col_types = "cccc--dd") |>
    mutate(modality = "latent_residual", .after = "tissue"),
  read_tsv("data/processed/gtex-pantry.twas_hits.tsv.gz", col_types = "ccccc--dd")
)

####################
## Panels a, b, c ## Coloc TWAS genes to xGenes ratio
####################

coloc_genes <- twas |>
  filter(coloc_pp > 0.8) |>
  distinct(tissue, modality, gene_id) |>
  count(tissue, modality, name = "coloc_n_genes") |>
  mutate(modality = factor(modalities[modality], levels = rev(modalities)))

coloc_genes_xgene_ratio <- qtls_hybrid |>
  distinct(tissue, modality, gene_id) |>
  count(tissue, modality, name = "n_xgenes") |>
  full_join(coloc_genes, by = c("tissue", "modality"), relationship = "one-to-one") |>
  mutate(coloc_xgene_ratio = coloc_n_genes / n_xgenes)

p1 <- coloc_genes_xgene_ratio |>
  ggplot(aes(x = coloc_xgene_ratio, y = modality, fill = modality)) +
  geom_boxplot(show.legend = FALSE, outlier.size = 0.3, linewidth = 0.4, median.linewidth = 0.5) +
  scale_fill_manual(values = modality_colors) +
  theme_classic() +
  xlab("Genes w. colocalizing TWAS hits / xGenes") +
  ylab("Modality")

p2 <- coloc_genes_xgene_ratio |>
  ggplot(aes(x = coloc_n_genes / 1000, y = modality, fill = modality)) +
  geom_boxplot(show.legend = FALSE, outlier.size = 0.3, linewidth = 0.4, median.linewidth = 0.5) +
  scale_fill_manual(values = modality_colors) +
  theme_classic() +
  xlab("Genes w. colocalizing TWAS hits (×1000)") +
  ylab(NULL)

p3 <- coloc_genes_xgene_ratio |>
  ggplot(aes(x = n_xgenes / 1000, y = modality, fill = modality)) +
  geom_boxplot(show.legend = FALSE, outlier.size = 0.3, linewidth = 0.4, median.linewidth = 0.5) +
  scale_fill_manual(values = modality_colors) +
  theme_classic() +
  xlab("xGenes (×1000)") +
  ylab(NULL)

####################
## Panels d, e, f ## Top coloc TWAS to top xQTL ratio
####################

coloc_top <- twas |>
  filter(coloc_pp > 0.8) |>
  slice_min(n = 1, order_by = twas_p, by = c("tissue", "trait", "gene_id")) |>
  count(tissue, modality, name = "coloc_top_n") |>
  mutate(modality = factor(modalities[modality], levels = rev(modalities)))

coloc_top_qtl_ratio <- qtls_hybrid |>
  filter(rank == 1) |>
  count(tissue, modality, name = "n_top_qtls") |>
  full_join(coloc_top, by = c("tissue", "modality"), relationship = "one-to-one") |>
  mutate(coloc_top_qtl_ratio = coloc_top_n / n_top_qtls)

p4 <- coloc_top_qtl_ratio |>
  ggplot(aes(x = coloc_top_qtl_ratio, y = modality, fill = modality)) +
  geom_boxplot(show.legend = FALSE, outlier.size = 0.3, linewidth = 0.4, median.linewidth = 0.5) +
  scale_fill_manual(values = modality_colors) +
  theme_classic() +
  xlab("Top colocalizing TWAS hits / Top xQTLs") +
  ylab("Modality")

p5 <- coloc_top_qtl_ratio |>
  ggplot(aes(x = coloc_top_n / 1000, y = modality, fill = modality)) +
  geom_boxplot(show.legend = FALSE, outlier.size = 0.3, linewidth = 0.4, median.linewidth = 0.5) +
  scale_fill_manual(values = modality_colors) +
  theme_classic() +
  xlab("Top colocalizing TWAS hits (×1000)") +
  ylab(NULL)

p6 <- coloc_top_qtl_ratio |>
  ggplot(aes(x = n_top_qtls / 1000, y = modality, fill = modality)) +
  geom_boxplot(show.legend = FALSE, outlier.size = 0.3, linewidth = 0.4, median.linewidth = 0.5) +
  scale_fill_manual(values = modality_colors) +
  theme_classic() +
  xlab("Top xQTLs (×1000)") +
  ylab(NULL)

(p1 + p2 + p3) / plot_spacer() / (p4 + p5 + p6) +
  plot_layout(heights = c(6, 1, 6)) +
  plot_annotation(tag_levels = "a")

ggsave("figures/figureS8.png", width = 10, height = 4, device = png)
