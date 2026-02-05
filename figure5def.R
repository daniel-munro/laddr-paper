## Colocalization to xQTL ratios

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

qtls_rddp <- read_tsv(
  "data/processed/gtex-residual-cross.qtls.tsv.gz",
  col_types = "ccicccd"
) |>
  mutate(modality = factor(modalities[modality], levels = rev(modalities)))

twas <- bind_rows(
  read_tsv("data/processed/gtex-residual.twas_hits.tsv.gz", col_types = "cccc--dd") |>
    mutate(modality = "latent_residual", .after = "tissue"),
  read_tsv("data/processed/gtex-pantry.twas_hits.tsv.gz", col_types = "ccccc--dd")
)

coloc <- twas |>
  summarise(
    coloc_n = sum(coloc_pp > 0.8),
    coloc_frac = mean(coloc_pp > 0.8),
    .by = c(tissue, modality)
  ) |>
  mutate(modality = factor(modalities[modality], levels = rev(modalities)))

coloc_qtl_ratio <- qtls_rddp |>
  count(tissue, modality, name = "n_qtls") |>
  full_join(coloc, by = c("tissue", "modality"), relationship = "one-to-one") |>
  mutate(coloc_qtl_ratio = coloc_n / n_qtls)

p1 <- coloc_qtl_ratio |>
  ggplot(aes(x = coloc_qtl_ratio, y = modality, fill = modality)) +
  geom_boxplot(show.legend = FALSE, outlier.size = 0.3, linewidth = 0.4, median.linewidth = 0.5) +
  scale_fill_manual(values = modality_colors) +
  theme_classic() +
  xlab("No. Colocalizing TWAS hits / No. xQTLs") +
  ylab("Modality") +
  ggtitle(label = NULL, subtitle = "Relative colocalization ratio")

p2 <- coloc_qtl_ratio |>
  ggplot(aes(x = coloc_n / 1000, y = modality, fill = modality)) +
  geom_boxplot(show.legend = FALSE, outlier.size = 0.3, linewidth = 0.4, median.linewidth = 0.5) +
  scale_fill_manual(values = modality_colors) +
  theme_classic() +
  xlab("Colocalizing TWAS hits (×1000)") +
  ylab(NULL) +
  ggtitle(label = NULL, subtitle = "Numerator")

p3 <- coloc_qtl_ratio |>
  ggplot(aes(x = n_qtls / 1000, y = modality, fill = modality)) +
  geom_boxplot(show.legend = FALSE, outlier.size = 0.3, linewidth = 0.4, median.linewidth = 0.5) +
  scale_fill_manual(values = modality_colors) +
  theme_classic() +
  xlab("xQTLs (×1000)") +
  ylab(NULL) +
  ggtitle(label = NULL, subtitle = "Denominator")

p1 + p2 + p3

ggsave("figures/figure5/figure5def.png", width = 10, height = 2, device = png)
