library(tidyverse)

#############
## Panel d ## Latent vs. explicit TWAS gene-trait pairs by category
#############

modalities <- c(
  expression = "Expression",
  isoforms = "Isoform ratio",
  splicing = "Intron excision",
  alt_TSS = "Alt. TSS",
  alt_polyA = "Alt. polyA",
  stability = "RNA stability",
  latent = "Data-driven"
)

# Use muted version of Pantry colors to deemphasize what is already known
modality_colors <- c(
  Expression = "#bf4042",
  `Isoform ratio` = "#6a90cd",
  `Intron excision` = "#59a257",
  `Alt. TSS` = "#896090",
  `Alt. polyA` = "#d97f26",
  `RNA stability` = "#ddb23c",
  `Data-driven` = "#13918d"
)

categories <- read_tsv("data/pantry/geuvadis/twas/gwas_metadata.txt",
                       col_types = cols(Tag = "c", Category = "c", .default = "-")) |>
  mutate(Category = fct_lump_min(Category, 10)) |>
  deframe()

twas <- read_tsv("data/processed/geuvadis-full.twas_hits.tsv.gz", col_types = "cccdddd")

twas_pantry <- read_tsv("data/processed/geuvadis-pantry.twas_hits.tsv.gz", col_types = "cccc--dd")

twas_both <- bind_rows(
  twas |>
    mutate(trait = as.character(trait)) |>
    mutate(phenos = "Data-driven", modality = "latent") |>
    select(phenos, trait, gene_id, modality, twas_p, coloc_pp),
  twas_pantry |>
    mutate(phenos = "Knowledge-driven") |>
    select(phenos, trait, gene_id, modality, twas_p, coloc_pp),
) |>
  mutate(category = categories[trait] |> fct_infreq(),
         phenos = phenos)

twas_topmod_count <- twas_both |>
  slice_min(twas_p, n = 1, with_ties = FALSE, by = c(phenos, trait, gene_id)) |>
  mutate(modality = factor(modalities[modality], levels = modalities)) |>
  count(phenos, category, modality)

twas_topmod_count |>
  ggplot(aes(y = category, x = n / 1000, fill = modality)) +
  facet_wrap(~ phenos, ncol = 1) +
  geom_col(width = 0.6) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.04))) +
  scale_fill_manual(values = modality_colors) +
  theme_bw() +
  theme(
    axis.text = element_text(color = "black"),
    legend.position = "inside",
    legend.position.inside = c(0.75, 0.25),
    legend.key.size = unit(8, "pt"),
    legend.margin = margin_auto(0),
    legend.spacing.y = unit(0, "pt"),
    legend.text = element_text(size = 8),
    legend.title = element_text(size = 10, margin = margin_auto(b = 2, unit = "pt")),
    panel.grid = element_blank(),
  ) +
  labs(fill = "Modality of pair's top hit") +
  ylab("Trait category") +
  xlab("Gene-trait pairs with TWAS hit(s)  (×1000) ")

ggsave("figures/figure3/figure3d.png", width = 5, height = 3.3, device = png)

#############
## Panel e ## Latent vs. explicit colocalizing TWAS hits
#############

category_colors <- c(
  Anthropometric = "#1b9e77",
  Blood = "#d95f02",
  Cardiometabolic = "#7570b3",
  Immune = "#e7298a",
  `Psychiatric-neurologic` = "#66a61e",
  Other = "#888888"
)

twas_coloc_pairs <- twas_both |>
  filter(coloc_pp >= 0.8) |>
  distinct(phenos, trait, gene_id) |>
  count(phenos, trait) |>
  pivot_wider(id_cols = trait, names_from = phenos, values_from = n, values_fill = 0L) |>
  mutate(category = categories[trait])

twas_coloc_pairs |>
  ggplot(aes(x = `Knowledge-driven`, y = `Data-driven`, color = category)) +
  geom_abline(slope = 1, intercept = 0, color = "black", linewidth = 0.3) +
  geom_point(alpha = 0.8) +
  annotate("text", x = 230, y = 140, label = "y = x", color = "black") +
  coord_fixed() +
  scale_x_log10() +
  scale_y_log10() +
  scale_color_manual(values = category_colors) +
  theme_bw() +
  theme(
    axis.text = element_text(color = "black"),
    legend.box.background = element_rect(),
    legend.margin = margin_auto(2, "pt"),
    legend.position = "inside",
    legend.position.inside = c(0.72, 0.2),
    legend.key.size = unit(10, "pt"),
    legend.key.spacing = unit(0, "pt"),
    panel.grid = element_blank(),
  ) +
  xlab("Coloc. gene-trait pairs, KDPs") +
  ylab("Coloc. gene-trait pairs, DDPs") +
  labs(color = "Trait category")

ggsave("figures/figure3/figure3e.png", width = 3.5, height = 3.5, device = png)
