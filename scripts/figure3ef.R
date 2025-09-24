library(tidyverse)

#############
## Panel e ## Latent vs. explicit TWAS gene-trait pairs by category
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
         phenos = fct_rev(phenos))

twas_topmod_count <- twas_both |>
  slice_min(twas_p, n = 1, with_ties = FALSE, by = c(phenos, trait, gene_id)) |>
  mutate(modality = factor(modalities[modality], levels = modalities)) |>
  count(phenos, category, modality)

twas_topmod_count |>
  ggplot(aes(x = category, y = n / 1000, fill = modality)) +
  facet_wrap(~ phenos) +
  geom_col() +
  scale_y_continuous(expand = expansion(mult = c(0, 0.04))) +
  scale_fill_manual(values = modality_colors) +
  theme_bw() +
  theme(
    axis.text = element_text(color = "black"),
    axis.text.x = element_text(hjust = 1, vjust = 1, angle = 45),
    legend.position = "inside",
    legend.position.inside = c(0.32, 0.7),
    legend.key.size = unit(9, "pt"),
    legend.text = element_text(size = 7),
    legend.title = element_text(size = 10),
    panel.grid = element_blank(),
  ) +
  labs(fill = "Modality of\npair's top hit") +
  xlab("Trait category") +
  ylab("Gene-trait pairs with TWAS hit(s)  (×1000) ")

ggsave("figures/figure3/figure3e.png", width = 4, height = 4, device = png)

#############
## Panel f ## Latent vs. explicit colocalizing TWAS hits
#############

twas_coloc_pairs <- twas_both |>
  filter(coloc_pp >= 0.8) |>
  distinct(phenos, trait, gene_id) |>
  count(phenos, trait) |>
  pivot_wider(id_cols = trait, names_from = phenos, values_from = n, values_fill = 0L) |>
  mutate(category = categories[trait])

twas_coloc_pairs |>
  ggplot(aes(x = `Knowledge-driven`, y = `Data-driven`, color = category)) +
  geom_abline(slope = 1, intercept = 0, color = "black", linewidth = 0.3) +
  geom_point() +
  annotate("text", x = 230, y = 150, label = "y = x", color = "black") +
  coord_fixed() +
  scale_x_log10() +
  scale_y_log10() +
  theme_bw() +
  theme(
    axis.text = element_text(color = "black"),
    legend.box.background = element_rect(),
    legend.position = "inside",
    legend.position.inside = c(0.745, 0.2),
    legend.key.size = unit(12, "pt"),
    legend.key.spacing = unit(3, "pt"),
    panel.grid = element_blank(),
  ) +
  xlab("Colocalizing gene-trait pairs, KPs") +
  ylab("Colocalizing gene-trait pairs, DPs") +
  labs(color = "Trait category")

ggsave("figures/figure3/figure3f.png", width = 4, height = 4, device = png)
