# Hybrid phenotype set TWAS

library(tidyverse)

#############
## Panel e ## Explicit + residual latent TWAS overlap
#############

modalities <- c(
  expression = "Expression",
  isoforms = "Isoform ratio",
  splicing = "Intron excision",
  alt_TSS = "Alt. TSS",
  alt_polyA = "Alt. polyA",
  stability = "RNA stability",
  latent = "Residual DD"
)

# Use muted version of Pantry colors to deemphasize what is already known
modality_colors <- c(
  Expression = "#bf4042",
  `Isoform ratio` = "#6a90cd",
  `Intron excision` = "#59a257",
  `Alt. TSS` = "#896090",
  `Alt. polyA` = "#d97f26",
  `RNA stability` = "#ddb23c",
  `Residual DD` = "#1ce6df"
)

categories <- read_tsv("data/pantry/geuvadis/twas/gwas_metadata.txt",
                       col_types = cols(Tag = "c", Category = "c", .default = "-")) |>
  mutate(Category = fct_lump_min(Category, 10)) |>
  deframe()

twas_pantry <- read_tsv("data/processed/geuvadis-pantry.twas_hits.tsv.gz", col_types = "cccc--dd")

twas_resid <- read_tsv("data/processed/geuvadis-residual.twas_hits.tsv.gz", col_types = "ccc--dd")

twas_panres <- bind_rows(
  twas_pantry |>
    select(trait, gene_id, modality, twas_p),
  twas_resid |>
    mutate(modality = "latent") |>
    select(trait, gene_id, modality, twas_p),
) |>
  mutate(category = categories[trait] |> fct_infreq(),
         modality = factor(modalities[modality], levels = modalities))

twas_panres_overlap <- twas_panres |>
  mutate(modality_type = if_else(modality == "Residual DD", "latent", "explicit")) |>
  distinct(trait, gene_id, modality_type) |>
  summarise(
    modality_hits = str_c(sort(unique(modality_type)), collapse = "_"),
    .by = c(trait, gene_id)
  ) |>
  mutate(
    modality_hits = c(explicit = "KP only",
                      explicit_latent = "KP & rDP",
                      latent = "rDP only")[modality_hits] |>
      fct_relevel("KP only", "KP & rDP", "rDP only"),
    category = fct_infreq(categories[trait]),
  )

twas_panres_overlap |>
  count(category, modality_hits) |>
  ggplot(aes(x = category, y = n / 1000, fill = modality_hits)) +
  geom_col(position = position_stack(reverse = TRUE), width = 0.7, color = "black") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.01))) +
  scale_fill_manual(values = c("white", "gray", "#444444")) +
  coord_flip() +
  theme_classic() +
  theme(
    axis.text = element_text(color = "black"),
    # axis.text.x = element_text(hjust = 1, vjust = 1, angle = 45),
    legend.position = "inside",
    legend.position.inside = c(0.8, 0.7),
    legend.key.size = unit(9, "pt"),
    legend.text = element_text(size = 7),
    legend.title = element_text(size = 9),
    panel.grid = element_blank(),
  ) +
  labs(fill = "Gene's xTWAS\nhits include") +
  xlab("Trait category") +
  ylab("Gene-trait pairs with TWAS hit(s) (×1000)")

ggsave("figures/figure4/figure4e.png", width = 5, height = 2, device = png)

#############
## Panel f ## Explicit + residual latent TWAS top hits
#############

twas_panres_tophit <- twas_panres |>
  slice_min(twas_p, n = 1, with_ties = FALSE, by = c(trait, gene_id))

twas_panres_tophit |>
  count(category, modality) |>
  ggplot(aes(x = category, y = n / 1000, fill = modality)) +
  geom_col(position = position_stack(reverse = TRUE), width = 0.7) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.01))) +
  scale_fill_manual(values = modality_colors) +
  coord_flip() +
  theme_classic() +
  theme(
    axis.text = element_text(color = "black"),
    # axis.text.x = element_text(hjust = 1, vjust = 1, angle = 45),
    legend.margin = margin_auto(0),
    legend.position = "inside",
    legend.position.inside = c(0.8, 0.6),
    legend.key.size = unit(8, "pt"),
    legend.text = element_text(size = 7),
    legend.title = element_text(size = 9),
    panel.grid = element_blank(),
  ) +
  labs(fill = "Modality of\npair's top hit") +
  xlab("Trait category") +
  ylab("Gene-trait pairs with TWAS hit(s) (×1000)")

ggsave("figures/figure4/figure4f.png", width = 5, height = 2, device = png)
