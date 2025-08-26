## Residual latent RNA phenotypes

library(tidyverse)

modalities <- c(
  expression = "Expression",
  isoforms = "Isoform ratio",
  splicing = "Intron excision",
  alt_TSS = "Alt. TSS",
  alt_polyA = "Alt. polyA",
  stability = "RNA stability",
  latent_residual = "Residual DD",
  latent_full = "Data-driven"
)

# Use muted version of Pantry colors to deemphasize what is already known
modality_colors <- c(
  Expression = "#bf4042",
  `Isoform ratio` = "#6a90cd",
  `Intron excision` = "#59a257",
  `Alt. TSS` = "#896090",
  `Alt. polyA` = "#d97f26",
  `RNA stability` = "#ddb23c",
  `Residual DD` = "#1ce6df",
  `Data-driven` = "#13918d"
)

#############
## Panel a ## Latent-explicit correlations
#############

latent_types = c(
  full = "Data-driven",
  residual = "Residual data-driven",
  null = "DP shuffled (control)"
)

latent_colors <- c(
  `Data-driven` = "#13918d",
  `Residual data-driven` = "#1ce6df",
  `DP shuffled (control)` = "#aaaaaa"
)

corrs_max <- read_tsv("data/processed/latent_explicit_corrs.tsv.gz", col_types = "ccccd") |>
  mutate(PC = PC |> str_replace("PC", "") |> fct_inorder(),
         r2_max = r^2)

corrs_max_stats <- corrs_max |>
  summarise(
    r2_max_05 = quantile(r2_max, 0.05),
    r2_max_median = median(r2_max),
    r2_max_95 = quantile(r2_max, 0.95),
    .by = c(latent, PC)
  )

corrs_max_stats |>
  mutate(latent = factor(latent_types[latent], levels = latent_types)) |>
  ggplot(aes(x = PC,
             y = r2_max_median,
             ymin = r2_max_05,
             ymax = r2_max_95,
             color = latent)) +
  geom_pointrange(linewidth = 0.8, size = 0.3, position = position_dodge(width = 0.6)) +
  scale_y_continuous(expand = c(0, 0)) +
  expand_limits(y = 0.85) +
  scale_color_manual(values = latent_colors) +
  theme_classic() +
  theme(
    axis.text = element_text(color = "black"),
    legend.position = "inside",
    legend.position.inside = c(0.8, 0.8),
  ) +
  xlab("Data-driven phenotype rank per gene") +
  ylab(expression("Maximum "*r^2*" to a KP")) +
  labs(color = "Phenotype type")

ggsave("figures/figure4/figure4a.png", width = 4.5, height = 3, device = png)

#############
## Panel b ## Explicit vs explicit + latent vs. full latent xQTLs
#############

versions <- c(
  `residual-cross_pantry` = "KP",
  `residual-cross_latent` = "Hybrid (KP + RP)",
  `full-latent` = "DP"
)

qtls_geuvadis <- read_tsv("data/processed/geuvadis.qtls.tsv.gz", col_types = "ccccdci") |>
  mutate(modality = factor(modalities[modality], levels = names(modality_colors)),
         version = factor(versions[version], levels = versions))

qtls_geuvadis |>
  count(version, modality) |>
  mutate(modality = fct_rev(modality),
         version = fct_rev(version)) |>
  ggplot(aes(x = n / 1000, y = version, fill = modality)) +
  geom_col(width = 0.8) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_fill_manual(values = modality_colors, guide = guide_legend(reverse = TRUE)) +
  theme_classic() +
  theme(
    axis.text = element_text(color = "black"),
    legend.key.size = unit(10, "pt"),
  ) +
  xlab("Independent cis-xQTLs (×1000)") +
  ylab("RNA phenotypes") +
  labs(fill = "Modality")

ggsave("figures/figure4/figure4b.png", width = 6, height = 1.8, device = png)

#############
## Panel c ## Held-out modality xQTLs
#############

qtls_held_out <- read_tsv(
  "data/processed/held_out-geuvadis.qtls.tsv.gz", col_types = "ccccdci"
) |>
  mutate(held_out = c(modalities, none = "None held out")[held_out] |>
           fct_inorder(),
         modality = factor(modalities[modality], levels = modalities))

# Convert held_out factor to numeric to add space between "none" and others
qtls_held_out |>
  count(held_out, modality) |>
  mutate(held_out = fct_rev(held_out) |> as.integer(),
         held_out = if_else(held_out == 7, 7.5, held_out),
         modality = fct_rev(modality)) |>
  ggplot(aes(x = n / 1000, y = held_out, fill = modality)) +
  geom_col(orientation = "y") +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(breaks = c(1:6, 7.5), labels = rev(levels(qtls_held_out$held_out))) +
  scale_fill_manual(values = modality_colors, guide = guide_legend(reverse = TRUE)) +
  theme_classic() +
  theme(
    axis.text = element_text(color = "black"),
    legend.key.size = unit(10, "pt"),
  ) +
  xlab("Independent cis-xQTLs (×1000)") +
  ylab("Modality held out        ") +
  labs(fill = "Modality")

ggsave("figures/figure4/figure4c.png", width = 6, height = 2.3, device = png)

#############
## Panel d ## Explicit + residual latent TWAS overlap
#############

modalities2 <- c(
  expression = "Expression",
  isoforms = "Isoform ratio",
  splicing = "Intron excision",
  alt_TSS = "Alt. TSS",
  alt_polyA = "Alt. polyA",
  stability = "RNA stability",
  latent = "Residual DD"
)

# Use muted version of Pantry colors to deemphasize what is already known
modality_colors2 <- c(
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
         modality = factor(modalities2[modality], levels = modalities2))

twas_panres_overlap <- twas_panres |>
  mutate(modality_type = if_else(modality == "Residual DD", "latent", "explicit")) |>
  distinct(trait, gene_id, modality_type) |>
  summarise(
    modality_hits = str_c(sort(unique(modality_type)), collapse = "_"),
    .by = c(trait, gene_id)
  ) |>
  mutate(
    modality_hits = c(explicit = "KP only",
                      explicit_latent = "KP & RP",
                      latent = "RP only")[modality_hits] |>
      fct_relevel("KP only", "KP & RP", "RP only"),
    category = fct_infreq(categories[trait]),
  )

twas_panres_overlap |>
  count(category, modality_hits) |>
  ggplot(aes(x = category, y = n / 1000, fill = modality_hits)) +
  geom_col(color = "black") +
  scale_y_continuous(expand = c(0, 0)) +
  expand_limits(y = 6.4) +
  scale_fill_manual(values = c("white", "gray", "#444444")) +
  theme_classic() +
  theme(
    axis.text = element_text(color = "black"),
    axis.text.x = element_text(hjust = 1, vjust = 1, angle = 45),
    legend.position = "inside",
    legend.position.inside = c(0.68, 0.7),
    legend.key.size = unit(9, "pt"),
    legend.text = element_text(size = 7),
    legend.title = element_text(size = 10),
    panel.grid = element_blank(),
  ) +
  labs(fill = "Gene's xTWAS\nhits include") +
  xlab("Trait category") +
  ylab("Gene-trait pairs with TWAS hit(s) (×1000)")

ggsave("figures/figure4/figure4d.png", width = 2.3, height = 4, device = png)

#############
## Panel e ## Explicit + residual latent TWAS top hits
#############

twas_panres_tophit <- twas_panres |>
  slice_min(twas_p, n = 1, with_ties = FALSE, by = c(trait, gene_id))

twas_panres_tophit |>
  count(category, modality) |>
  ggplot(aes(x = category, y = n / 1000, fill = modality)) +
  geom_col() +
  scale_y_continuous(expand = c(0, 0)) +
  expand_limits(y = 6.4) +
  scale_fill_manual(values = modality_colors2) +
  theme_classic() +
  theme(
    axis.text = element_text(color = "black"),
    axis.text.x = element_text(hjust = 1, vjust = 1, angle = 45),
    legend.position = "inside",
    legend.position.inside = c(0.7, 0.7),
    legend.key.size = unit(9, "pt"),
    legend.text = element_text(size = 7),
    legend.title = element_text(size = 10),
    panel.grid = element_blank(),
  ) +
  labs(fill = "Modality of\npair's top hit") +
  xlab("Trait category") +
  ylab("Gene-trait pairs with TWAS hit(s) (×1000)")

ggsave("figures/figure4/figure4e.png", width = 2.3, height = 4, device = png)
