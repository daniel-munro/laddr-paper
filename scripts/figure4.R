## Residual latent RNA phenotypes

library(tidyverse)

modalities <- c(
  expression = "Expression",
  isoforms = "Isoform ratio",
  splicing = "Intron excision",
  alt_TSS = "Alt. TSS",
  alt_polyA = "Alt. polyA",
  stability = "RNA stability",
  latent_full = "Latent (full)",
  latent_residual = "Latent (residual)"
)

# Use muted version of Pantry colors to deemphasize what is already known
modality_colors <- c(
  Expression = "#bf4042",
  `Isoform ratio` = "#6a90cd",
  `Intron excision` = "#59a257",
  `Alt. TSS` = "#896090",
  `Alt. polyA` = "#d97f26",
  `RNA stability` = "#ddb23c",
  `Latent (full)` = "#13918d",
  `Latent (residual)` = "#1ce6df"
)

## Panel a: Latent-explicit correlations

latent_types = c(
  full = "Full",
  residual = "Residual",
  null = "Null"
)

latent_colors <- c(
  Full = "#13918d",
  Residual = "#1ce6df",
  Null = "white"
)

corrs_max <- read_tsv("data/processed/latent_explicit_corrs.tsv.gz", col_types = "ccccd") |>
  mutate(PC = PC |> str_replace("PC", "") |> fct_inorder(),
         r2_max = r^2)

corrs_max |>
  mutate(latent = factor(latent_types[latent], levels = latent_types)) |>
  ggplot(aes(x = PC, y = r2_max, fill = latent)) +
  geom_boxplot(outlier.size = 0.1, linewidth = 0.3) +
  scale_fill_manual(values = latent_colors) +
  theme_classic() +
  theme(
    axis.text = element_text(color = "black"),
    legend.position = "inside",
    legend.position.inside = c(0.8, 0.8),
  ) +
  xlab("Latent RNA phenotype rank per gene") +
  ylab(expression("Maximum "*r^2*" to explicit phenotype")) +
  labs(fill = "Latent type")

ggsave("figures/figure4/figure4a.png", width = 4.5, height = 3, device = png)

## Panel b: Explicit vs explicit + latent vs. full latent xQTLs

versions <- c(
  `residual-cross_pantry` = "Explicit",
  `residual-cross_latent` = "Explicit + Latent (residual)",
  `full-latent` = "Latent (full)"
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
  scale_fill_manual(values = modality_colors, guide = guide_legend(reverse = TRUE)) +
  theme_classic() +
  theme(
    axis.text = element_text(color = "black"),
    legend.key.size = unit(10, "pt"),
  ) +
  xlab("Independent cis-QTLs (×1000)") +
  ylab("RNA phenotypes") +
  labs(fill = "Modality")

ggsave("figures/figure4/figure4b.png", width = 6, height = 1.8, device = png)

## Panel c: Held-out modality xQTLs

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
  scale_y_continuous(breaks = c(1:6, 7.5), labels = rev(levels(qtls_held_out$held_out))) +
  scale_fill_manual(values = modality_colors, guide = guide_legend(reverse = TRUE)) +
  theme_classic() +
  theme(
    axis.text = element_text(color = "black"),
    legend.key.size = unit(10, "pt"),
  ) +
  xlab("Independent cis-QTLs (×1000)") +
  ylab("Modality held out        ") +
  labs(fill = "xQTL modality")

ggsave("figures/figure4/figure4c.png", width = 6, height = 2.3, device = png)
