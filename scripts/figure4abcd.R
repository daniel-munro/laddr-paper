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
  null = "DD shuffled (control)"
)

latent_colors <- c(
  `Data-driven` = "#13918d",
  `Residual data-driven` = "#1ce6df",
  `DD shuffled (control)` = "#aaaaaa"
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
  scale_y_continuous(expand = expansion(mult = c(0, 0.02))) +
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
  `residual-cross_latent` = "Hybrid (KP + rDP)",
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
## Panel c ## GTEx xQTLs
#############

qtls_gtex_dp <- read_tsv("data/processed/gtextcga-full.qtls.tsv.gz", col_types = "cciccd")
qtls_gtex_hybrid <- read_tsv("data/processed/gtex-residual-cross.qtls.tsv.gz", col_types = "ccicccd")
qtls_gtex_kp <- read_tsv("data/processed/gtex-pantry.qtls.tsv.gz", col_types = "ccicccd")

qtls_gtex_all <- bind_rows(
  qtls_gtex_dp |> mutate(modality = "latent_full", mapping = "Data-driven"),
  qtls_gtex_hybrid |> mutate(mapping = "Hybrid (KP + rDP)"),
  qtls_gtex_kp |> mutate(mapping = "Knowledge-driven"),
)

tissue_order <- qtls_gtex_all |>
  count(tissue, sort = TRUE) |>
  pull(tissue)

qtls_gtex_all |>
  count(mapping, tissue, modality) |>
  mutate(modality = factor(modalities[modality], levels = modalities) |> fct_rev(),
         mapping = fct_rev(mapping),
         tissue = factor(tissue, levels = tissue_order)) |>
  ggplot(aes(x = tissue, y = n / 1000, fill = modality)) +
  facet_wrap(~mapping) +
  geom_col(width = 1, show.legend = FALSE) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.04))) +
  scale_fill_manual(values = modality_colors) +
  theme_classic() +
  theme(
    axis.text.x = element_blank(),
    axis.text.y = element_text(color = "black"),
    axis.ticks.x = element_blank(),
  ) +
  xlab("GTEx tissues ordered by total count across panels") +
  ylab("Independent cis-xQTLs (×1000)")

ggsave("figures/figure4/figure4c.png", width = 5.5, height = 2.3, device = png)

#############
## Panel d ## Held-out modality xQTLs
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

ggsave("figures/figure4/figure4d.png", width = 6, height = 2.3, device = png)

