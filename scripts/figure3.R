## Latent xQTLs and TWAS

library(tidyverse)

## Panel a: Full latent vs. explicit xQTLs bar plot

qtls_gtextcga_full <- read_tsv("data/processed/gtextcga-full.qtls.tsv.gz", col_types = "cciccd")

qtls_pantry <- read_tsv("data/pantry/processed/gtex.comb.qtls.tsv.gz", col_types = "ccicccid")

qtls <- bind_rows(
  qtls_gtextcga_full |>
    mutate(type = "Latent"),
  qtls_pantry |>
    mutate(type = "Explicit")
) |>
  count(type, tissue) |>
  arrange(desc(type), desc(n)) |>
  mutate(tissue = fct_inorder(tissue),
         type = fct_inorder(type) |> fct_rev())

ggplot(qtls, aes(x = tissue, y = n / 1000, fill = type)) +
  geom_col(position = "dodge") +
  scale_y_continuous(expand = c(0, 0)) +
  scale_fill_manual(values = c("coral", "#13918d")) +
  expand_limits(y = max(qtls$n * 1.03 / 1000)) +
  theme_bw() +
  theme(
    axis.text.x = element_blank(),
    axis.text.y = element_text(color = "black"),
    axis.ticks.x = element_blank(),
    legend.key.size = unit(12, "pt"),
    legend.position = "inside",
    legend.position.inside = c(0.8, 0.8),
    panel.grid = element_blank(),
  ) +
  xlab("GTEx tissues") +
  ylab("xQTLs (×1000)") +
  labs(fill = NULL)

ggsave("figures/figure3/figure3a.png", width = 3.5, height = 3, device = png)

## Panel b: Full latent vs. explicit xQTLs scatter plot

gtex_colors <- read_tsv(
  "data/pantry/gtex/tissueInfo.tsv",
  col_types = cols(tissueSiteDetailAbbr = "c", colorHex = "c", .default = "-")
) |>
  mutate(colorHex = str_c("#", colorHex)) |>
  deframe()

qtls_wide <- qtls |>
  pivot_wider(id_cols = tissue, names_from = "type", values_from = "n")

ggplot(qtls_wide, aes(x = Explicit / 1000, y = Latent / 1000, color = tissue)) +
  geom_abline(slope = 1, intercept = 0, linetype = 3) +
  geom_point(show.legend = FALSE) +
  expand_limits(
    x = c(0, max(qtls_wide$Explicit) * 1.06 / 1000),
    y = c(0, max(qtls_wide$Latent) * 1.03 / 1000),
  ) +
  # scale_x_continuous(breaks = c(0, 0.5, 1)) +
  scale_y_continuous(breaks = c(0, 10, 20, 30, 40, 50, 60)) +
  scale_color_manual(values = gtex_colors) +
  coord_fixed(expand = 0) +
  theme_bw() +
  theme(
    axis.text = element_text(color = "black"),
    panel.grid = element_blank(),
  ) +
  xlab("Explicit xQTLs (×1000)") +
  ylab("Latent xQTLs (×1000)") +
  labs(fill = NULL)

ggsave("figures/figure3/figure3b.png", width = 2.5, height = 4, device = png)

## Panel c: xQTLs by PC number

qtls_gtextcga_full |>
  mutate(PC = str_split_i(phenotype_id, "__PC", 2) |>
           as.integer()) |>
  count(PC) |>
  mutate(mean_per_tissue = n / n_distinct(qtls_gtextcga_full$tissue)) |>
  ggplot(aes(x = PC, y = mean_per_tissue / 1000)) +
  geom_col(width = 0.75) +
  scale_x_continuous(breaks = 1:16) +
  theme_classic() +
  theme(
    axis.text = element_text(color = "black"),
  ) +
  ylab("Mean xQTLs per tissue (×1000)")

ggsave("figures/figure3/figure3c.png", width = 3.5, height = 3, device = png)

## Panel d: Latent vs. explicit TWAS hits by category

modalities <- c(
  expression = "Expression",
  isoforms = "Isoform ratio",
  splicing = "Intron excision",
  alt_TSS = "Alt. TSS",
  alt_polyA = "Alt. polyA",
  stability = "RNA stability",
  latent = "Latent"
)

# Use muted version of Pantry colors to deemphasize what is already known
modality_colors <- c(
  Expression = "#bf4042",
  `Isoform ratio` = "#6a90cd",
  `Intron excision` = "#59a257",
  `Alt. TSS` = "#896090",
  `Alt. polyA` = "#d97f26",
  `RNA stability` = "#ddb23c",
  `Latent` = "#13918d"
)

categories <- read_tsv("data/pantry/geuvadis/twas/gwas_metadata.txt",
                       col_types = cols(Tag = "c", Category = "c", .default = "-")) |>
  mutate(Category = fct_lump_min(Category, 10)) |>
  deframe()

twas <- read_tsv("data/processed/geuvadis-full.twas_hits.tsv.gz", col_types = "cccdddd")

twas_pantry <- read_tsv("data/pantry/processed/geuvadis.twas.tsv.gz", col_types = "cccc---d") |>
  filter(TWAS.P < 5e-8)

df <- bind_rows(
  twas |>
    select(trait, gene_id) |>
    mutate(trait = as.character(trait)) |>
    mutate(phenos = "Latent", modality = "latent", .before = 1),
  twas_pantry |>
    select(trait, gene_id, modality) |>
    mutate(phenos = "Explicit", .before = 1),
) |>
  mutate(
    category = categories[trait] |>
      fct_infreq(),
    modality = factor(modalities[modality], levels = modalities),
  )

df |>
  ggplot(aes(x = category, fill = modality)) +
  facet_wrap(~ phenos) +
  geom_bar() +
  scale_fill_manual(values = modality_colors) +
  theme_bw() +
  theme(
    axis.text.x = element_text(hjust = 1, vjust = 1, angle = 45),
    legend.position = "inside",
    legend.position.inside = c(0.32, 0.72),
    legend.key.size = unit(9, "pt"),
    legend.text = element_text(size = 7),
    legend.title = element_text(size = 10),
    panel.grid = element_blank(),
  ) +
  labs(fill = "Modality") +
  xlab("Trait category") +
  ylab("TWAS hits")

ggsave("figures/figure3/figure3d.png", width = 4, height = 4, device = png)
