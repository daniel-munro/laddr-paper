## Latent xQTLs and TWAS

library(tidyverse)

#############
## Panel a ## Full latent vs. explicit xQTLs bar plot
#############

qtls_gtextcga_full <- read_tsv("data/processed/gtextcga-full.qtls.tsv.gz", col_types = "cciccd")

qtls_pantry <- read_tsv("data/processed/gtex-pantry.qtls.tsv.gz", col_types = "ccicccd")

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

#############
## Panel b ## Full latent vs. explicit xQTLs scatter plot
#############

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

ggsave("figures/figure3/figure3b.png", width = 2.5, height = 3.5, device = png)

#############
## Panel c ## xQTLs by PC number
#############

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

twas_pantry <- read_tsv("data/processed/geuvadis-pantry.twas_hits.tsv.gz", col_types = "cccc--dd")

twas_both <- bind_rows(
  twas |>
    mutate(trait = as.character(trait)) |>
    mutate(phenos = "Latent", modality = "latent") |>
    select(phenos, trait, gene_id, modality, twas_p, coloc_pp),
  twas_pantry |>
    mutate(phenos = "Explicit") |>
    select(phenos, trait, gene_id, modality, twas_p, coloc_pp),
) |>
  mutate(category = categories[trait] |> fct_infreq())

twas_topmod <- twas_both |>
  slice_min(twas_p, n = 1, with_ties = FALSE, by = c(phenos, trait, gene_id)) |>
  mutate(modality = factor(modalities[modality], levels = modalities))

twas_topmod |>
  ggplot(aes(x = category, fill = modality)) +
  facet_wrap(~ phenos) +
  geom_bar() +
  scale_fill_manual(values = modality_colors) +
  theme_bw() +
  theme(
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
  ylab("Gene-trait pairs with TWAS hit(s)")

ggsave("figures/figure3/figure3d.png", width = 4, height = 4, device = png)

#############
## Panel e ## Latent vs. explicit colocalizing TWAS hits
#############

twas_coloc_pairs <- twas_both |>
df |>
  filter(coloc_pp >= 0.8) |>
  distinct(phenos, trait, gene_id) |>
  count(phenos, trait) |>
  pivot_wider(id_cols = trait, names_from = phenos, values_from = n, values_fill = 0L) |>
  mutate(category = categories[trait]) |>
  ggplot(aes(x = Explicit, y = Latent, color = category)) +
  geom_abline(slope = 1, intercept = 0, color = "#cccccc") +
  geom_point() +
  annotate("text", x = 230, y = 150, label = "y = x", color = "#999999") +
  coord_fixed() +
  scale_x_log10() +
  scale_y_log10() +
  theme_bw() +
  theme(
    legend.box.background = element_rect(),
    legend.position = "inside",
    legend.position.inside = c(0.745, 0.2),
    legend.key.size = unit(12, "pt"),
    legend.key.spacing = unit(3, "pt"),
    panel.grid = element_blank(),
  ) +
  xlab("Colocalizing gene-trait pairs, explicit phenos") +
  ylab("Colocalizing gene-trait pairs, latent phenos") +
  labs(color = "Trait category")

ggsave("figures/figure3/figure3e.png", width = 4, height = 4, device = png)

#############
## Panel f ## Explicit + residual latent TWAS top hits
#############

modalities2 <- c(
  expression = "Expression",
  isoforms = "Isoform ratio",
  splicing = "Intron excision",
  alt_TSS = "Alt. TSS",
  alt_polyA = "Alt. polyA",
  stability = "RNA stability",
  latent = "Latent (residual)"
)

# Use muted version of Pantry colors to deemphasize what is already known
modality_colors2 <- c(
  Expression = "#bf4042",
  `Isoform ratio` = "#6a90cd",
  `Intron excision` = "#59a257",
  `Alt. TSS` = "#896090",
  `Alt. polyA` = "#d97f26",
  `RNA stability` = "#ddb23c",
  `Latent (residual)` = "#13918d"
)

twas_resid <- read_tsv("data/processed/geuvadis-residual.twas_hits.tsv.gz", col_types = "ccccdddd")

twas_panres <- bind_rows(
  twas_pantry |>
    select(trait, gene_id, modality, twas_p),
  twas_resid |>
    mutate(modality = "latent") |>
    select(trait, gene_id, modality, twas_p),
) |>
  mutate(category = categories[trait] |> fct_infreq(),
         modality = factor(modalities2[modality], levels = modalities2))

twas_panres_tophit <- twas_panres |>
  slice_min(twas_p, n = 1, with_ties = FALSE, by = c(trait, gene_id))

twas_panres_tophit |>
  count(category, modality) |>
  ggplot(aes(x = category, y = n / 1000, fill = modality)) +
  geom_col() +
  scale_fill_manual(values = modality_colors2) +
  theme_bw() +
  theme(
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

ggsave("figures/figure3/figure3f.png", width = 2.3, height = 4, device = png)

#############
## Panel g ## Explicit + residual latent TWAS overlap
#############

twas_panres_overlap <- twas_panres |>
  mutate(modality_type = if_else(modality == "Latent (residual)", "Latent", "Explicit")) |>
  distinct(trait, gene_id, modality_type) |>
  summarise(
    modality_hits = str_c(sort(unique(modality_type)), collapse = "_"),
    .by = c(trait, gene_id)
  ) |>
  mutate(
    modality_hits = c(`Explicit` = "Explicit only",
                      `Explicit_Latent` = "Explicit & Latent",
                      `Latent` = "Latent only")[modality_hits] |>
      fct_relevel("Explicit only", "Explicit & Latent", "Latent only"),
    category = fct_infreq(categories[trait]),
  )

twas_panres_overlap |>
  count(category, modality_hits) |>
  ggplot(aes(x = category, y = n / 1000, fill = modality_hits)) +
  geom_col(color = "black") +
  scale_fill_manual(values = c("white", "gray", "#444444")) +
  theme_bw() +
  theme(
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

ggsave("figures/figure3/figure3g.png", width = 2.3, height = 4, device = png)
