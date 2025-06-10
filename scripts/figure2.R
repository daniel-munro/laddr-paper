## Latent QTLs and effects of residualization

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

## Panel a: Full latent vs. explicit xQTLs
# (protein-coding only or all?)

qtls_gtextcga_full <- read_tsv("data/processed/gtextcga-full.qtls.tsv.gz", col_types = "cciccd")

qtls_pantry <- read_tsv("data/pantry/processed/gtex.comb.qtls.tsv.gz", col_types = "ccicccid")

pantry_pcg <- rtracklayer::import("data/pantry/Homo_sapiens.GRCh38.106.gtf.gz") |>
  as_tibble() |>
  filter(type == "gene",
         gene_biotype == "protein_coding")

pcg_shared <- read_tsv("data/processed/protein_coding_genes.tsv", col_types = "ccciic") |>
  filter(gene_id %in% pantry_pcg$gene_id)

qtls_per_gene <- bind_rows(
  qtls_gtextcga_full |>
    filter(gene_id %in% pcg_shared$gene_id) |>
    mutate(type = "Latent"),
  qtls_pantry |>
    filter(gene_id %in% pcg_shared$gene_id) |>
    mutate(type = "Explicit")
) |>
  summarise(qtls_per_gene = n() / nrow(genes),
            .by = c(type, tissue)) |>
  arrange(desc(type), desc(qtls_per_gene)) |>
  mutate(tissue = fct_inorder(tissue),
         type = fct_inorder(type) |> fct_rev())

ggplot(qtls_per_gene, aes(x = tissue, y = qtls_per_gene, fill = type)) +
  geom_col(position = "dodge") +
  scale_y_continuous(expand = c(0, 0)) +
  scale_fill_manual(values = c("coral", "#13918d")) +
  expand_limits(y = max(qtls_per_gene$qtls_per_gene * 1.03)) +
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
  ylab("xQTLs per gene") +
  labs(fill = NULL)

ggsave("figures/figure2/figure2a.png", width = 3.5, height = 3, device = png)

## Alternative: scatter plot

gtex_colors <- read_tsv(
  "data/pantry/gtex/tissueInfo.tsv",
  col_types = cols(tissueSiteDetailAbbr = "c", colorHex = "c", .default = "-")
) |>
  mutate(colorHex = str_c("#", colorHex)) |>
  deframe()

qtls_per_gene_wide <- qtls_per_gene |>
  pivot_wider(id_cols = tissue, names_from = "type", values_from = "qtls_per_gene")

ggplot(qtls_per_gene_wide, aes(x = Explicit, y = Latent, color = tissue)) +
  geom_abline(slope = 1, intercept = 0, linetype = 3) +
  geom_point(show.legend = FALSE) +
  expand_limits(
    x = c(0, max(qtls_per_gene_wide$Explicit) * 1.06),
    y = c(0, max(qtls_per_gene_wide$Latent) * 1.03),
  ) +
  scale_x_continuous(breaks = c(0, 0.5, 1)) +
  scale_color_manual(values = gtex_colors) +
  coord_fixed(expand = 0) +
  theme_bw() +
  theme(
    axis.text = element_text(color = "black"),
    panel.grid = element_blank(),
  ) +
  xlab("Explicit xQTLs per gene") +
  ylab("Latent xQTLs per gene") +
  labs(fill = NULL)

ggsave("figures/figure2/figure2a_scatter.png", width = 2.5, height = 4, device = png)

## Panel b: xQTLs by PC number

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

ggsave("figures/figure2/figure2b.png", width = 3.5, height = 3, device = png)

## Panel c: Latent-explicit correlations

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

ggsave("figures/figure2/figure2c.png", width = 4.5, height = 3, device = png)

## Panel d: Explicit vs explicit + latent vs. full latent xQTLs

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

ggsave("figures/figure2/figure2d.png", width = 6, height = 1.8, device = png)

## Panel e: Held-out modality xQTLs

load_qtls_held_out <- function(held_out) {
  fname <- str_glue("data/held_out/cross-no_{held_out}.cis_independent_qtl.txt.gz")
  df <- read_tsv(fname, col_types = "c-----c---------dc-i") |>
    rename(gene_id = group_id) |>
    separate_wider_delim(phenotype_id, ":",
                         names = c("modality", "phenotype_id"),
                         too_many = "merge"
    )
}

qtls_held_out <- tibble(held_out = names(modalities)[1:6]) |>
  reframe(load_qtls_held_out(held_out), .by = held_out) |>
  mutate(held_out = modalities[held_out], .before = 1) |>
  mutate(modality = factor(modalities[modality], levels = modalities))

qtls_held_out <- bind_rows(
  qtls_geuvadis |>
    filter(version == "Explicit + Latent (residual)") |>
    select(-version) |>
    mutate(held_out = "None held out", .before = 1),
  qtls_held_out
) |>
  mutate(held_out = fct_inorder(held_out))

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

ggsave("figures/figure2/figure2e.png", width = 6, height = 2.3, device = png)
