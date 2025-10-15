## Latent RNA phenotype characteristics

library(tidyverse)
library(patchwork)

modalities <- c(
  expression = "Expression",
  isoforms = "Isoform ratio",
  splicing = "Intron excision",
  alt_TSS = "Alt. TSS",
  alt_polyA = "Alt. polyA",
  stability = "RNA stability",
  latent = "Latent"
)

modality_colors <- c(
  Expression = "#bf4042",
  `Isoform ratio` = "#6a90cd",
  `Intron excision` = "#59a257",
  `Alt. TSS` = "#896090",
  `Alt. polyA` = "#d97f26",
  `RNA stability` = "#ddb23c",
  Latent = "#13918d"
)

#############
## Panel a ## Correlation heatmap example
#############

# gene <- "ENSG00000170291"  # ELP5
# gene <- "ENSG00000090238"  # YPEL3
gene <- "ENSG00000223547"  # ZNF844
gene_names <- read_tsv("data/processed/pcg_and_lncrna.tsv", col_types = "cc-----") |>
  deframe()
gene_name <- gene_names[gene]

modalities <- c(
  expression = "Expression",
  isoforms = "Isoform ratio",
  splicing = "Intron excision",
  alt_TSS = "Alt. TSS",
  alt_polyA = "Alt. polyA",
  stability = "RNA stability"
)

phenos_latent <- read_tsv(
  "data/phenos/gtextcga-full/NERVET-latent.bed.gz",
  col_types = cols(`#chr` = "-", start = "-", end = "-", phenotype_id = "c", .default = "d")
) |>
  separate_wider_delim(phenotype_id, "__", names = c("gene_id", "PC")) |>
  filter(gene_id == gene) |>
  select(-gene_id) |>
  column_to_rownames("PC")

phenos_pantry <- tibble(modality = names(modalities)) |>
  reframe(
    read_tsv(
      str_glue("data/pantry_phenos/NERVET/{modality}.bed.gz"),
      col_types = cols(`#chr` = "-", start = "-", end = "-", phenotype_id = "c", .default = "d")
    ) |>
      mutate(gene_id = str_replace(phenotype_id, "__.+", ""), .before = 1) |>
      filter(gene_id == gene) |>
      select(-gene_id),
    .by = modality
  ) |>
  mutate(
    phenotype_id = if (n() == 1) {
      modalities[modality]
    } else {
      str_glue("{modalities[modality]} {seq_len(n())}")
    },
    .by = modality
  ) |>
  select(-modality) |>
  column_to_rownames("phenotype_id")

stopifnot(identical(colnames(phenos_latent), colnames(phenos_pantry)))

phenos_corr <- cor(t(phenos_latent), t(phenos_pantry), method = "spearman") |>
  as_tibble(rownames = "PC") |>
  pivot_longer(-PC, names_to = "pantry_pheno", values_to = "rho") |>
  mutate(PC = str_replace(PC, "PC", "DDP") |> fct_inorder() |> fct_rev(),
         pantry_pheno = fct_inorder(pantry_pheno),
         rho_abs = abs(rho))

ggplot(phenos_corr, aes(x = pantry_pheno, y = PC, fill = rho_abs)) +
  geom_tile() +
  coord_fixed(expand = 0) +
  scale_x_discrete(position = "top") +
  scale_fill_gradient(low = "white", high = "black", breaks = c(0, 1)) +
  expand_limits(fill = c(0, 1)) +
  theme_bw() +
  theme(
    axis.text.x = element_text(hjust = 0, angle = 60, color = "black"),
    axis.text.y = element_text(color = "black"),
    legend.background = element_rect(color = "black", linewidth = 0.3),
    legend.direction = "horizontal",
    legend.justification = c(1, 0),
    legend.margin = margin_auto(4),
    legend.key.size = unit(12, "pt"),
    legend.position = "inside",
    legend.position.inside = c(0.95, 0.05),
    panel.grid = element_blank(),
    plot.margin = margin_auto(3),
  ) +
  xlab(str_glue("Knowledge-driven phenotypes for {gene_name}")) +
  ylab(str_glue("Data-driven phenotypes for {gene_name}")) +
  labs(fill = expression("|"*rho*"|"))

ggsave("figures/figure2/figure2a.png", width = 4, height = 4, device = png)

#############
## Panel b ## Number of heritable phenotypes per PC
#############

hsq_latent <- read_tsv("data/twas/Geuvadis-latent.profile", col_types = "cidddddd-dddd-d") |>
  separate_wider_delim(id, "__", names = c("gene_id", "PC"), cols_remove = FALSE) |>
  mutate(PC = str_replace(PC, "PC", "") |> as.integer())

hsq_pantry <- read_tsv("data/pantry/processed/geuvadis.hsq.tsv.gz", col_types = "ccciddddddddddd") |>
  mutate(modality = factor(modalities[modality], levels = modalities) |>
           fct_infreq())

p1 <- hsq_latent |>
  count(PC) |>
  ggplot(aes(x = PC, y = n / 1000)) +
  geom_col(width = 0.7, fill = modality_colors["Latent"]) +
  scale_x_continuous(expand = c(0, 0.3)) +
  scale_y_continuous(expand = c(0, 0)) +
  expand_limits(y = 10.5) +
  theme_bw() +
  theme(
    axis.text = element_text(color = "black"),
    panel.grid = element_blank(),
    plot.margin = margin_auto(0),
    plot.title = element_text(hjust = 0.5, size = 11),
  ) +
  xlab("DDP rank in gene") +
  ylab(expression("Phenotypes with cis-"*h^2*" P<0.01 (×1000)")) +
  ggtitle("Data-driven")

p2 <- hsq_pantry |>
  count(modality) |>
  ggplot(aes(x = modality, y = n / 1000, fill = modality)) +
  geom_col(width = 0.4, show.legend = FALSE) +
  scale_y_continuous(expand = c(0, 0)) +
  expand_limits(y = 10.5) +
  scale_fill_manual(values = modality_colors) +
  theme_bw() +
  theme(
    axis.text.y = element_blank(),
    axis.text.x = element_text(hjust = 1, vjust = 1, angle = 45, color = "black"),
    axis.ticks.y = element_blank(),
    panel.grid = element_blank(),
    plot.margin = margin_auto(l = 0, r = 2),
    plot.title = element_text(hjust = 0.5, size = 11),
  ) +
  xlab("Modality") +
  ylab(NULL) +
  ggtitle("Knowledge-driven")

p1 + p2 + plot_layout(widths = c(8, 6))

ggsave("figures/figure2/figure2b.png", width = 3.25, height = 4, device = png)

#############
## Panel c ## Latent vs explicit phenotype cis-heritability
#############

p3 <- hsq_latent |>
  mutate(PC = as.character(PC) |>
           fct_reorder(PC) |>
           fct_lump_n(n = 9, other_level = "10+")) |>
  ggplot(aes(x = PC, y = hsq)) +
  geom_boxplot(width = 0.55, outlier.size = 0.2, outlier.alpha = 0.5,
               fill = fill_alpha(modality_colors["Latent"], 0.4)) +
  scale_y_log10(expand = c(0, 0), minor_breaks = c(0.01, 0.02, 0.03)) +
  expand_limits(y = c(0.0099, 1.01)) +
  theme_bw() +
  theme(
    axis.text = element_text(color = "black"),
    panel.grid = element_blank(),
    plot.margin = margin_auto(0),
    plot.title = element_text(hjust = 0.5, size = 11),
  ) +
  xlab("DDP rank in gene") +
  ylab(expression("cis-"*h^2)) +
  ggtitle("Data-driven")

p4 <- hsq_pantry |>
  ggplot(aes(x = modality, y = hsq, fill = modality)) +
  geom_boxplot(width = 0.5, outlier.size = 0.2, outlier.alpha = 0.5, show.legend = FALSE) +
  scale_y_log10(expand = c(0, 0)) +
  scale_fill_manual(values = alpha(modality_colors, 0.4)) +
  expand_limits(y = c(0.0099, 1.01)) +
  theme_bw() +
  theme(
    axis.text.x = element_text(hjust = 1, vjust = 1, angle = 45, color = "black"),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    panel.grid = element_blank(),
    plot.margin = margin_auto(l = 0, r = 2),
    plot.title = element_text(hjust = 0.5, size = 11),
  ) +
  xlab("Modality") +
  ylab(NULL) +
  ggtitle("Knowledge-driven")

p3 + p4 + plot_layout(widths = c(9, 6))

ggsave("figures/figure2/figure2c.png", width = 3.25, height = 4, device = png)

