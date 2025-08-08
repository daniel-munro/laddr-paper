## Latent RNA phenotype characteristics

library(tidyverse)
library(patchwork)

hsq_latent <- read_tsv("data/twas/Geuvadis-latent.profile", col_types = "cidddddd-dddd-d") |>
  separate_wider_delim(id, "__", names = c("gene_id", "PC"), cols_remove = FALSE) |>
  mutate(PC = str_replace(PC, "PC", "") |> as.integer())

#############
## Panel a ## Correlation heatmap example
#############

# gene <- "ENSG00000170291"  # ELP5
gene <- "ENSG00000090238"  # YPEL3
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
  "data/phenos/geuvadis-full/Geuvadis-latent.bed.gz",
  col_types = cols(`#chr` = "-", start = "-", end = "-", phenotype_id = "c", .default = "d")
) |>
  separate_wider_delim(phenotype_id, "__", names = c("gene_id", "PC")) |>
  filter(gene_id == gene) |>
  select(-gene_id) |>
  column_to_rownames("PC")

phenos_pantry <- tibble(modality = names(modalities)) |>
  reframe(
    read_tsv(
      str_glue("data/pantry_phenos/geuvadis/{modality}.bed.gz"),
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

cor(t(phenos_latent), t(phenos_pantry), method = "spearman") |>
  as_tibble(rownames = "PC") |>
  pivot_longer(-PC, names_to = "pantry_pheno", values_to = "rho") |>
  mutate(PC = fct_inorder(PC) |> fct_rev(),
         pantry_pheno = fct_inorder(pantry_pheno)) |>
  ggplot(aes(x = pantry_pheno, y = PC, fill = rho)) +
  geom_tile() +
  coord_fixed(expand = 0) +
  scale_x_discrete(position = "top") +
  scale_fill_gradient2() +
  theme_bw() +
  theme(
    axis.text.x = element_text(hjust = 0, angle = 60, color = "black"),
    axis.text.y = element_text(color = "black"),
    legend.box.spacing = unit(30, "pt"),
    panel.grid = element_blank(),
  ) +
  xlab(str_glue("Explicit phenotypes for {gene_name}")) +
  ylab(str_glue("Latent phenotypes for {gene_name}")) +
  labs(fill = expression("Corr. "*(rho)))

ggsave("figures/figure2/figure2a.png", width = 4.8, height = 4, device = png)

## Show significant heritability values to the right

hsq_latent |>
  filter(gene_id == gene) |>
  arrange(PC) |>
  select(gene_id, PC, hsq, hsq.se, hsq.pv)

#############
## Panel b ## Number of heritable phenotypes per PC
#############

hsq_latent |>
  count(PC) |>
  ggplot(aes(x = PC, y = n / 1000)) +
  geom_col(width = 0.7, fill = "black") +
  scale_x_continuous(expand = c(0, 0.2)) +
  scale_y_continuous(expand = c(0, 0)) +
  expand_limits(y = 10.4) +
  theme_bw() +
  theme(
    axis.text = element_text(color = "black"),
    panel.grid = element_blank(),
  ) +
  ylab(expression("Genes with significant cis-"*h^2*" (×1000)"))

ggsave("figures/figure2/figure2b.png", width = 2, height = 3, device = png)

#############
## Panel c ## Latent vs explicit phenotype cis-heritability
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

hsq_pantry <- read_tsv("data/pantry/processed/geuvadis.hsq.tsv.gz", col_types = "ccciddddddddddd") |>
  mutate(modality = factor(modalities[modality], levels = modalities) |>
           fct_reorder(hsq, .desc = TRUE))

p1 <- hsq_latent |>
  mutate(PC = as.character(PC) |>
           fct_reorder(PC) |>
           fct_lump_n(n = 10, other_level = "11+")) |>
  ggplot(aes(x = PC, y = hsq)) +
  geom_boxplot() +
  scale_y_log10(expand = c(0, 0), minor_breaks = c(0.01, 0.02, 0.03)) +
  expand_limits(y = c(0.0099, 1.01)) +
  theme_bw() +
  theme(
    axis.text = element_text(color = "black"),
    panel.grid = element_blank(),
    plot.margin = margin(5.5, 0, 5.5, 5.5),
  ) +
  ylab(expression(h^2)) +
  ggtitle("Latent phenotypes")

p2 <- hsq_pantry |>
  ggplot(aes(x = modality, y = hsq)) +
  geom_boxplot() +
  scale_y_log10(expand = c(0, 0)) +
  expand_limits(y = c(0.0099, 1.01)) +
  theme_bw() +
  theme(
    axis.text.x = element_text(hjust = 1, vjust = 1, angle = 45, color = "black"),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    panel.grid = element_blank(),
    plot.margin = margin(5.5, 5.5, 5.5, 0),
  ) +
  xlab("Modality") +
  ylab(NULL) +
  ggtitle("Explicit phenotypes")

p1 + p2 + plot_layout(widths = c(11, 6))

ggsave("figures/figure2/figure2c.png", width = 5, height = 4, device = png)

