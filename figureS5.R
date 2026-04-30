## MAJIQ retained intron matches

library(tidyverse)
library(patchwork)
library(scales)

modalities <- c(
  expression = "Expression",
  isoforms = "Isoform ratio",
  splicing = "Intron excision",
  alt_TSS = "Alt. TSS",
  alt_polyA = "Alt. polyA",
  stability = "RNA stability",
  latent_residual = "Latent (residual)"
)

modality_colors <- c(
  Expression = "#bf4042",
  `Isoform ratio` = "#6a90cd",
  `Intron excision` = "#59a257",
  `Alt. TSS` = "#896090",
  `Alt. polyA` = "#d97f26",
  `RNA stability` = "#ddb23c",
  `Latent (residual)` = "#1ce6df"
)

match_group_colors <- c(modality_colors, `Multiple modalities` = "#7a4517", None = "#cccccc")

summary_metrics <- read_tsv("data/majiq/summary.tsv", show_col_types = FALSE) |>
  mutate(value = as.integer(value))

retained_introns <- summary_metrics |>
  filter(metric == "retained_introns") |>
  pull(value)

laddr_matches <- read_tsv("data/majiq/high_correlation_ir_laddr.tsv.gz", show_col_types = FALSE) |>
  transmute(
    majiq_ir_id,
    phenotype_id = paste0("latent_residual:", phenotype_id),
    gene_id = laddr_gene_id,
    modality = "latent_residual"
  )

kdp_matches <- read_tsv("data/majiq/high_correlation_ir_kdp.tsv.gz", show_col_types = FALSE) |>
  transmute(
    majiq_ir_id,
    phenotype_id,
    gene_id = phenotype_gene_id,
    modality = str_extract(phenotype_id, "^[^:]+")
  )

ir_matches <- bind_rows(laddr_matches, kdp_matches) |>
  mutate(modality_label = modalities[modality])

single_modality_counts <- ir_matches |>
  distinct(majiq_ir_id, modality_label) |>
  summarise(
    n_modalities = n(),
    modality_label = first(sort(modality_label)),
    .by = majiq_ir_id
  ) |>
  mutate(
    match_group = if_else(n_modalities == 1, modality_label, "Multiple modalities")
  ) |>
  count(match_group, name = "retained_introns")

match_breakdown <- bind_rows(
  single_modality_counts,
  tibble(
    match_group = "None",
    retained_introns = retained_introns - n_distinct(ir_matches$majiq_ir_id)
  )
) |>
  mutate(
    match_group = match_group |>
      fct_reorder(-retained_introns) |>
      fct_relevel("Multiple modalities", "None", after = Inf)
  ) |>
  complete(match_group, fill = list(retained_introns = 0))

phenotype_totals <- read_tsv(
  "data/phenos/geuvadis-residual/geuvadis-residual-Geuvadis-cross_latent.phenotype_groups.txt.gz",
  col_names = c("phenotype_id", "gene_id"),
  show_col_types = FALSE
) |>
  mutate(modality = str_extract(phenotype_id, "^[^:]+")) |>
  count(modality, name = "total_phenotypes")

phenotype_matches <- ir_matches |>
  summarise(matched_genes = n_distinct(gene_id), .by = modality)

phenotype_counts <- phenotype_totals |>
  left_join(phenotype_matches, by = "modality") |>
  mutate(
    matched_genes = replace_na(matched_genes, 0L),
    modality_label = modalities[modality] |> fct_reorder(-matched_genes)
  )

p_ir <- match_breakdown |>
  mutate(match_group = fct_rev(match_group)) |>
  ggplot(aes(x = "All retained introns", y = retained_introns, fill = match_group)) +
  geom_col(width = 0.55, color = "black", linewidth = 0.25) +
  scale_fill_manual(values = match_group_colors) +
  scale_y_continuous(
    labels = label_number(big.mark = ","),
    expand = expansion(mult = c(0, 0.04))
  ) +
  theme_bw() +
  theme(
    axis.text = element_text(color = "black"),
    legend.key.size = unit(9, "pt"),
    legend.text = element_text(size = 7),
    legend.title = element_text(size = 9, margin = margin(b = 2, unit = "pt")),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  xlab(NULL) +
  ylab("MAJIQ retained introns") +
  labs(fill = "KDP/rDDP match")

p_phenotypes <- phenotype_counts |>
  ggplot(aes(x = modality_label, y = matched_genes, fill = modality_label)) +
  geom_col(
    width = 0.72,
    color = "black",
    linewidth = 0.25
  ) +
  scale_fill_manual(values = modality_colors, guide = "none") +
  scale_y_continuous(
    labels = label_number(big.mark = ","),
    expand = expansion(mult = c(0, 0.04))
  ) +
  theme_bw() +
  theme(
    axis.text = element_text(color = "black"),
    axis.text.x = element_text(angle = 35, hjust = 1, vjust = 1),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  xlab(NULL) +
  ylab("Genes with IR-matched phenotypes")

p_ir + p_phenotypes +
  plot_layout(widths = c(1, 2)) +
  plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(face = "bold"))

ggsave("figures/figureS5.png", width = 6, height = 4, device = png)
