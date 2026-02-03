## Functional enrichment of finemapped hybrid xQTLs

library(tidyverse)

categories <- c(
  enhancer_d = "Enhancer",
  promoter_d = "Promoter",
  open_chromatin_region_d = "Open chromatin",
  promoter_flanking_region_d = "Promoter-flanking",
  CTCF_binding_site_d = "CTCF binding site",
  TF_binding_site_d = "TF binding site",
  `3_prime_UTR_variant_d` = "3' UTR",
  `5_prime_UTR_variant_d` = "5' UTR",
  frameshift_variant_d = "Frameshift",
  intron_variant_d = "Intron",
  missense_variant_d = "Missense",
  non_coding_transcript_exon_variant_d = "NC transcript",
  splice_acceptor_variant_d = "Splice acceptor",
  splice_donor_variant_d = "Splice donor",
  splice_region_variant_d = "Splice region",
  stop_gained_d = "Stop gained",
  synonymous_variant_d = "Synonymous",
  splicing = "Splicing",
  truncating = "Truncating"
)

modalities <- c(
  expression = "Expression",
  isoforms = "Isoform ratio",
  splicing = "Intron excision",
  alt_TSS = "Alt. TSS",
  alt_polyA = "Alt. polyA",
  stability = "RNA stability",
  latent_residual = "Residual DD"
)

modality_colors <- c(
  Expression = "#bf4042",
  `Isoform ratio` = "#6a90cd",
  `Intron excision` = "#59a257",
  `Alt. TSS` = "#896090",
  `Alt. polyA` = "#d97f26",
  `RNA stability` = "#ddb23c",
  `Residual DD` = "#1ce6df"
)

enrich <- read_tsv("data/analyses/enrich.finemap.tsv", col_types = "ccciiddddd") |>
  filter(method == "pip_weighted") |>
  mutate(modality = factor(modalities[modality], levels = modalities),
         category = factor(categories[category]) |>
           fct_reorder(if_else(modality == "Expression", log2_enrich, NA_real_), .na_rm = TRUE)
  )

enrich_coloc <- read_tsv("data/analyses/enrich.finemap-coloc.tsv", col_types = "ccciiddddd") |>
  filter(method == "pip_weighted") |>
  mutate(modality = factor(modalities[modality], levels = modalities),
         category = factor(categories[category]) |>
           fct_reorder(if_else(modality == "Expression", log2_enrich, NA_real_), .na_rm = TRUE)
  )

stripes <- tibble(y = seq(1, length(levels(enrich$category)), by = 2) - 0.5)

#############
## Panel b ## All xQTLs
#############

enrich |>
  mutate(category = as.integer(category) + (1/4) - (1/20) * as.integer(modality)) |>
  ggplot(aes(x = log2_enrich,
             y = category, color = modality, shape = modality)) +
  geom_rect(aes(ymin = y, ymax = y + 1, xmin = -Inf, xmax = Inf, color = NULL,
                shape = NULL, x = NULL, y = NULL),
            fill = "#eeeeee",
            data = stripes, show.legend = FALSE) +
  geom_vline(xintercept = 0, lty = 2) +
  geom_point(size = 1.5, stroke = 0.75, fill = modality_colors["Residual DD"]) +
  scale_color_manual(values = modality_colors) +
  scale_shape_manual(values = c(16, 17, 15, 4, 8, 5, 23)) +
  scale_y_continuous(breaks = 1:length(levels(enrich$category)),
                     labels = levels(enrich$category),
                     expand = c(0, 0)) +
  coord_cartesian(ylim = c(1 - 0.6, length(levels(enrich$category)) + 0.6)) +
  theme_classic() +
  theme(
    axis.text = element_text(color = "black"),
    legend.background = element_rect(color = "black", linewidth = 0.25),
    legend.justification = c(1, 0),
    legend.key.size = unit(10, "pt"),
    legend.position = "inside",
    legend.position.inside = c(0.98, 0.05),
    # plot.margin = unit(c(5.5, 1, 5.5, 5.5), "pt"),
  ) +
  xlab(expression(log[2]*" fold enrichment in xVariants")) +
  ylab(NULL) +
  labs(color = "Modality", shape = "Modality")

ggsave("figures/figure5/figure5b.png", width = 4.5, height = 3.5, device = png)

#############
## Panel c ## xQTLs for tissue-phenotype pairs with colocalizing TWAS hits
#############

enrich_coloc |>
  mutate(category = as.integer(category) + (1/4) - (1/20) * as.integer(modality)) |>
  ggplot(aes(x = log2_enrich,
             y = category, color = modality, shape = modality)) +
  geom_rect(aes(ymin = y, ymax = y + 1, xmin = -Inf, xmax = Inf, color = NULL,
                shape = NULL, x = NULL, y = NULL),
            fill = "#eeeeee",
            data = stripes, show.legend = FALSE) +
  geom_vline(xintercept = 0, lty = 2) +
  geom_point(size = 1.5, stroke = 0.75, fill = modality_colors["Residual DD"]) +
  scale_color_manual(values = modality_colors) +
  scale_shape_manual(values = c(16, 17, 15, 4, 8, 5, 23)) +
  scale_y_continuous(breaks = 1:length(levels(enrich_coloc$category)),
                     labels = levels(enrich_coloc$category),
                     expand = c(0, 0)) +
  coord_cartesian(ylim = c(1 - 0.6, length(levels(enrich_coloc$category)) + 0.6)) +
  theme_classic() +
  theme(
    axis.text = element_text(color = "black"),
    legend.background = element_rect(color = "black", linewidth = 0.25),
    legend.justification = c(1, 0),
    legend.key.size = unit(10, "pt"),
    legend.position = "inside",
    legend.position.inside = c(0.98, 0.05),
    plot.margin = unit(c(5.5, 1, 5.5, 5.5), "pt"),
  ) +
  xlab(expression(log[2]*" fold enrichment in xVariants")) +
  ylab(NULL) +
  labs(color = "Modality", shape = "Modality")

ggsave("figures/figure5/figure5c.png", width = 4.5, height = 3.5, device = png)
