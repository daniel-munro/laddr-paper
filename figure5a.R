## Gene-relative positions of KP, DP, and RP xQTLs

library(tidyverse)

modality_colors <- c(
  latent_full = "#13918d",
  latent_residual = "#1ce6df",
  expression = "#bf4042",
  isoforms = "#6a90cd",
  splicing = "#59a257",
  alt_TSS = "#896090",
  alt_polyA = "#d97f26",
  stability = "#ddb23c"
)

modality_labels <- c(
  latent_full = "Data-driven",
  latent_residual = "Residual data-driven",
  expression = "Expression",
  isoforms = "Isoform ratio",
  splicing = "Intron excision",
  alt_TSS = "Alternative TSS",
  alt_polyA = "Alternative polyA",
  stability = "RNA stability"
) |>
  enframe(name = "modality", value = "label") |>
  mutate(label = factor(label, levels = label))

genes <- read_tsv("data/processed/pcg_and_lncrna.tsv", col_types = "c-cciic") |>
  mutate(TSS = if_else(strand == "-", end, start),
         TES = if_else(strand == "-", start, end)) |>
  select(gene_id, chrom, TSS, TES)

qtls <- read_tsv("data/processed/geuvadis.qtls.tsv.gz", col_types = "ccccdci") |>
  filter(version %in% c("full-latent", "residual-cross_latent"))

qtls_pos <- qtls |>
  mutate(chrom = str_split_i(variant_id, "_", 1),
         pos = str_split_i(variant_id, "_", 2) |> as.integer()) |>
  left_join(genes, by = "gene_id") |>
  mutate(
    rel_pos_gene = (pos - TSS) / (TES - TSS),
  ) |>
  left_join(modality_labels, by = "modality") |>
  filter(rel_pos_gene >= -1,
         rel_pos_gene <= 2)
stopifnot(all(qtls_pos$chrom.x == qtls_pos$chrom.y))

ggplot(qtls_pos, aes(x = rel_pos_gene, fill = modality)) +
  # facet_grid(rows = vars(modality), scales = "free_y", drop = TRUE) +
  facet_wrap(~label, scales = "free_y", drop = TRUE, ncol = 1) +
  geom_histogram(bins = 100, show.legend = FALSE) +
  geom_vline(xintercept = c(0, 1), alpha = 1, linewidth = 0.3) +
  # geom_text(mapping = aes(label = label), data = modality_labels,
  #           x = -0.95, y = Inf, hjust = 0, vjust = 1, fontface = 1) +
  scale_fill_manual(values = modality_colors) +
  scale_x_continuous(expand = c(0, 0), breaks = c(0, 1),
                     labels = c("Gene start", "Gene end")) +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 2)) +
  theme_classic() +
  theme(
    axis.text = element_text(color = "black"),
    # strip.text = element_blank(),
  ) +
  xlab("xVariant position normalized to xGene length") +
  ylab("No. xQTLs")

ggsave("figures/figure5/figure5a.png", width = 4.5, height = 7, device = png)
