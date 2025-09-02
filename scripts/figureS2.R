## Gene-relative positions of held-out latent QTLs

library(tidyverse)

modalities <- c(
  expression = "Expression",
  isoforms = "Isoform ratio",
  splicing = "Intron excision",
  alt_TSS = "Alt. TSS",
  alt_polyA = "Alt. polyA",
  stability = "RNA stability"
)

modality_colors <- c(
  Expression = "#bf4042",
  `Isoform ratio` = "#6a90cd",
  `Intron excision` = "#59a257",
  `Alt. TSS` = "#896090",
  `Alt. polyA` = "#d97f26",
  `RNA stability` = "#ddb23c"
)

genes <- read_tsv("data/processed/pcg_and_lncrna.tsv", col_types = "c-cciic") |>
  mutate(TSS = if_else(strand == "-", end, start),
         TES = if_else(strand == "-", start, end)) |>
  select(gene_id, chrom, TSS, TES)

qtls_rp <- read_tsv(
  "data/processed/held_out-geuvadis.qtls.tsv.gz", col_types = "ccccdci"
) |>
  filter(modality == "latent_residual")

original_genes <- qtls_rp |>
  filter(held_out == "none") |>
  distinct(gene_id) |>
  pull()

qtls_rp_new_genes <- qtls_rp |>
  filter(held_out != "none",
         !(gene_id %in% original_genes)) |>
  mutate(held_out = c(modalities)[held_out] |> fct_inorder())

qtls_rp_pos <- qtls_rp_new_genes |>
  mutate(chrom = str_split_i(variant_id, "_", 1),
         pos = str_split_i(variant_id, "_", 2) |> as.integer()) |>
  left_join(genes, by = "gene_id") |>
  mutate(
    rel_pos_gene = (pos - TSS) / (TES - TSS),
    rel_pos_TSS = if_else(TSS < TES, pos - TSS, TSS - pos),
    rel_pos_TES = if_else(TSS < TES, pos - TES, TES - pos),
  )
stopifnot(all(qtls_rp_pos$chrom.x == qtls_rp_pos$chrom.y))

df <- qtls_rp_pos |>
  filter(rel_pos_gene >= -1,
         rel_pos_gene <= 2) |>
  mutate(modality = held_out)

ggplot(df, aes(x = rel_pos_gene, fill = modality)) +
  facet_grid(rows = vars(modality), scales = "free_y", drop = TRUE) +
  geom_histogram(bins = 80, show.legend = FALSE) +
  scale_fill_manual(values = modality_colors) +
  scale_x_continuous(expand = c(0, 0), breaks = c(0, 1),
                     labels = c("Gene start", "Gene end")) +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 2)) +
  geom_vline(xintercept = c(0, 1), alpha = 0.5) +
  geom_text(mapping = aes(label = modality), data = distinct(df, modality),
            x = -0.95, y = Inf, hjust = 0, vjust = 1, fontface = 1) +
  theme_classic() +
  theme(
    axis.text = element_text(color = "black"),
    strip.text = element_blank(),
  ) +
  xlab("xVariant position normalized to xGene length") +
  ylab("No. xQTLs") +
  ggtitle("New rDP xQTLs when holding out the indicated modality")

ggsave("figures/figureS2.png", width = 6, height = 6, device = png)
