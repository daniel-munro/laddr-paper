## Colocalization analysis between KDP/rDDP xQTLs and FinnGen GWAS traits

library(locuszoomr)
library(EnsDb.Hsapiens.v86)
library(tidyverse)
library(patchwork)
library(scales)

traits <- read_tsv("data/finngen_coloc/finngen_R12_manifest.tsv", col_types = "c-c----") |>
  mutate(
    category = category |>
      str_replace("^[XVI]+ ", "") |>
      str_replace(" \\(.+\\)$", "")
  )
  
coloc <- read_tsv(
  "data/finngen_coloc/coloc.significant.tsv.gz",
  col_types = cols(phenotype_id = "c", phenotype_class = "c", finngen_trait = "c", clpp = "d", .default = "-")
) |>
  filter(clpp >= 0.01) |>
  mutate(
    qtl_gene_id = phenotype_id |>
      str_replace("^.*:", "") |>
      str_replace("__.*$", ""),
    modality_group = if_else(phenotype_class == "latent_residual", "rDDP", "KDP"),
  ) |>
  left_join(traits, by = c("finngen_trait" = "phenocode"))

pairs <- coloc |>
  summarise(
    modality_groups = str_c(sort(unique(modality_group)), collapse = "_"),
    .by = c(qtl_gene_id, finngen_trait, category)
  ) |>
  mutate(
    category = category |>
      fct_infreq() |>
      fct_lump_n(n = 25, other_level = "Other categories"),
    modality_groups = c(
      KDP = "KDP only",
      KDP_rDDP = "KDP & rDDP",
      rDDP = "rDDP only"
    )[modality_groups] |>
      fct_relevel(c("KDP only", "KDP & rDDP", "rDDP only"))
  )

# Panel a: bar plot of gene-trait pairs, for each trait category, with rDDP-only, rDDP + KDP, and KDP-only colocs (CLPP > 0.1?)

p_barplot <- pairs |>
  ggplot(aes(fill = modality_groups, y = category)) +
  geom_bar(width = 0.8, color = "black") +
  scale_x_continuous(expand = expansion(mult = c(0, 0.01))) +
  scale_fill_manual(values = c("white", "gray", "#444444")) +
  theme_classic() +
  theme(
    legend.position = "inside",
    legend.position.inside = c(0.9, 0.9),
    legend.justification.inside = c(1, 1),
  ) +
  labs(fill = "xQTL types") +
  xlab("Colocalized gene-trait pairs") +
  ylab("FinnGen trait category")

# Panels b-d: locuszoom-style plots for selected example colocalizations

example_ids_to_plot <- str_glue("example_{c(3, 4, 5)}")
plot_flank <- 500000
ens_db <- "EnsDb.Hsapiens.v86"
ld_cache_file <- "data/finngen_coloc/coloc_examples.ld_cache.tsv"

stopifnot(file.exists(ld_cache_file))

examples <- read_tsv(
  "data/finngen_coloc/coloc_examples.annotated.tsv",
  show_col_types = FALSE
) |>
  mutate(
    example_id = factor(example_id, levels = unique(example_id)),
    gene_id = str_extract(phenotype_id, "ENSG[0-9]+"),
    pc = str_extract(phenotype_id, "PC[0-9]+$") |> str_replace("PC", "rDDP"),
    region_chr = str_extract(finngen_region, "^chr[^:]+"),
    region_chr_no_prefix = str_remove(region_chr, "^chr"),
    region_start = str_match(finngen_region, ":([0-9]+)-")[, 2] |> as.integer(),
    region_end = str_match(finngen_region, "-([0-9]+)$")[, 2] |> as.integer(),
    shared_pos = str_match(top_shared_variant, "^[^:]+:([0-9]+):")[, 2] |> as.integer(),
    plot_start = pmax(region_start, shared_pos - plot_flank),
    plot_end = pmin(region_end, shared_pos + plot_flank),
    shared_snp = as.character(str_glue("{region_chr}:{shared_pos}"))
  ) |>
  filter(example_id %in% example_ids_to_plot)

gene_names <- read_tsv("data/processed/pcg_and_lncrna.tsv", col_types = "cc-----") |>
  rename(gene_label = gene_name)

examples <- examples |>
  left_join(gene_names, by = "gene_id", relationship = "many-to-one") |>
  mutate(
    gene_label = coalesce(na_if(gene_label, "NA"), gene_id),
    example_label = str_glue("{gene_label} {pc} in {tissue}")
  )

ld_cache <- read_tsv(ld_cache_file, show_col_types = FALSE) |>
  mutate(
    example_id = as.character(example_id),
    snp = as.character(snp),
    ld = as.numeric(ld)
  ) |>
  filter(example_id %in% as.character(examples$example_id))

missing_ld <- setdiff(as.character(examples$example_id), unique(ld_cache$example_id))
if (length(missing_ld) > 0) {
  stop(
    "Missing LD cache for: ",
    str_c(missing_ld, collapse = ", "),
    ". Render analyses/finngen_coloc/coloc_locuszoom.qmd first."
  )
}

qtl_region <- read_tsv(
  "data/finngen_coloc/coloc_examples.qtl_region.tsv.gz",
  show_col_types = FALSE,
  guess_max = 200000
) |>
  mutate(
    example_id = as.character(example_id),
    chrom = str_remove(variant_chrom, "^chr"),
    pos = str_match(variant_key, "^[^:]+:([0-9]+):")[, 2] |> as.integer(),
    p = pmax(as.numeric(pval_nominal), .Machine$double.xmin),
    snp = as.character(str_glue("chr{chrom}:{pos}"))
  ) |>
  filter(example_id %in% as.character(examples$example_id)) |>
  dplyr::select(example_id, chrom, pos, p, snp, variant_key) |>
  left_join(ld_cache, by = c("example_id", "snp")) |>
  left_join(
    examples |>
      transmute(example_id = as.character(example_id), shared_snp),
    by = "example_id"
  ) |>
  mutate(ld = if_else(snp == shared_snp, 1, ld))

finngen_region <- read_tsv(
  "data/finngen_coloc/coloc_examples.finngen_region.tsv.gz",
  show_col_types = FALSE,
  guess_max = 200000
) |>
  mutate(
    example_id = as.character(example_id),
    chrom = str_remove(as.character(`#chrom`), "^chr"),
    pos = as.integer(pos),
    p = pmax(as.numeric(pval), .Machine$double.xmin),
    snp = as.character(str_glue("chr{chrom}:{pos}"))
  ) |>
  filter(example_id %in% as.character(examples$example_id)) |>
  dplyr::select(example_id, chrom, pos, p, snp, variant_key) |>
  left_join(ld_cache, by = c("example_id", "snp")) |>
  left_join(
    examples |>
      transmute(example_id = as.character(example_id), shared_snp),
    by = "example_id"
  ) |>
  mutate(ld = if_else(snp == shared_snp, 1, ld))

make_locus <- function(df, info) {
  df <- df |>
    filter(!is.na(pos), !is.na(p))

  locus(
    data = as.data.frame(df),
    xrange = c(info$plot_start, info$plot_end),
    seqname = info$region_chr_no_prefix,
    ens_db = ens_db,
    chrom = "chrom",
    pos = "pos",
    p = "p",
    labs = "snp",
    index_snp = info$shared_snp,
    LD = "ld"
  )
}

plot_locus_example <- function(example) {
  info <- examples |>
    filter(example_id == example) |>
    slice_head(n = 1)

  qtl_loc <- make_locus(
    qtl_region |> filter(example_id == info$example_id),
    info
  )

  fg_loc <- make_locus(
    finngen_region |> filter(example_id == info$example_id),
    info
  )

  p_qtl <- gg_scatter(
    qtl_loc,
    index_snp = info$shared_snp,
    pcutoff = 0,
    labels = NULL,
    xticks = FALSE,
    ylab = expression(xQTL~-log[10](italic(P))),
    showLD = TRUE,
    legend_pos = "topright",
    size = 2,
    LD_scheme = c("grey75", "royalblue", "cyan2",
                  "green3", "orange", "red", "purple")
  ) +
    scale_x_continuous(expand = 0) +
    labs(title = info$example_label) +
    theme(
      plot.title = element_text(size = 9, face = "bold"),
      legend.position = "none",
      axis.title.x = element_blank(),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
    )

  p_fg <- gg_scatter(
    fg_loc,
    index_snp = info$shared_snp,
    pcutoff = 5e-8,
    labels = NULL,
    xticks = FALSE,
    ylab = expression(FinnGen~-log[10](italic(P))),
    showLD = TRUE,
    legend_pos = "topright",
    size = 2,
    LD_scheme = c("grey75", "royalblue", "cyan2",
                  "green3", "orange", "red", "purple")
  ) +
    scale_x_continuous(expand = 0) +
    labs(title = info$finngen_trait) +
    theme(
      plot.title = element_text(size = 8),
      legend.position = "none",
      axis.title.x = element_blank(),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
    )

  p_genes <- gg_genetracks(
    fg_loc,
    xticks = TRUE,
    xlab = str_glue("{info$region_chr} position (Mb)"),
    showExons = TRUE,
    maxrows = 5,
    highlight = info$gene_label,
    highlight_col = "red3",
    blanks = "fill",
    gene_col = "grey35",
    exon_col = "grey35",
    exon_border = "grey35"
  ) +
    scale_x_continuous(expand = 0)

  p_qtl / p_fg / p_genes +
    plot_layout(heights = c(2.1, 2.1, 1.4))
}

locus_plots <- map(as.character(examples$example_id), plot_locus_example)
locus_panels <- map(locus_plots, wrap_elements)

combined_plot <- free(p_barplot) + locus_panels +
  plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(face = "bold"))

combined_plot

output_dir <- Sys.getenv("FIGURES_DIR", unset = "figures")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

ggsave(
  file.path(output_dir, "figureS7.png"),
  width = 13,
  height = 10,
  device = png
)
