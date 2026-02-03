## Compute enrichment of functional categories in xVariants of fine-mapped xQTLs
## For rDDPs and KDP phenotypes. Also filter to tissue-phenotype pairs with
## colocalizing TWAS hits, for enrichment of true causal SNPs.

library(GenomicRanges)
library(tidyverse)

load_finemap <- function(path) {
  read_tsv(
    path,
    col_types = cols(
      tissue = col_character(),
      phenotype_id = col_character(),
      variant_id = col_character(),
      pip = col_double(),
      af = col_double(),
      cs_id = col_integer()
    )
  ) |>
    separate_wider_delim(
      phenotype_id,
      ":",
      names = c("modality", "phenotype_id"),
      too_many = "merge"
    ) |>
    mutate(gene_id = str_replace(phenotype_id, "__.+$", ""))
}

snps_in_windows <- function(genes, snps, snps_rng) {
  cis_rng <- with(genes, GRanges(chrom, IRanges(TSS - 1e6, TSS + 1e6)))
  snps[countOverlaps(snps_rng, cis_rng) > 0]
}

background_snps <- function(gene_ids, genes, snps, snps_rng) {
  genes |>
    filter(gene_id %in% gene_ids) |>
    snps_in_windows(snps, snps_rng)
}

count_in_snps <- function(snps, anno) {
  map_int(anno, \(anno_snps) sum(snps %in% anno_snps)) |>
    enframe(name = "category", value = "count")
}

count_in_snps_weighted <- function(variants, anno) {
  weights <- set_names(variants$weight, variants$variant_id)
  map_dbl(anno, \(anno_snps) sum(weights[anno_snps], na.rm = TRUE)) |>
    enframe(name = "category", value = "count")
}

background_counts <- function(tested, genes, anno, anno_snps) {
  snp_info <- tibble(SNP = anno_snps) |>
    mutate(
      chrom = str_extract(SNP, "^(chr[^_]+)_", group = 1),
      pos = str_extract(SNP, "^[^_]+_([^_]+)_", group = 1) |> as.integer()
    )
  snps_rng <- with(snp_info, GRanges(chrom, IRanges(pos, pos)))
  tested |>
    reframe({
      bg_snps <- background_snps(gene_id, genes, snp_info$SNP, snps_rng)
      count_in_snps(bg_snps, anno) |>
        rename(count_bg = count) |>
        mutate(total_bg = length(bg_snps))
    }, .by = modality)
}

enrichment_table_binary <- function(qtls, anno, count_bg) {
  count_qtls <- qtls |>
    distinct(modality, variant_id) |>
    reframe(
      count_in_snps(variant_id, anno) |>
        rename(count_qtls = count) |>
        mutate(total_qtls = length(variant_id)),
      .by = modality
    )

  count_bg |>
    left_join(count_qtls,
              by = c("modality", "category")) |>
    mutate(
      count_qtls = replace_na(count_qtls, 0L),
      total_qtls = replace_na(total_qtls, 0L),
      frac_qtls = if_else(total_qtls > 0, count_qtls / total_qtls, NA_real_),
      frac_bg = count_bg / total_bg,
      log2_enrich = log2(frac_qtls / frac_bg)
    )
}

enrichment_table_weighted <- function(qtls, anno, count_bg) {
  count_qtls <- qtls |>
    reframe(
      count_in_snps_weighted(pick(variant_id, weight), anno) |>
        rename(count_qtls = count) |>
        mutate(total_qtls = sum(weight)),
      .by = modality
    )

  count_bg |>
    left_join(count_qtls,
              by = c("modality", "category")) |>
    mutate(
      count_qtls = replace_na(count_qtls, 0),
      total_qtls = replace_na(total_qtls, 0),
      frac_qtls = if_else(total_qtls > 0, count_qtls / total_qtls, NA_real_),
      frac_bg = count_bg / total_bg,
      log2_enrich = log2(frac_qtls / frac_bg)
    )
}

finemap <- load_finemap("data/finemap/cis_susie.all_tissues.tsv.gz")

tested_genes <- read_tsv(
  "data/phenos/gtex-residual/genes.tsv",
  col_types = cols(
    gene_id = col_character(),
    seqname = col_character(),
    window_start = col_integer(),
    window_end = col_integer(),
    strand = col_character(),
    tss = col_integer(),
    batch = col_integer()
  )
)

tested <- finemap |>
  distinct(modality) |>
  crossing(tested_genes |>
             distinct(gene_id))

genes <- tested_genes |>
  transmute(gene_id, chrom = seqname, TSS = tss)

anno <- read_tsv(
  "data/pantry/published_data/gtex/WGS_Feature_overlap_collapsed_VEP_short_4torus.MAF01.txt.gz",
  col_types = cols(SNP = "c", .default = "l")
)
anno_snps <- anno$SNP
anno <- colnames(anno)[-1] |>
  set_names() |>
  map(\(category) anno$SNP[anno[[category]]])

count_bg <- background_counts(tested, genes, anno, anno_snps)

## All xQTLs

variant_pips <- finemap |>
  group_by(modality, variant_id) |>
  summarize(weight = max(pip), .groups = "drop")

enrich_weighted <- enrichment_table_weighted(variant_pips, anno, count_bg) |>
  mutate(method = "pip_weighted")

enrich_binary <- variant_pips |>
  filter(weight > 0.5) |>
  enrichment_table_binary(anno = anno, count_bg = count_bg) |>
  mutate(method = "pip_gt_0.5")

enrich_all <- bind_rows(enrich_weighted, enrich_binary) |>
  select(method, modality, category, everything())

write_tsv(enrich_all, "data/analyses/enrich.finemap.tsv")

## xQTLs for tissue-phenotype pairs with colocalizing TWAS hits

coloc <- bind_rows(
  read_tsv("data/processed/gtex-residual.twas_hits.tsv.gz", col_types = "cc-c---d") |>
    mutate(modality = "latent_residual", .after = "tissue"),
  read_tsv("data/processed/gtex-pantry.twas_hits.tsv.gz", col_types = "ccc-c---d")
) |>
  filter(coloc_pp > 0.8) |>
  distinct(tissue, modality, phenotype_id)

variant_pips <- finemap |>
  semi_join(coloc, by = c("tissue", "modality", "phenotype_id"))
  group_by(modality, variant_id) |>
  summarize(weight = max(pip), .groups = "drop")

enrich_weighted <- enrichment_table_weighted(variant_pips, anno, count_bg) |>
  mutate(method = "pip_weighted")

enrich_binary <- variant_pips |>
  filter(weight > 0.5) |>
  enrichment_table_binary(anno = anno, count_bg = count_bg) |>
  mutate(method = "pip_gt_0.5")

enrich_all <- bind_rows(enrich_weighted, enrich_binary) |>
  select(method, modality, category, everything())

write_tsv(enrich_all, "data/analyses/enrich.finemap-coloc.tsv")
