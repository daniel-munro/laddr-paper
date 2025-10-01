# Number of DDP, rDDP, and KDP TWAS hits per tissue per trait

library(tidyverse)

modalities <- c(
  expression = "KDP_expression",
  isoforms = "KDP_isoforms",
  splicing = "KDP_splicing",
  alt_TSS = "KDP_alt_TSS",
  alt_polyA = "KDP_alt_polyA",
  stability = "KDP_stability",
  latent_residual = "rDDP"
)

traits <- read_tsv("data/pantry/geuvadis/twas/gwas_metadata.txt",
                   col_types = cols(Tag = "c", Phenotype = "c", .default = "-")) |>
  # deframe()
  rename(trait = Tag, trait_name = Phenotype)

tissue_info <- read_tsv(
  "data/info/tissueInfo.tsv",
  col_types = cols(tissueSiteDetailId = "c", tissueSiteDetailAbbr = "c", hasEGenes = "l", .default = "-")
) |>
  filter(hasEGenes) |>
  select(tissue = tissueSiteDetailAbbr,
         tissue_name = tissueSiteDetailId)

twas_geuv_ddp <- read_tsv("data/processed/geuvadis-full.twas_hits.tsv.gz", col_types = "cccdddd")
twas_geuv_rddp <- read_tsv("data/processed/geuvadis-residual.twas_hits.tsv.gz", col_types = "cccdddd")
twas_geuv_kdp <- read_tsv("data/processed/geuvadis-pantry.twas_hits.tsv.gz", col_types = "ccccdddd")

n_ddp <- twas_geuv_ddp |>
  summarise(
    genes_DDP = n_distinct(gene_id),
    hits_DDP = n(),
    .by = trait
  )

n_rddp <- twas_geuv_rddp |>
  summarise(
    genes_rDDP = n_distinct(gene_id),
    hits_rDDP = n(),
    .by = trait
  )

n_kdp <- twas_geuv_kdp |>
  summarise(
    genes_KDP = n_distinct(gene_id),
    hits_KDP = n(),
    .by = trait
  )

n_kdp_andor_rddp <- bind_rows(
  twas_geuv_kdp |> mutate(type = "KDP"),
  twas_geuv_rddp |> mutate(type = "rDDP")
) |>
  distinct(trait, gene_id, type) |>
  summarise(
    types = str_c(sort(type), collapse = "_"),
    .by = c(trait, gene_id)
  ) |>
  count(trait, types) |>
  complete(trait, types, fill = list(n = 0)) |>
  pivot_wider(names_from = types, values_from = n) |>
  rename(genes_KDP_no_rDDP = KDP,
         genes_KDP_and_rDDP = KDP_rDDP,
         genes_rDDP_no_KDP = rDDP) |>
  mutate(genes_KDP_or_rDDP = genes_KDP_no_rDDP + genes_KDP_and_rDDP + genes_rDDP_no_KDP)

counts <- traits |>
  full_join(n_ddp, by = "trait", relationship = "one-to-one") |>
  full_join(n_rddp, by = "trait", relationship = "one-to-one") |>
  full_join(n_kdp, by = "trait", relationship = "one-to-one") |>
  full_join(n_kdp_andor_rddp, by = "trait", relationship = "one-to-one") |>
  replace_na(
    list(genes_DDP = 0, hits_DDP = 0, genes_rDDP = 0, hits_rDDP = 0, genes_KDP = 0,
         hits_KDP = 0, genes_KDP_no_rDDP = 0, genes_KDP_and_rDDP = 0,
         genes_rDDP_no_KDP = 0, genes_KDP_or_rDDP = 0)
  ) |>
  select(trait, trait_name, hits_DDP, hits_rDDP, hits_KDP, genes_DDP, genes_rDDP,
         genes_KDP, genes_KDP_or_rDDP, genes_KDP_no_rDDP, genes_KDP_and_rDDP,
         genes_rDDP_no_KDP) |>
  arrange(desc(genes_KDP_or_rDDP), trait)

write_tsv(counts, "tables/tableS2.tsv")
