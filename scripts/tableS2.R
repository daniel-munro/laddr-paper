# Number of DP, rDP, and KP TWAS hits per tissue per trait

library(tidyverse)

modalities <- c(
  expression = "KP_expression",
  isoforms = "KP_isoforms",
  splicing = "KP_splicing",
  alt_TSS = "KP_alt_TSS",
  alt_polyA = "KP_alt_polyA",
  stability = "KP_stability",
  latent_residual = "rDP"
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

twas_geuv_dp <- read_tsv("data/processed/geuvadis-full.twas_hits.tsv.gz", col_types = "cccdddd")
twas_geuv_rdp <- read_tsv("data/processed/geuvadis-residual.twas_hits.tsv.gz", col_types = "cccdddd")
twas_geuv_kp <- read_tsv("data/processed/geuvadis-pantry.twas_hits.tsv.gz", col_types = "ccccdddd")

n_dp <- twas_geuv_dp |>
  summarise(
    genes_DP = n_distinct(gene_id),
    hits_DP = n(),
    .by = trait
  )

n_rdp <- twas_geuv_rdp |>
  summarise(
    genes_rDP = n_distinct(gene_id),
    hits_rDP = n(),
    .by = trait
  )

n_kp <- twas_geuv_kp |>
  summarise(
    genes_KP = n_distinct(gene_id),
    hits_KP = n(),
    .by = trait
  )

n_kp_andor_rdp <- bind_rows(
  twas_geuv_kp |> mutate(type = "KP"),
  twas_geuv_rdp |> mutate(type = "rDP")
) |>
  distinct(trait, gene_id, type) |>
  summarise(
    types = str_c(sort(type), collapse = "_"),
    .by = c(trait, gene_id)
  ) |>
  count(trait, types) |>
  complete(trait, types, fill = list(n = 0)) |>
  pivot_wider(names_from = types, values_from = n) |>
  rename(genes_KP_no_rDP = KP,
         genes_KP_and_rDP = KP_rDP,
         genes_rDP_no_KP = rDP) |>
  mutate(genes_KP_or_rDP = genes_KP_no_rDP + genes_KP_and_rDP + genes_rDP_no_KP)

counts <- traits |>
  full_join(n_dp, by = "trait", relationship = "one-to-one") |>
  full_join(n_rdp, by = "trait", relationship = "one-to-one") |>
  full_join(n_kp, by = "trait", relationship = "one-to-one") |>
  full_join(n_kp_andor_rdp, by = "trait", relationship = "one-to-one") |>
  replace_na(
    list(genes_DP = 0, hits_DP = 0, genes_rDP = 0, hits_rDP = 0, genes_KP = 0,
         hits_KP = 0, genes_KP_no_rDP = 0, genes_KP_and_rDP = 0,
         genes_rDP_no_KP = 0, genes_KP_or_rDP = 0)
  ) |>
  select(trait, trait_name, hits_DP, hits_rDP, hits_KP, genes_DP, genes_rDP,
         genes_KP, genes_KP_or_rDP, genes_KP_no_rDP, genes_KP_and_rDP,
         genes_rDP_no_KP) |>
  arrange(desc(genes_KP_or_rDP), trait)

write_tsv(counts, "tables/tableS2.tsv")
