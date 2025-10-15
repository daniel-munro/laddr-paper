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
  rename(trait = Tag, trait_name = Phenotype)

tissue_info <- read_tsv(
  "data/info/tissueInfo.tsv",
  col_types = cols(tissueSiteDetailId = "c", tissueSiteDetailAbbr = "c", hasEGenes = "l", .default = "-")
) |>
  filter(hasEGenes) |>
  select(tissue = tissueSiteDetailAbbr,
         tissue_name = tissueSiteDetailId) |>
  bind_rows(
    tibble(tissue = "GEUVADIS",
           tissue_name = "Geuvadis_LCL")
  )

twas_gtex_ddp <- read_tsv("data/processed/gtextcga-full.twas_hits.tsv.gz", col_types = "ccccdddd")
twas_gtex_rddp <- read_tsv("data/processed/gtex-residual.twas_hits.tsv.gz", col_types = "ccccdddd")
twas_gtex_kdp <- read_tsv("data/processed/gtex-pantry.twas_hits.tsv.gz", col_types = "cccccdddd")

twas_geuv_ddp <- read_tsv("data/processed/geuvadis-full.twas_hits.tsv.gz", col_types = "cccdddd")
twas_geuv_rddp <- read_tsv("data/processed/geuvadis-residual.twas_hits.tsv.gz", col_types = "cccdddd")
twas_geuv_kdp <- read_tsv("data/processed/geuvadis-pantry.twas_hits.tsv.gz", col_types = "ccccdddd")

twas_ddp <- bind_rows(
  twas_gtex_ddp,
  twas_geuv_ddp |> mutate(tissue = "GEUVADIS"),
)

twas_rddp <- bind_rows(
  twas_gtex_rddp,
  twas_geuv_rddp |> mutate(tissue = "GEUVADIS"),
)

twas_kdp <- bind_rows(
  twas_gtex_kdp,
  twas_geuv_kdp |> mutate(tissue = "GEUVADIS"),
)

n_ddp <- twas_ddp |>
  summarise(
    genes_DDP = n_distinct(gene_id),
    hits_DDP = n(),
    .by = c(trait, tissue)
  )

n_rddp <- twas_rddp |>
  summarise(
    genes_rDDP = n_distinct(gene_id),
    hits_rDDP = n(),
    .by = c(trait, tissue)
  )

n_kdp <- twas_kdp |>
  summarise(
    genes_KDP = n_distinct(gene_id),
    hits_KDP = n(),
    .by = c(trait, tissue)
  )

n_kdp_andor_rddp <- bind_rows(
  twas_kdp |> mutate(type = "KDP"),
  twas_rddp |> mutate(type = "rDDP")
) |>
  distinct(trait, tissue, gene_id, type) |>
  summarise(
    types = str_c(sort(type), collapse = "_"),
    .by = c(trait, tissue, gene_id)
  ) |>
  count(trait, tissue, types) |>
  complete(trait, tissue, types, fill = list(n = 0)) |>
  pivot_wider(names_from = types, values_from = n) |>
  rename(genes_KDP_no_rDDP = KDP,
         genes_KDP_and_rDDP = KDP_rDDP,
         genes_rDDP_no_KDP = rDDP) |>
  mutate(genes_KDP_or_rDDP = genes_KDP_no_rDDP + genes_KDP_and_rDDP + genes_rDDP_no_KDP)

counts <- crossing(trait = traits$trait,
                   tissue = tissue_info$tissue) |>
  full_join(n_ddp, by = c("trait", "tissue"), relationship = "one-to-one") |>
  full_join(n_rddp, by = c("trait", "tissue"), relationship = "one-to-one") |>
  full_join(n_kdp, by = c("trait", "tissue"), relationship = "one-to-one") |>
  full_join(n_kdp_andor_rddp, by = c("trait", "tissue"), relationship = "one-to-one") |>
  left_join(traits, by = "trait", relationship = "many-to-one") |>
  left_join(tissue_info, by = "tissue", relationship = "many-to-one") |>
  replace_na(
    list(genes_DDP = 0, hits_DDP = 0, genes_rDDP = 0, hits_rDDP = 0, genes_KDP = 0,
         hits_KDP = 0, genes_KDP_no_rDDP = 0, genes_KDP_and_rDDP = 0,
         genes_rDDP_no_KDP = 0, genes_KDP_or_rDDP = 0)
  ) |>
  select(trait, trait_name, tissue, tissue_name, hits_DDP, hits_rDDP, hits_KDP,
         genes_DDP, genes_rDDP, genes_KDP, genes_KDP_or_rDDP, genes_KDP_no_rDDP,
         genes_KDP_and_rDDP, genes_rDDP_no_KDP) |>
  arrange(desc(genes_KDP_or_rDDP), trait, tissue)

write_tsv(counts, "tables/tableS2.tsv")
