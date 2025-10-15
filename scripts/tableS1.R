# Number of DDP, hybrid (KDP + rDDP), and KDP xQTLs per tissue

library(tidyverse)

modalities <- c(
  expression = "expression",
  isoforms = "isoforms",
  splicing = "splicing",
  alt_TSS = "alt_TSS",
  alt_polyA = "alt_polyA",
  stability = "stability",
  latent_residual = "rDDP",
  latent_full = "DDP"
)

tissue_info <- read_tsv(
  "data/info/tissueInfo.tsv",
  col_types = cols(tissueSiteDetailId = "c", tissueSiteDetailAbbr = "c", hasEGenes = "l", .default = "-")
) |>
  filter(hasEGenes) |>
  select(tissue = tissueSiteDetailAbbr,
         tissue_name = tissueSiteDetailId)

qtls_ddp <- read_tsv("data/processed/gtextcga-full.qtls.tsv.gz", col_types = "cciccd")

qtls_hybrid <- read_tsv("data/processed/gtex-residual-cross.qtls.tsv.gz", col_types = "ccicccd") |>
  mutate(modality = factor(modalities[modality], levels = modalities[1:7]))

qtls_kdp <- read_tsv("data/processed/gtex-pantry.qtls.tsv.gz", col_types = "ccicccd") |>
  mutate(modality = factor(modalities[modality], levels = modalities[1:6]))

qtls_geuv <- read_tsv("data/processed/geuvadis.qtls.tsv.gz", col_types = "ccccdci") |>
  mutate(modality = factor(modalities[modality], levels = modalities))

n_ddp <- qtls_ddp |>
  count(tissue, name = "QTLs_DDP")

n_hybrid <- qtls_hybrid |>
  count(tissue, modality) |>
  pivot_wider(id_cols = tissue, names_from = modality, values_from = n,
              names_prefix = "QTLs_hybrid_")

n_kdp <- qtls_kdp |>
  count(tissue, modality) |>
  pivot_wider(id_cols = tissue, names_from = modality, values_from = n,
              names_prefix = "QTLs_KDP_")

n_geuv <- qtls_geuv |>
  count(version, modality) |>
  mutate(
    prefix = c(
      `full-latent` = "QTLs",
      `residual-cross_latent` = "QTLs_hybrid",
      `residual-cross_pantry` = "QTLs_KDP"
    )[version],
    name = str_glue("{prefix}_{modality}"),
    tissue = "GEUVADIS"
  ) |>
  select(tissue, name, n) |>
  pivot_wider(id_cols = tissue, names_from = name, values_from = n) |>
  mutate(tissue_name = "Geuvadis_LCL", .after = tissue)

counts <- tissue_info |>
  full_join(n_ddp, by = "tissue", relationship = "one-to-one") |>
  full_join(n_hybrid, by = "tissue", relationship = "one-to-one") |>
  full_join(n_kdp, by = "tissue", relationship = "one-to-one") |>
  bind_rows(n_geuv) |>
  arrange(tissue)

write_tsv(counts, "tables/tableS1.tsv")
