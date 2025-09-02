# Number of DP, hybrid (KP + rDP), and KP xQTLs per tissue

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

tissue_info <- read_tsv(
  "data/info/tissueInfo.tsv",
  col_types = cols(tissueSiteDetailId = "c", tissueSiteDetailAbbr = "c", hasEGenes = "l", .default = "-")
) |>
  filter(hasEGenes) |>
  select(tissue = tissueSiteDetailAbbr,
         tissue_name = tissueSiteDetailId)

qtls_dp <- read_tsv("data/processed/gtextcga-full.qtls.tsv.gz", col_types = "cciccd")

qtls_hybrid <- read_tsv("data/processed/gtex-residual-cross.qtls.tsv.gz", col_types = "ccicccd") |>
  mutate(modality = factor(modalities[modality], levels = modalities))

qtls_kp <- read_tsv("data/processed/gtex-pantry.qtls.tsv.gz", col_types = "ccicccd") |>
  mutate(modality = factor(modalities[modality], levels = modalities[1:6]))

# qtls_geuv <- read_tsv("data/processed/geuvadis.qtls.tsv.gz", col_types = "ccccdci")

n_dp <- qtls_dp |>
  count(tissue, name = "QTLs_DP")

n_hybrid <- qtls_hybrid |>
  count(tissue, modality) |>
  pivot_wider(id_cols = tissue, names_from = modality, values_from = n,
              names_prefix = "QTLs_hybrid_")

n_kp <- qtls_kp |>
  count(tissue, modality) |>
  pivot_wider(id_cols = tissue, names_from = modality, values_from = n,
              names_prefix = "QTLs_")

# n_geuv <- qtls_geuv |>
#   count(version, modality) |>
#   mutate(
#     version = c(
#       `full-latent` = "DP",
#       `residual-cross_latent` = "hybrid",
#       `residual-cross_pantry` = "KP"
#     )[version],
#     name = str_glue("QTLs_{version}_{modality}"),
#     tissue = "Geuvadis"
#   ) |>
#   select(name, n) |>
#   pivot_wider(tissue, names_from = name, values_from = n)

counts <- tissue_info |>
  full_join(n_dp, by = "tissue", relationship = "one-to-one") |>
  full_join(n_hybrid, by = "tissue", relationship = "one-to-one") |>
  full_join(n_kp, by = "tissue", relationship = "one-to-one")

write_tsv(counts, "tables/tableS1.tsv")
