# FinnGen colocalizations

library(tidyverse)

finngen_coloc <- read_tsv(
  "data/finngen_coloc/coloc.significant.tsv.gz",
  col_types = cols(tissue = "c", phenotype_id = "c", finngen_trait = "c", finngen_region = "c",
                   cs_overlap = "c", clpp = "d", clpa = "d",
                   top_shared_variant = "c", .default = "-")
) |>
  filter(clpp >= 0.01) |>
  separate_wider_delim(phenotype_id, delim = ":", names = c("modality", "phenotype_id")) |>
  mutate(
    gene_id = phenotype_id |>
      str_replace("^.*:", "") |>
      str_replace("__.*$", ""),
    .before = modality
  )

write_tsv(finngen_coloc, "tables/tableS3.tsv")
