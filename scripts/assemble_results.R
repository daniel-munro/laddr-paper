library(tidyverse)

tissues_gtex5 <- read_lines("data/info/tissues.gtex5.txt")
tissues_gtex <- read_lines("data/info/tissues.gtex.txt")

modalities <- c("expression", "isoforms", "splicing", "alt_TSS", "alt_polyA", "stability")

###########
## Genes ##
###########

genes <- rtracklayer::import("data/ref/Homo_sapiens.GRCh38.113.chr.gtf.gz") |>
    as_tibble() |>
    filter(type == "gene",
           gene_biotype %in% c("protein_coding", "lncRNA")) |>
    mutate(chrom = str_c("chr", seqnames)) |>
    select(gene_id, gene_name, gene_biotype, chrom, start, end, strand)

write_tsv(genes, "data/processed/pcg_and_lncrna.tsv")

##########
## QTLs ##
##########

qtls_geuvadis <- bind_rows(
  read_tsv(
    "data/qtl/geuvadis-residual/Geuvadis-cross_pantry.cis_independent_qtl.txt.gz",
    col_types = "c-----c---------cc-i"
  ) |>
    rename(gene_id = group_id) |>
    separate_wider_delim(
      phenotype_id, ":", names = c("modality", "phenotype_id"), too_many = "merge"
    ) |>
    mutate(version = "residual-cross_pantry", .before = 1),
  read_tsv(
    "data/qtl/geuvadis-residual/Geuvadis-cross_latent.cis_independent_qtl.txt.gz",
    col_types = "c-----c---------cc-i"
  ) |>
    rename(gene_id = group_id) |>
    separate_wider_delim(
      phenotype_id, ":", names = c("modality", "phenotype_id"), too_many = "merge"
    ) |>
    mutate(version = "residual-cross_latent", .before = 1),
  read_tsv(
    "data/qtl/geuvadis-full/Geuvadis-latent.cis_independent_qtl.txt.gz",
    col_types = "c-----c---------cc-i"
  ) |>
    rename(gene_id = group_id) |>
    mutate(modality = "latent_full", .before = phenotype_id) |>
    mutate(version = "full-latent", .before = 1),
) |>
  filter(gene_id %in% genes$gene_id)

write_tsv(qtls_geuvadis, "data/processed/geuvadis.qtls.tsv.gz")

# qtlassoc_gtex5_full <- tibble(tissue = tissues_gtex5) |>
#     reframe(
#         read_tsv(
#             str_glue("data/qtl/gtex5-full/{tissue}-latent.cis_qtl.txt.gz"),
#             col_types = "c-----c---------cc-c-"
#         ),
#         .by = tissue
#     ) |>
#     rename(gene_id = group_id) |>
#     relocate(gene_id, .before = 2)
# 
# write_tsv(qtlassoc_gtex5_full, "data/processed/gtex5-full.qtl_assoc.tsv.gz")

qtls_gtex5_full <- tibble(tissue = tissues_gtex) |>
    reframe(
        read_tsv(
            str_glue("data/qtl/gtex5-full/{tissue}-latent.cis_independent_qtl.txt.gz"),
            col_types = "c-----c---------cc-i"
        ),
        .by = tissue
    ) |>
    select(tissue, gene_id = group_id, rank, phenotype_id, variant_id, pval_beta)

write_tsv(qtls_gtex5_full, "data/processed/gtex5-full.qtls.tsv.gz")

qtls_gtex_full <- tibble(tissue = tissues_gtex) |>
    reframe(
        read_tsv(
            str_glue("data/qtl/gtex-full/{tissue}-latent.cis_independent_qtl.txt.gz"),
            col_types = "c-----c---------cc-i"
        ),
        .by = tissue
    ) |>
    select(tissue, gene_id = group_id, rank, phenotype_id, variant_id, pval_beta)

write_tsv(qtls_gtex_full, "data/processed/gtex-full.qtls.tsv.gz")

qtls_gtextcga_full <- tibble(tissue = tissues_gtex) |>
    reframe(
        read_tsv(
            str_glue("data/qtl/gtextcga-full/{tissue}-latent.cis_independent_qtl.txt.gz"),
            col_types = "c-----c---------cc-i"
        ),
        .by = tissue
    ) |>
    select(tissue, gene_id = group_id, rank, phenotype_id, variant_id, pval_beta)

write_tsv(qtls_gtextcga_full, "data/processed/gtextcga-full.qtls.tsv.gz")

qtls_gtex5_cross <- tibble(tissue = tissues_gtex5) |>
    reframe(
        read_tsv(
            str_glue("data/qtl/gtex5-residual/{tissue}-cross_latent.cis_independent_qtl.txt.gz"),
            col_types = "c-----c---------cc-i"
        ),
        .by = tissue
    ) |>
    select(tissue, gene_id = group_id, rank, phenotype_id, variant_id, pval_beta) |>
    separate_wider_delim(phenotype_id, names = c("modality", "phenotype_id"), delim = ":", too_many = "merge")

write_tsv(qtls_gtex5_cross, "data/processed/gtex5-residual-cross.qtls.tsv.gz")

qtls_gtex_cross <- tibble(tissue = tissues_gtex) |>
    reframe(
        read_tsv(
            str_glue("data/qtl/gtex-residual/{tissue}-cross_latent.cis_independent_qtl.txt.gz"),
            col_types = "c-----c---------cc-i"
        ),
        .by = tissue
    ) |>
    select(tissue, gene_id = group_id, rank, phenotype_id, variant_id, pval_beta) |>
    separate_wider_delim(phenotype_id, names = c("modality", "phenotype_id"), delim = ":", too_many = "merge")

write_tsv(qtls_gtex_cross, "data/processed/gtex-residual-cross.qtls.tsv.gz")

##################
## QTLs special ##
##################

qtls_prune <- crossing(
  map_group = c("latent", "pantry", modalities),
  pruning = c(0, 20, 40, 60, 80, 100)
) |>
  reframe(
    {
      fname <- str_glue("data/prune_anno/{map_group}-{pruning}.cis_independent_qtl.txt.gz")
      if (map_group %in% c("expression", "stability")) {
        df <- read_tsv(fname, col_types = "c-----c---------ci")
      } else {
        df <- read_tsv(fname, col_types = "c-----c---------cc-i")
      }
      if (map_group == "pantry") {
        separate_wider_delim(df, phenotype_id, ":", names = c("modality", "phenotype_id"),
                             too_many = "merge")
      } else {
        df |>
          mutate(modality = map_group, .before = phenotype_id)
      }
    },
    .by = c(map_group, pruning)
  ) |>
  mutate(gene_id = if_else(is.na(group_id), phenotype_id, group_id)) |>
  select(map_group, pruning, modality, gene_id, rank, phenotype_id, variant_id, pval_beta)

write_tsv(qtls_prune, "data/processed/prune-BRNCTXB.qtls.tsv.gz")

qtls_held_out <- bind_rows(
  qtls_geuvadis |>
    filter(version == "residual-cross_latent") |>
    select(-version) |>
    mutate(held_out = "none", .before = 1),
  tibble(held_out = modalities) |>
    reframe(
      read_tsv(
        str_glue("data/held_out/cross-no_{held_out}.cis_independent_qtl.txt.gz"),
        col_types = "c-----c---------cc-i"
      ) |>
        rename(gene_id = group_id) |>
        separate_wider_delim(
          phenotype_id, ":", names = c("modality", "phenotype_id"), too_many = "merge"
        ),
      .by = held_out
    )
)

write_tsv(qtls_held_out, "data/processed/held_out-geuvadis.qtls.tsv.gz")

qtls_seqsim <- tibble(
  reads = c("pe-50bp-100pct", "pe-75bp-100pct", "se-50bp-100pct", "se-75bp-100pct", "se1-50bp-100pct", "se1-75bp-100pct", "se2-50bp-100pct", "se2-75bp-100pct", "pe-75bp-50pct", "pe-75bp-25pct")
) |>
  reframe(
    read_tsv(
      str_glue("data/seqsim/{reads}.cis_independent_qtl.txt.gz"),
      col_types = "c-----c---------cc-i"
    ),
    .by = reads
  ) |>
  rename(gene_id = group_id)

write_tsv(qtls_seqsim, "data/processed/seqsim.qtls.tsv.gz")

##########
## TWAS ##
##########

twas_geuv_full <- read_tsv(
  str_glue("data/twas/twas_hits.geuvadis-full-Geuvadis.tsv"),
  col_types = "cccccccccccccccccccccccc"
) |>
  mutate(gene_id = str_split_i(ID, "__", i = 1)) |>
  select(trait = TRAIT, gene_id, phenotype_id = ID, hsq = HSQ, twas_z = TWAS.Z, twas_p = TWAS.P, coloc_pp = COLOC.PP4) |>
  arrange(trait, gene_id, phenotype_id)

write_tsv(twas_geuv_full, "data/processed/geuvadis-full.twas_hits.tsv.gz")

twas_gtextcga_full <- tibble(tissue = tissues_gtex5) |>
    reframe(
        read_tsv(
            str_glue("data/twas/twas_hits.gtextcga-full-{tissue}.tsv"),
            col_types = "cccccccccccccccccccccccc"
        ),
        .by = tissue
    ) |>
    mutate(gene_id = str_split_i(ID, "__", i = 1)) |>
    select(tissue, trait = TRAIT, gene_id, phenotype_id = ID, hsq = HSQ, twas_z = TWAS.Z, twas_p = TWAS.P, coloc_pp = COLOC.PP4) |>
    arrange(tissue, trait, gene_id, phenotype_id)

write_tsv(twas_gtextcga_full, "data/processed/gtextcga-full.twas_hits.tsv.gz")
