library(tidyverse)

tissues_gtex5 <- read_lines("data/info/tissues.gtex5.txt")
tissues_gtex <- read_lines("data/info/tissues.gtex.txt")

###########
## Genes ##
###########

genes <- rtracklayer::import("data/ref/Homo_sapiens.GRCh38.113.chr.gtf.gz") |>
    as_tibble() |>
    filter(type == "gene",
           gene_biotype == "protein_coding") |>
    mutate(chrom = str_c("chr", seqnames)) |>
    select(gene_id, gene_name, chrom, start, end, strand)

write_tsv(genes, "data/processed/protein_coding_genes.tsv")

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
)

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

##########
## TWAS ##
##########

twas_geuvadis <- read_tsv(
  str_glue("data/twas/twas_hits.geuvadis-full-Geuvadis.tsv"),
  col_types = "cccccccccccccccccccccccc"
) |>
  mutate(gene_id = str_split_i(ID, "__", i = 1)) |>
  select(trait = TRAIT, gene_id, phenotype_id = ID, hsq = HSQ, twas_z = TWAS.Z, twas_p = TWAS.P, coloc_pp = COLOC.PP4) |>
  arrange(trait, gene_id, phenotype_id)

write_tsv(twas_geuvadis, "data/processed/geuvadis.twas_hits.tsv.gz")

twas_gtex5 <- tibble(tissue = tissues_gtex5) |>
    reframe(
        read_tsv(
            str_glue("data/twas/twas_hits.gtex5-full-{tissue}.tsv"),
            col_types = "cccccccccccccccccccccccc"
        ),
        .by = tissue
    ) |>
    mutate(gene_id = str_split_i(ID, "__", i = 1)) |>
    select(tissue, trait = TRAIT, gene_id, phenotype_id = ID, hsq = HSQ, twas_z = TWAS.Z, twas_p = TWAS.P, coloc_pp = COLOC.PP4) |>
    arrange(tissue, trait, gene_id, phenotype_id)

write_tsv(twas_gtex5, "data/processed/gtex5.twas_hits.tsv.gz")
