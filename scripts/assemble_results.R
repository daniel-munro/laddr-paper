library(tidyverse)

tissues_gtex5 <- read_lines("data/info/tissues.gtex5.txt")

##########
## QTLs ##
##########

qtl_assoc <- tibble(tissue = tissues_gtex5) |>
    reframe(
        read_tsv(
            str_glue("data/qtl/gtex5-full-{tissue}.cis_qtl.txt.gz"),
            col_types = "c-----c---------cc-c-"
        ),
        .by = tissue
    ) |>
    rename(gene_id = group_id) |>
    relocate(gene_id, .before = 2)

write_tsv(qtl_assoc, "data/processed/gtex5.qtl_assoc.tsv.gz")

qtls <- tibble(tissue = tissues_gtex5) |>
    reframe(
        read_tsv(
            str_glue("data/qtl/gtex5-full-{tissue}.cis_independent_qtl.txt.gz"),
            col_types = "c-----c---------cc-i"
        ),
        .by = tissue
    ) |>
    select(tissue, gene_id = group_id, rank, phenotype_id, variant_id, pval_beta)

write_tsv(qtls, "data/processed/gtex5.qtls.tsv.gz")

##########
## TWAS ##
##########

twas <- tibble(tissue = tissues_gtex5) |>
    reframe(
        read_tsv(
            str_glue("data/twas/twas_hits.gtex5-full-{tissue}.tsv"),
            col_types = "cccccccccccccccccccccccc"
        ),
        .by = tissue
    ) |>
    mutate(gene_id = str_replace(ID, ":.*$", "")) |>
    select(tissue, trait = TRAIT, gene_id, phenotype_id = ID, HSQ, MODEL, TWAS.Z, TWAS.P, COLOC.PP0, COLOC.PP1, COLOC.PP2, COLOC.PP3, COLOC.PP4) |>
    arrange(tissue, trait, gene_id, phenotype_id)

write_tsv(twas, "data/processed/gtex5.twas_hits.tsv.gz")
