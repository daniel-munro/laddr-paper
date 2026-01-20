## For each latent phenotype, full and residual, look at the maximum Pearson r2
## between each latent phenotype and an explicit phenotype of the same gene. For
## comparison, also shuffle the values for each latent-full phenotype to create
## a latent-null version.

library(tidyverse)

regress_out_pcs <- function(phenos, n_pcs) {
    phenos <- t(as.matrix(phenos))
    pca <- prcomp(phenos, center = FALSE, scale. = FALSE)
    pcs <- pca$x[, seq_len(n_pcs)]
    
    ## Project the data onto the PC sub-space and subtract to get residuals
    ## (projection matrix H = PCs %*% (t(PCs) PCs)^-1 %*% t(PCs))
    XtX_inv <- solve(crossprod(pcs))          # (t(PCs) %*% PCs)^-1
    hat_mat <- pcs %*% XtX_inv %*% t(pcs)     # projection onto PC space
    residual <- phenos - hat_mat %*% phenos
    t(residual)
}

all_corrs <- function(df, gene_id) {
    df1 <- df |>
        filter(modality == "latent")
    x1 <- as.data.frame(select(df1, -phenotype_id, -modality))
    rownames(x1) <- df1$phenotype_id
    df2 <- df |>
        filter(modality != "latent")
    x2 <- as.data.frame(select(df2, -phenotype_id, -modality))
    rownames(x2) <- df2$phenotype_id
    if (nrow(df1) * nrow(df2) == 0) { return(tibble()) }
    cor(t(x1), t(x2)) |>
        as_tibble(rownames = "PC") |>
        pivot_longer(-PC, names_to = "pheno", values_to = "r")
}

pantry <- tibble(modality = c("expression", "isoforms", "splicing", "alt_TSS", "alt_polyA", "stability")) |>
    reframe(
        read_tsv(
            str_glue("data/pantry_phenos/geuvadis/{modality}.bed.gz"),
            col_types = cols(`#chr` = "-", start = "-", end = "-", phenotype_id = "c", .default = "d")
        ),
        .by = modality
    ) |>
    mutate(gene_id = str_split_i(phenotype_id, "__", 1), .before = phenotype_id) |>
    mutate(phenotype_id = str_c(modality, ":", phenotype_id)) |>
    slice_sample(prop = 1, replace = FALSE) |> # In case of ties for max corr later, choose randomly
    arrange(gene_id)

latent_full <- read_tsv(
    "data/phenos/geuvadis-full/Geuvadis-latent.bed.gz",
    col_types = cols(`#chr` = "-", start = "-", end = "-", phenotype_id = "c", .default = "d")
) |>
    separate_wider_delim(phenotype_id, "__", names = c("gene_id", "phenotype_id")) |>
    # mutate(PC = str_replace(phenotype_id, "PC", "") |> as.integer(), .after = phenotype_id) |>
    # filter(PC <= 8) |>
    # select(-PC) |>
    mutate(modality = "latent", .before = 1)

latent_residual <- read_tsv(
    "data/phenos/geuvadis-residual/Geuvadis-latent.bed.gz",
    col_types = cols(`#chr` = "-", start = "-", end = "-", phenotype_id = "c", .default = "d")
) |>
    separate_wider_delim(phenotype_id, "__", names = c("gene_id", "phenotype_id")) |>
    mutate(modality = "latent", .before = 1)

## Regress out covariates before comparing
pantry[, 4:ncol(pantry)] <- regress_out_pcs(pantry[, 4:ncol(pantry)], n_pcs = 4)
latent_full[, 4:ncol(latent_full)] <- regress_out_pcs(latent_full[, 4:ncol(latent_full)], n_pcs = 4)
latent_residual[, 4:ncol(latent_residual)] <- regress_out_pcs(latent_residual[, 4:ncol(latent_residual)], n_pcs = 4)

# ## To subset to a number of random genes:
# genes <- bind_rows(
#     pantry |>
#         distinct(gene_id) |>
#         mutate(group = "pantry"),
#     latent_full |>
#         distinct(gene_id) |>
#         mutate(group = "full"),
#     latent_residual |>
#         distinct(gene_id) |>
#         mutate(group = "residual")
# ) |>
#     filter(n() == 3, .by = gene_id) |>
#     distinct(gene_id) |>
#     sample_n(1000, replace = FALSE) |>
#     pull() |>
#     sort()
# pantry <- pantry |> filter(gene_id %in% genes)
# latent_full <- latent_full |> filter(gene_id %in% genes)
# latent_residual <- latent_residual |> filter(gene_id %in% genes)

latent_null <- latent_full
latent_null[, 4:ncol(latent_null)] <- t(apply(latent_null[, 4:ncol(latent_null)], 1, sample))

corrs_full <- bind_rows(pantry, latent_full) |>
    group_by(gene_id) |>
    group_modify(all_corrs) |>
    ungroup()

corrs_residual <- bind_rows(pantry, latent_residual) |>
    group_by(gene_id) |>
    group_modify(all_corrs) |>
    ungroup()

corrs_null <- bind_rows(pantry, latent_null) |>
    group_by(gene_id) |>
    group_modify(all_corrs) |>
    ungroup()

corrs_max <- bind_rows(
    corrs_full |> mutate(latent = "full", .before = 1),
    corrs_residual |> mutate(latent = "residual", .before = 1),
    corrs_null |> mutate(latent = "null", .before = 1),
) |>
    mutate(r2 = r^2) |>
    slice_max(r2, n = 1, by = c("latent", "gene_id", "PC"), with_ties = FALSE) |>
    select(-r2)

write_tsv(corrs_max, "data/processed/latent_explicit_corrs.tsv.gz")

## For full latent, also save max per modality

corrs_max_mod <- corrs_full |>
    mutate(modality = str_replace(pheno, ":.*$", ""),
           r2 = r^2) |>
    slice_max(r2, n = 1, by = c("gene_id", "PC", "modality"), with_ties = FALSE) |>
    select(-r2, -modality)

write_tsv(corrs_max_mod, "data/processed/latent_explicit_corrs_mod.tsv.gz")
