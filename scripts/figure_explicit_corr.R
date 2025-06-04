## For each latent phenotype, full and residual, maximum Pearson r2 between each
## latent phenotype and an explicit phenotype of the same gene.

library(tidyverse)

modalities <- c(
    expression = "Expression",
    isoforms = "Isoform ratio",
    splicing = "Intron excision",
    alt_TSS = "Alt. TSS",
    alt_polyA = "Alt. polyA",
    stability = "RNA stability",
    latent_full = "Latent (full)",
    latent_residual = "Latent (residual)"
)

modality_colors <- c(
    Expression = "#e41a1c",
    `Isoform ratio` = "#377eb8",
    `Intron excision` = "#4daf4a",
    `Alt. TSS` = "#984ea3",
    `Alt. polyA` = "#ff7f00",
    `RNA stability` = "#fdc11c",
    `Latent (full)` = "#1bf2eb",
    `Latent (residual)` = "#a4dedc"
)

latent_types = c(
    full = "Full",
    residual = "Residual",
    null = "Null"
)

latent_colors <- c(
    Full = "#1bf2eb",
    Residual = "#a4dedc",
    Null = "white"
)

corrs_max <- read_tsv("data/processed/latent_explicit_corrs.tsv.gz", col_types = "ccccd") |>
    separate(pheno, c("modality", "pheno"), sep = ":", extra = "merge") |>
    mutate(modality = factor(modalities[modality], levels = modalities),
           PC = fct_inorder(PC),
           r2_max = r^2)

corrs_max |>
    mutate(latent = factor(latent_types[latent], levels = latent_types)) |>
    ggplot(aes(x = PC, y = r2_max, fill = latent)) +
    geom_boxplot(outlier.size = 0.25) +
    scale_fill_manual(values = latent_colors) +
    theme_classic() +
    theme(
        legend.position = "inside",
        legend.position.inside = c(0.8, 0.8),
    ) +
    xlab("Latent RNA phenotype rank per gene") +
    ylab(expression("Maximum "*r^2*" to explicit phenotype")) +
    labs(fill = "Latent type")

ggsave("figures/figure_explicit_corr.png", width = 7, height = 5, device = png)

## Split by modality of max correlation

corrs_max |>
    filter(latent == "full",
           PC %in% str_c("PC", 1:8)) |>
    ggplot(aes(x = PC, y = r2_max, fill = modality)) +
    geom_boxplot(outlier.size = 0.25) +
    scale_fill_manual(values = modality_colors) +
    theme_classic() +
    theme(
        legend.position = "inside",
        legend.position.inside = c(0.8, 0.8),
    ) +
    xlab("Latent RNA phenotype rank per gene") +
    ylab(expression("Maximum "*r^2*" to explicit phenotype")) +
    labs(fill = "Modality")

ggsave("figures/figure_explicit_corr_mod.png", width = 5, height = 3, device = png)

## Modality median lines overlaid on latent types

corrs_max_mod <- corrs_max |>
    filter(latent == "full") |>
    summarise(med_r2_max = median(r2_max),
              .by = c(PC, modality))

corrs_max |>
    mutate(latent = factor(latent_types[latent], levels = latent_types)) |>
    ggplot(aes(x = PC, y = r2_max, fill = latent)) +
    geom_boxplot(outlier.size = 0.25) +
    geom_line(aes(x = as.integer(PC) - 0.25, y = med_r2_max, color = modality, group = modality, fill = NULL), data = corrs_max_mod) +
    geom_point(aes(x = as.integer(PC) - 0.25, y = med_r2_max, color = modality, fill = NULL), data = corrs_max_mod) +
    scale_fill_manual(values = latent_colors) +
    scale_color_manual(values = modality_colors) +
    theme_classic() +
    theme(
        legend.position = "inside",
        legend.position.inside = c(0.8, 0.8),
    ) +
    xlab("Latent RNA phenotype rank per gene") +
    ylab(expression("Maximum "*r^2*" to explicit phenotype")) +
    labs(fill = "Modality")

ggsave("figures/figure_explicit_corr_mod2.png", width = 7, height = 5, device = png)

## As heatmap

corrs_max_mod |>
    filter(PC %in% str_c("PC", 1:8)) |>
    ggplot(aes(x = PC, y = modality, fill = med_r2_max)) +
    geom_tile() +
    scale_fill_gradient(low = "white", high = "#555555") +
    expand_limits(fill = 0) +
    theme_classic() +
    theme(
        legend.position = "top",
    ) +
    xlab("Latent RNA phenotype rank per gene") +
    ylab("Modality") +
    labs(fill = expression("Median max "*r^2*"\nto explicit phenotype"))

ggsave("figures/figure_explicit_corr_mod3.png", width = 5, height = 3, device = png)

###################################
## Calculate max r2 per modality ##
###################################

corrs_max_mod2 <- read_tsv("data/processed/latent_explicit_corrs_mod.tsv.gz", col_types = "cccd") |>
    separate(pheno, c("modality", "pheno"), sep = ":", extra = "merge") |>
    mutate(modality = factor(modalities[modality], levels = modalities),
           PC = fct_inorder(PC),
           r2_max = r^2)

corrs_max_mod2 |>
    filter(PC %in% str_c("PC", 1:8)) |>
    ggplot(aes(x = PC, y = r2_max, fill = modality)) +
    facet_wrap(~ modality, ncol = 3) +
    geom_boxplot(outlier.size = 0.25, show.legend = FALSE) +
    scale_fill_manual(values = modality_colors) +
    theme_classic() +
    theme(
        legend.position = "inside",
        legend.position.inside = c(0.8, 0.8),
    ) +
    xlab("Latent RNA phenotype rank") +
    ylab(expression("Maximum "*r^2*" to explicit phenotype")) +
    ggtitle("Max corr to each modality (calculated separately)")

ggsave("figures/figure_explicit_corr_mod4.png", width = 7, height = 5, device = png)
