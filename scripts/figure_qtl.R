library(tidyverse)

load_qtls <- function(version) {
    fname <- str_glue("data/residual/{version}.cis_independent_qtl.txt.gz")
    df <- read_tsv(fname, col_types = "c-----c---------dc-i") |>
        rename(gene_id = group_id)
    if (version == "latent-full") {
        mutate(df, modality = "latent_full", .before = phenotype_id)
    } else {
        separate_wider_delim(df, phenotype_id, ":", names = c("modality", "phenotype_id"),
                             too_many = "merge")
    }
}

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

# Use muted version of Pantry colors to deemphasize what is already known
modality_colors <- c(
    Expression = "#bf4042",
    `Isoform ratio` = "#6a90cd",
    `Intron excision` = "#59a257",
    `Alt. TSS` = "#896090",
    `Alt. polyA` = "#d97f26",
    `RNA stability` = "#ddb23c",
    `Latent (full)` = "#13918d",
    `Latent (residual)` = "#1ce6df"
)

versions <- c(
    `all-pantry` = "Explicit",
    `all-residual` = "Explicit + Latent (residual)",
    `latent-full` = "Latent (full)"
)

qtls <- tibble(version = c("all-pantry", "all-residual", "latent-full")) |>
    reframe(load_qtls(version), .by = version) |>
    mutate(modality = factor(modalities[modality], levels = names(modality_colors)),
           version = factor(versions[version], levels = versions))

qtls |>
    count(version, modality) |>
    mutate(modality = fct_rev(modality),
           version = fct_rev(version)) |>
    ggplot(aes(x = n / 1000, y = version, fill = modality)) +
    geom_col(width = 0.8) +
    scale_fill_manual(values = modality_colors) +
    theme_classic() +
    theme(
        legend.key.size = unit(10, "pt"),
    ) +
    xlab("Independent cis-QTLs (×1000)") +
    ylab("RNA phenotypes") +
    labs(fill = "Modality")

ggsave("figures/figure_qtl.png", width = 6, height = 1.8, device = png)

## Number of additional xGenes by adding latent residual to Pantry

qtls |>
    distinct(version, gene_id) |>
    count(version, name = "xGenes")

## Actual overlaps of the xGene sets

qtls |>
    distinct(version, gene_id) |>
    summarise(versions = str_c(version, collapse = " & "),
              .by = gene_id) |>
    count(versions)

qtls |>
    filter(version %in% versions[1:2]) |>
    distinct(version, gene_id) |>
    summarise(versions = str_c(version, collapse = " & "),
              .by = gene_id) |>
    count(versions)
