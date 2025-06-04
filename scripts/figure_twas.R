library(tidyverse)

modalities <- c(
    expression = "Expression",
    isoforms = "Isoform ratio",
    splicing = "Intron excision",
    alt_TSS = "Alt. TSS",
    alt_polyA = "Alt. polyA",
    stability = "RNA stability",
    latent = "Latent"
)

# Use muted version of Pantry colors to deemphasize what is already known
modality_colors <- c(
    Expression = "#bf4042",
    `Isoform ratio` = "#6a90cd",
    `Intron excision` = "#59a257",
    `Alt. TSS` = "#896090",
    `Alt. polyA` = "#d97f26",
    `RNA stability` = "#ddb23c",
    `Latent` = "#13918d"
)

categories <- read_tsv("data/pantry/geuvadis/twas/gwas_metadata.txt",
                   col_types = cols(Tag = "c", Category = "c", .default = "-")) |>
    mutate(Category = fct_lump_min(Category, 10)) |>
    deframe()

twas <- read_tsv("data/twas/twas_hits.geuvadis-full-Geuvadis.tsv",
                 col_types = "cc---d-----------ddddddd") |>
    separate_wider_delim(ID, "__", names = c("gene_id", "PC"), cols_remove = FALSE) |>
    rename(trait = TRAIT)

twas_pantry <- read_tsv("data/pantry/processed/geuvadis.twas.tsv.gz", col_types = "ccccdcdd") |>
    filter(TWAS.P < 5e-8)

df <- bind_rows(
    twas |>
        select(trait, gene_id) |>
        mutate(trait = as.character(trait)) |>
        mutate(phenos = "Latent", modality = "latent", .before = 1),
    twas_pantry |>
        select(trait, gene_id, modality) |>
        mutate(phenos = "Explicit", .before = 1),
) |>
    mutate(
        category = categories[trait] |>
            fct_infreq(),
        modality = factor(modalities[modality], levels = modalities),
    )

df |>
    ggplot(aes(x = category, fill = modality)) +
    facet_wrap(~ phenos) +
    geom_bar() +
    scale_fill_manual(values = modality_colors) +
    theme_bw() +
    theme(
        axis.text.x = element_text(hjust = 1, vjust = 1, angle = 45),
        legend.position = "inside",
        legend.position.inside = c(0.32, 0.75),
        legend.key.size = unit(10, "pt"),
        panel.grid = element_blank(),
    ) +
    labs(fill = "Modality") +
    xlab("Trait category") +
    ylab("TWAS hits")

ggsave("figures/figure_twas.png", width = 4.5, height = 4.5, device = png)
