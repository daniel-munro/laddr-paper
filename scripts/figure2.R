## Latent RNA phenotype characteristics

library(tidyverse)
library(patchwork)

## Panel a: Correlation heatmap example

## Panel b: Latent phenotype cis-heritability

modalities <- c(
  expression = "Expression",
  isoforms = "Isoform ratio",
  splicing = "Intron excision",
  alt_TSS = "Alt. TSS",
  alt_polyA = "Alt. polyA",
  stability = "RNA stability",
  latent = "Latent"
)

hsq_latent <- read_tsv("data/twas/Geuvadis-latent.profile", col_types = "cidddddd-dddd-d") |>
  separate_wider_delim(id, "__", names = c("gene_id", "PC"), cols_remove = FALSE) |>
  mutate(PC = str_replace(PC, "PC", "") |> as.integer())

hsq_pantry <- read_tsv("data/pantry/processed/geuvadis.hsq.tsv.gz", col_types = "ccciddddddddddd") |>
  mutate(modality = factor(modalities[modality], levels = modalities) |>
           fct_reorder(hsq, .desc = TRUE))

p1 <- hsq_latent |>
  mutate(PC = as.character(PC) |>
           fct_reorder(PC) |>
           fct_lump_n(n = 10, other_level = "11+")) |>
  ggplot(aes(x = PC, y = hsq)) +
  geom_boxplot() +
  scale_y_log10(minor_breaks = c(0.01, 0.02, 0.03)) +
  expand_limits(y = 0.01) +
  theme_bw() +
  theme(
    axis.text = element_text(color = "black"),
    panel.grid = element_blank(),
    plot.margin = margin(5.5, 0, 5.5, 5.5),
  ) +
  ylab(expression(h^2)) +
  ggtitle("Latent phenotypes")

p2 <- hsq_pantry |>
  ggplot(aes(x = modality, y = hsq)) +
  geom_boxplot() +
  scale_y_log10() +
  expand_limits(y = 0.01) +
  theme_bw() +
  theme(
    axis.text.x = element_text(hjust = 1, vjust = 1, angle = 45, color = "black"),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    panel.grid = element_blank(),
    plot.margin = margin(5.5, 5.5, 5.5, 0),
  ) +
  xlab("Modality") +
  ylab(NULL) +
  ggtitle("Explicit phenotypes")

p1 + p2 + plot_layout(widths = c(11, 6))

ggsave("figures/figure2/figure2b.png", width = 5, height = 4, device = png)
