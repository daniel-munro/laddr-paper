## Unique significant TWAS gene-trait pairs for DDP vs KDP phenotypes

library(tidyverse)
library(scales)

summarize_dataset <- function(dat) {
  n_tests <- nrow(dat)
  n_gene_trait_pairs <- dat |>
    distinct(gene_id, trait) |>
    nrow()

  dat <- dat |>
    mutate(
      bh_fdr = p.adjust(TWAS.P, method = "BH"),
      bonferroni_p = p.adjust(TWAS.P, method = "bonferroni")
    )

  thresholds <- tribble(
    ~method, ~threshold_label, ~threshold_value,
    "fixed_p", "p <= 5e-8", 5e-8,
    "bonferroni", "Bonferroni <= 0.05", 0.05,
    "bh_fdr", "BH FDR <= 0.05", 0.05
  )

  pmap_dfr(
    thresholds,
    function(method, threshold_label, threshold_value) {
      if (method == "fixed_p") {
        is_significant <- dat$TWAS.P <= threshold_value
      } else if (method == "bonferroni") {
        is_significant <- dat$bonferroni_p <= threshold_value
      } else {
        is_significant <- dat$bh_fdr <= threshold_value
      }

      hits <- dat |>
        filter(is_significant %in% TRUE)

      tibble(
        threshold_label = threshold_label,
        n_total = n_gene_trait_pairs,
        n_hits = hits |>
          distinct(gene_id, trait) |>
          nrow(),
        proportion = (hits |>
          distinct(gene_id, trait) |>
          nrow()) / n_gene_trait_pairs
      )
    }
  )
}

genes <- read_tsv("data/processed/pcg_and_lncrna.tsv", col_types = "c------") |>
  pull()

df_ddp <- read_tsv("data/twas/twas_pvalues_ddp.tsv.gz", col_types = "ccd") |>
  mutate(gene_id = str_remove(ID, "__.*$")) |>
  filter(gene_id %in% genes)

df_kdp <- read_tsv("data/twas/twas_pvalues_kdp.tsv.gz", col_types = "ccd") |>
  mutate(gene_id = str_remove(ID, "__.*$")) |>
  filter(gene_id %in% genes)

summary_dat <- bind_rows(
  summarize_dataset(df_ddp) |>
    mutate(phenotype_set = "Data-driven phenotypes"),
  summarize_dataset(df_kdp) |>
    mutate(phenotype_set = "Knowledge-driven phenotypes")
) |>
  mutate(
    phenotype_set = factor(phenotype_set),
    threshold_label = factor(
      threshold_label,
      levels = c(
        "p <= 5e-8",
        "Bonferroni <= 0.05",
        "BH FDR <= 0.05"
      )
    )
  )

summary_dat |>
  ggplot(aes(x = threshold_label, y = n_hits / 1000, fill = phenotype_set)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.6) +
  scale_fill_manual(values = c("#13918d", "#ad611f")) +
  scale_y_continuous(
    labels = label_number(big.mark = ","),
    expand = expansion(mult = c(0, 0.04))
  ) +
  theme_bw() +
  theme(
    axis.text = element_text(color = "black"),
    axis.text.x = element_text(angle = 30, hjust = 1),
    legend.background = element_rect(color = "black", linewidth = 0.2),
    legend.position = "inside",
    legend.position.inside = c(0.35, 0.8),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  xlab(NULL) +
  ylab("Unique xTWAS hit gene-trait pairs (×1000)") +
  labs(fill = NULL)

output_dir <- Sys.getenv("FIGURES_DIR", unset = "figures")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

ggsave(
  file.path(output_dir, "figureS3.png"),
  width = 4,
  height = 4,
  device = png
)
