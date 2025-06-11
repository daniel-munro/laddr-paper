## Comparison of binning and model fitting methods

library(tidyverse)
library(patchwork)

n_genes <- 932

stats <- read_tsv("data/compare/bin_pheno_qtl_counts.tsv", col_types = "ccciiii") |>
  mutate(
    bins_per_gene = n_bins / n_genes,
    phenos_per_gene = n_phenos / n_genes,
    phenos_tested_per_gene = n_phenos_tested / n_genes,
    qtls_per_gene = n_qtls / n_genes
  )

stats_avg <- stats |>
  summarise(
    bins_per_gene_mean = mean(bins_per_gene),
    bins_per_gene_sd = sd(bins_per_gene),
    phenos_per_gene_mean = mean(phenos_per_gene),
    phenos_per_gene_sd = sd(phenos_per_gene),
    qtls_per_gene_mean = mean(qtls_per_gene),
    qtls_per_gene_sd = sd(qtls_per_gene),
    .by = c(type, param)
  )

## Panel a: Binning

bin_types <- c(
  n = "N per coding/non-\ncoding region",
  w = "Fixed coding/non-\ncoding bin widths",
  a1 = "Adaptive: cumulative\ncovg. / max. corr.",
  a2 = "Adaptive: cumulative\nvariance of covg.",
  a3 = "Adaptive: cumulative\nvariance of covg. diff."
)

bin_colors <- c(
  `N per coding/non-\ncoding region` = "#169740",
  `Fixed coding/non-\ncoding bin widths` = "#3974d5",
  `Adaptive: cumulative\ncovg. / max. corr.` = "#fa6410",
  `Adaptive: cumulative\nvariance of covg.` = "#eb3499",
  `Adaptive: cumulative\nvariance of covg. diff.` = "#b264ed"
)

stats_avg_bin <- stats_avg |>
  filter(type == "binning") |>
  select(-type) |>
  separate_wider_delim(param, "_", names = c("bin_type", "param_number"),
                       too_many = "merge") |>
  mutate(bin_type = factor(bin_types[bin_type], levels = bin_types),
         param_number = str_replace(param_number, "_", ","))

p1 <- stats_avg_bin |>
  ggplot(aes(x = bins_per_gene_mean,
             xmin = bins_per_gene_mean - bins_per_gene_sd,
             xmax = bins_per_gene_mean + bins_per_gene_sd,
             y = qtls_per_gene_mean,
             ymin = qtls_per_gene_mean - qtls_per_gene_sd,
             ymax = qtls_per_gene_mean + qtls_per_gene_sd,
             color = bin_type,
             label = param_number)) +
  geom_linerange(orientation = "x") +
  geom_linerange(orientation = "y") +
  geom_point() +
  ggrepel::geom_text_repel(size = 3, show.legend = FALSE) +
  scale_color_manual(values = bin_colors) +
  expand_limits(x = 0) +
  theme_bw() +
  theme(
    legend.key.spacing.y = unit(8, "pt"),
    panel.grid = element_blank(),
  ) +
  xlab("Bins per gene (mean)") +
  ylab("xQTLs per gene (mean)") +
  labs(color = "Binning method") +
  ggtitle("Binning variations")
p1

## Panel b: Model fitting

model_types = c(
  pca = "PCA (max PCs)",
  fpca_xpos_disc = "FPCA discrete x=pos",
  fpca_xpos_spline = "FPCA spline x=pos",
  fpca_xbin_spline = "FPCA spline x=bin"
)

model_colors <- c(
  "PCA (max PCs)" = "#ff5b4f",
  "FPCA discrete x=pos" = "#b3e36b",
  "FPCA spline x=pos" = "#9ea1ff",
  "FPCA spline x=bin" = "#f2d968"
)

model_shapes <- c(
  "PCA (max PCs)" = 19,
  "FPCA discrete x=pos" = 3,
  "FPCA spline x=pos" = 0,
  "FPCA spline x=bin" = 5
)

stats_avg_model <- stats_avg |>
  filter(type == "model") |>
  select(-type, -bins_per_gene_mean, -bins_per_gene_sd) |>
  separate_wider_delim(param, "_n", names = c("model_type", "param_number"), too_few = "align_start") |>
  mutate(model_type = factor(model_types[model_type], levels = model_types)) |>
  replace_na(list(param_number = ""))

p2 <- stats_avg_model |>
  ggplot(aes(x = phenos_per_gene_mean,
             xmin = phenos_per_gene_mean - phenos_per_gene_sd,
             xmax = phenos_per_gene_mean + phenos_per_gene_sd,
             y = qtls_per_gene_mean,
             ymin = qtls_per_gene_mean - qtls_per_gene_sd,
             ymax = qtls_per_gene_mean + qtls_per_gene_sd,
             color = model_type,
             shape = model_type,
             label = param_number)) +
  geom_linerange(orientation = "x") +
  geom_point(stroke = 1) +
  ggrepel::geom_text_repel(size = 3, show.legend = FALSE) +
  scale_color_manual(values = model_colors) +
  scale_shape_manual(values = model_shapes) +
  expand_limits(x = 0) +
  theme_bw() +
  theme(
    axis.text = element_text(color = "black"),
    panel.grid = element_blank(),
  ) +
  xlab("Latent phenotypes per gene (mean)") +
  ylab("xQTLs per gene (mean)") +
  labs(color = "Model type", shape = "Model type") +
  ggtitle("Model variations")
p2

p1 / p2 +
  plot_layout(heights = c(3, 2)) +
  plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(face = "bold"))

ggsave("figures/figureS1.png", width = 8, height = 8, device = png)
