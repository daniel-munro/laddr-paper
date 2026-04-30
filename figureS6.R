## TWAS example loci: genotype-stratified RNA-seq coverage

library(tidyverse)
library(patchwork)

examples <- read_tsv(
  "data/analyses/twas_examples/twas_example_variants.tsv",
  col_types = cols(
    tissue = "c",
    gene_id = "c",
    variant_id = "c"
  )
)

gene_names <- read_tsv("data/processed/pcg_and_lncrna.tsv", col_types = "cc-----") |>
  tibble::deframe()

variant_label <- function(variant_id) {
  parts <- str_split_fixed(variant_id, "_", 5)
  str_glue("{parts[, 1]}:{parts[, 2]}")
}

plot_example <- function(tissue, gene_id, variant_id) {
  example_dir <- file.path(
    "data/analyses/twas_examples",
    str_glue("{tissue}.{gene_id}.{variant_id}")
  )
  
  coverage <- read_tsv(
    file.path(example_dir, "coverage.tsv.gz"),
    col_types = cols(
      sample_id = "c",
      individual_id = "c",
      tissue = "c",
      gene_id = "c",
      .default = "d"
    )
  )
  
  bins <- read_tsv(
    file.path(example_dir, "bins.tsv.gz"),
    col_types = cols(
      bin_id = "c",
      tissue = "c",
      gene_id = "c",
      pos = "i",
      chrom = "c",
      chrom_start = "i",
      chrom_end = "i"
    )
  )
  
  genotypes <- read_tsv(
    file.path(example_dir, "individual_genotypes.tsv.gz"),
    col_types = cols(
      individual_id = "c",
      variant_id = "c",
      genotype = "c",
      genotype_label = "c",
      .default = "-"
    )
  ) |>
    select(individual_id, genotype, genotype_label)
  
  genotype_levels <- genotypes |>
    distinct(genotype, genotype_label) |>
    mutate(genotype = factor(genotype, levels = c("0/0", "0/1", "1/1"))) |>
    arrange(genotype) |>
    pull(genotype_label)
  
  genotype_group_levels <- c(genotype_levels[1], genotype_levels[3], genotype_levels[2]) |>
    discard(is.na)
  
  genotype_counts <- coverage |>
    left_join(genotypes, by = "individual_id") |>
    filter(!is.na(genotype_label)) |>
    distinct(individual_id, genotype_label) |>
    count(genotype_label) |>
    mutate(genotype_label = factor(genotype_label, levels = genotype_levels)) |>
    arrange(genotype_label)

  coverage_summary <- coverage |>
    left_join(genotypes, by = "individual_id") |>
    filter(!is.na(genotype_label)) |>
    pivot_longer(
      starts_with("bin_"),
      names_to = "bin",
      values_to = "coverage"
    ) |>
    mutate(
      coverage = log10(coverage + 1),
      bin = parse_number(bin) + 1,
      genotype_label = factor(genotype_label, levels = genotype_levels),
      genotype_group = factor(genotype_label, levels = genotype_group_levels)
    ) |>
    summarise(
      covg_25th = quantile(coverage, 0.25),
      covg_50th = median(coverage),
      covg_75th = quantile(coverage, 0.75),
      covg_mean = mean(coverage),
      .by = c(bin, genotype_label, genotype_group)
    )

  gene_name <- gene_names[[gene_id]]

  genotype_colors <- set_names(
    c("#2c7fb8", "#6d6d6d", "#d95f0e")[seq_along(genotype_levels)],
    genotype_levels
  )
  
  genotype_labels <- genotype_counts |>
    mutate(label = str_glue("{genotype_label} (n={n})")) |>
    with(set_names(label, genotype_label))

  coverage_plot <- coverage_summary |>
    ggplot(aes(
      x = bin,
      y = covg_mean,
      color = genotype_label,
      fill = genotype_label,
      group = genotype_group
    )) +
    geom_line(linewidth = 0.5) +
    scale_color_manual(values = genotype_colors, labels = genotype_labels) +
    scale_fill_manual(values = genotype_colors) +
    theme_bw() +
    theme(
      axis.text = element_text(color = "black"),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      legend.key.height = unit(10, "pt"),
      legend.position = "inside",
      legend.position.inside = c(0.8, 0.7),
      legend.title = element_blank(),
      panel.grid = element_blank(),
    ) +
    xlab("Bins along gene start to end") +
    ylab("log10(coverage+1)") +
    labs(
      title = str_glue("{gene_name} in {tissue} grouped by {variant_label(variant_id)}")
    )
  
  coverage_plot
}

plots <- pmap(examples, plot_example)

combined_plot <- wrap_plots(plots, ncol = 1) +
  plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(face = "bold"))

if (interactive()) {
  print(combined_plot)
}

output_dir <- Sys.getenv("FIGURES_DIR", unset = "figures")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

ggsave(
  file.path(output_dir, "figureS6.png"),
  width = 8,
  height = 8,
  device = png
)
