## TWAS example loci: genotype-stratified RNA-seq coverage

library(tidyverse)
library(patchwork)

examples <- read_tsv(
  "data/analyses/twas_examples/twas_example_variants.tsv",
  col_types = "ccc--"
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
    file.path(example_dir, "raw_coverage.tsv.gz"),
    col_types = cols(
      sample_id = "c",
      individual_id = "c",
      tissue = "c",
      gene_id = "c",
      .default = "d"
    )
  )
  
  bins <- read_tsv(
    file.path(example_dir, "raw_coverage_bins.tsv.gz"),
    col_types = cols(
      raw_bin_id = "c",
      tissue = "c",
      seqname = "c",
      start = "i",
      end = "i",
      gene_id = "c"
    )
  ) |>
    transmute(
      bin_id = raw_bin_id,
      tissue,
      gene_id,
      chrom = seqname,
      chrom_start = start,
      chrom_end = end
    )
  
  zoom_coverage <- read_tsv(
    file.path(example_dir, "zoom_raw_coverage.tsv.gz"),
    col_types = cols(
      sample_id = "c",
      individual_id = "c",
      tissue = "c",
      gene_id = "c",
      .default = "d"
    )
  )
  
  zoom_bins <- read_tsv(
    file.path(example_dir, "zoom_raw_coverage_bins.tsv.gz"),
    col_types = cols(
      raw_bin_id = "c",
      tissue = "c",
      seqname = "c",
      start = "i",
      end = "i",
      gene_id = "c"
    )
  ) |>
    transmute(
      bin_id = raw_bin_id,
      tissue,
      gene_id,
      chrom = seqname,
      chrom_start = start,
      chrom_end = end
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
      starts_with("raw_bin_"),
      names_to = "bin",
      values_to = "coverage"
    ) |>
    left_join(
      bins |>
        mutate(
          bin = bin_id,
          x = (chrom_start + chrom_end) / 2e6
        ) |>
        select(bin, x),
      by = "bin"
    ) |>
    mutate(
      coverage = log10(coverage + 1),
      genotype_label = factor(genotype_label, levels = genotype_levels),
      genotype_group = factor(genotype_label, levels = genotype_group_levels)
    ) |>
    summarise(
      covg_25th = quantile(coverage, 0.25),
      covg_50th = median(coverage),
      covg_75th = quantile(coverage, 0.75),
      covg_mean = mean(coverage),
      .by = c(bin, x, genotype_label, genotype_group)
    )

  gene_name <- gene_names[[gene_id]]

  genotype_colors <- set_names(
    c("#2c7fb8", "#6d6d6d", "#d95f0e")[seq_along(genotype_levels)],
    genotype_levels
  )
  
  genotype_labels <- genotype_counts |>
    mutate(label = str_glue("{genotype_label} (n={n})")) |>
    with(set_names(label, genotype_label))

  plot_coverage <- function(coverage_summary, x, xlab, title, show_x_axis = FALSE,
                            show_legend = TRUE, x_limits = NULL) {
    coverage_summary |>
    ggplot(aes(
      x = {{ x }},
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
      axis.text.x = if (show_x_axis) element_text(color = "black") else element_blank(),
      axis.ticks.x = if (show_x_axis) element_line() else element_blank(),
      legend.key.height = unit(10, "pt"),
      legend.position = if (show_legend) "inside" else "none",
      legend.position.inside = c(0.78, 0.7),
      legend.title = element_blank(),
      panel.grid = element_blank(),
    ) +
    coord_cartesian(xlim = x_limits, expand = FALSE) +
    xlab(xlab) +
    ylab("log10(coverage+1)") +
    labs(title = title)
  }
  
  coverage_plot <- plot_coverage(
    coverage_summary,
    x,
    str_glue("{bins$chrom[[1]]} position (Mb)"),
    str_glue("{gene_name} in {tissue} grouped by {variant_label(variant_id)}"),
    show_x_axis = TRUE,
    x_limits = range(c(bins$chrom_start, bins$chrom_end)) / 1e6
  ) +
    scale_x_continuous(labels = scales::label_number()) +
    ylab("log10(raw coverage+1)")
  
  zoom_summary <- zoom_coverage |>
    left_join(genotypes, by = "individual_id") |>
    filter(!is.na(genotype_label)) |>
    pivot_longer(
      starts_with("raw_bin_"),
      names_to = "bin",
      values_to = "coverage"
    ) |>
    left_join(
      zoom_bins |>
        mutate(
          bin = bin_id,
          x = (chrom_start + chrom_end) / 2e6
        ) |>
        select(bin, x),
      by = "bin"
    ) |>
    mutate(
      coverage = log10(coverage + 1),
      genotype_label = factor(genotype_label, levels = genotype_levels),
      genotype_group = factor(genotype_label, levels = genotype_group_levels)
    ) |>
    summarise(
      covg_mean = mean(coverage),
      .by = c(bin, x, genotype_label, genotype_group)
    )
  
  zoom_start <- min(zoom_bins$chrom_start)
  zoom_end <- max(zoom_bins$chrom_end)
  zoom_plot <- plot_coverage(
    zoom_summary,
    x,
    str_glue("{zoom_bins$chrom[[1]]} position (Mb)"),
    str_glue("Zoomed to {zoom_bins$chrom[[1]]}:{zoom_start}-{zoom_end}"),
    show_x_axis = TRUE,
    show_legend = FALSE,
    x_limits = c(zoom_start, zoom_end) / 1e6
  ) +
    scale_x_continuous(labels = scales::label_number()) +
    ylab("log10(raw coverage+1)")
  
  coverage_plot + zoom_plot +
    plot_layout(widths = c(2, 1))
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
  file.path(output_dir, "figureS5.png"),
  width = 11,
  height = 8,
  device = png
)
