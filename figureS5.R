## TWAS example loci: genotype-stratified RNA-seq coverage

library(tidyverse)
library(patchwork)

examples <- read_tsv(
  "data/analyses/twas_examples/twas_example_variants.tsv",
  col_types = "cccc--"
)

gene_names <- read_tsv("data/processed/pcg_and_lncrna.tsv", col_types = "cc-----") |>
  tibble::deframe()

trait_names <- read_tsv(
  "data/twas/gwas_metadata.txt",
  col_types = cols(Tag = "c", Phenotype = "c", .default = "-")
) |>
  select(Tag, Phenotype) |>
  mutate(Phenotype = str_replace(Phenotype, "Triglycerids", "Triglycerides")) |>
  tibble::deframe()

variant_label <- function(variant_id) {
  parts <- str_split_fixed(variant_id, "_", 5)
  str_glue("{parts[, 1]}:{parts[, 2]}")
}

plot_example <- function(tissue, gene_id, variant_id, trait_id) {
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
      covg_mean = mean(coverage),
      .by = c(bin, x, genotype_label, genotype_group)
    )
  
  gene_name <- gene_names[[gene_id]]
  trait_name <- trait_names[[trait_id]]
  
  genotype_colors <- set_names(
    c("#2c7fb8", "#6d6d6d", "#d95f0e")[seq_along(genotype_levels)],
    genotype_levels
  )
  
  genotype_labels <- genotype_counts |>
    mutate(label = str_glue("{genotype_label} (n={n})")) |>
    with(set_names(label, genotype_label))
  
  zoom_start <- min(zoom_bins$chrom_start)
  zoom_end <- max(zoom_bins$chrom_end)
  zoom_x_limits <- c(zoom_start, zoom_end) / 1e6

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
  
  y_limits <- c(0, max(coverage_summary$covg_mean, zoom_summary$covg_mean))
  
  coverage_plot <- coverage_summary |>
    ggplot(aes(
      x = x,
      y = covg_mean,
      color = genotype_label,
      fill = genotype_label,
      group = genotype_group
    )) +
    annotate(
      "rect",
      xmin = zoom_x_limits[[1]],
      xmax = zoom_x_limits[[2]],
      ymin = -Inf,
      ymax = Inf,
      color = "black",
      fill = "#eeeeee",
      linewidth = 0.1
    ) +
    geom_line(linewidth = 0.5, show.legend = FALSE) +
    scale_x_continuous(expand = FALSE) +
    scale_y_continuous(limits = y_limits, expand = expansion(mult = c(0, 0.04))) +
    scale_color_manual(values = genotype_colors, labels = genotype_labels) +
    scale_fill_manual(values = genotype_colors) +
    theme_bw() +
    theme(
      axis.text = element_text(color = "black"),
      panel.grid = element_blank(),
      plot.margin = margin(0, 8, 0, 0, unit = "pt"),
      plot.title = element_text(size = 12),
    ) +
    labs(
      x = str_glue("{bins$chrom[[1]]} position (Mb)"),
      y ="log10(coverage+1)",
      title = str_glue("{gene_name} in {tissue} grouped by {variant_label(variant_id)}"),
      subtitle = str_glue("Associated with {trait_name}"),
    )
  
  zoom_plot <- zoom_summary |>
    ggplot(aes(
      x = x,
      y = covg_mean,
      color = genotype_label,
      fill = genotype_label,
      group = genotype_group
    )) +
    geom_line(linewidth = 0.5) +
    scale_x_continuous(expand = FALSE) +
    scale_y_continuous(limits = y_limits, expand = expansion(mult = c(0, 0.04))) +
    scale_color_manual(values = genotype_colors, labels = genotype_labels) +
    scale_fill_manual(values = genotype_colors) +
    theme_bw() +
    theme(
      axis.text = element_text(color = "black"),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      legend.key.height = unit(10, "pt"),
      legend.box.spacing = unit(0, "pt"),
      legend.position = "right",
      legend.title = element_blank(),
      panel.grid = element_blank(),
      plot.margin = margin_auto(0, unit = "pt"),
    ) +
    labs(
      subtitle = "Zoomed to highest-difference region",
      x = str_glue("{zoom_bins$chrom[[1]]} position (Mb)"),
      y = NULL,
    )
  
  wrap_elements(
    coverage_plot + zoom_plot +
      plot_layout(widths = c(1, 1))
  )
}

plots <- pmap(examples, plot_example)

combined_plot <- wrap_plots(plots, ncol = 1) +
  plot_annotation(tag_levels = "A", theme = theme(plot.margin = margin(0, 0, 0, 0))) &
  theme(
    plot.margin = margin_auto(5, unit = "pt"),
    plot.tag = element_text(face = "bold", margin = margin(0, 1, 0, 0, unit = "pt"))
  )

if (interactive()) {
  print(combined_plot)
}

output_dir <- Sys.getenv("FIGURES_DIR", unset = "figures")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

ggsave(
  file.path(output_dir, "figureS5.png"),
  width = 9.5,
  height = 8,
  device = png
)
