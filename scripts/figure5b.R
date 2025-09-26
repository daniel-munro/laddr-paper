## Percentage of xTWAS hits shared with nearby genes

library(tidyverse)
library(GenomicRanges)

has_nearby_hits <- function(chrom, start, end, distance) {
  # Create GRanges object from input vectors
  gr <- GRanges(seqnames = chrom,
                ranges   = IRanges(start, end))
  
  # Find overlaps within the specified distance
  # First, expand each range by the distance on both sides
  gr_expanded <- gr + distance
  
  # Find overlaps between expanded ranges
  overlaps <- findOverlaps(gr_expanded, gr_expanded)
  
  # Remove self-overlaps
  overlaps <- overlaps[queryHits(overlaps) != subjectHits(overlaps)]
  
  has_nearby <- logical(length(gr))
  if (length(overlaps) > 0) {
    has_nearby[queryHits(overlaps)] <- TRUE
  }
  
  return(has_nearby)
}

genes <- read_tsv("data/processed/pcg_and_lncrna.tsv", col_types = "cccciic") |>
  select(gene_id, gene_biotype, chrom, start, end)

twas_dp <- read_tsv("data/processed/geuvadis-full.twas_hits.tsv.gz", col_types = "cccdddd")
twas_kp <- read_tsv("data/processed/geuvadis-pantry.twas_hits.tsv.gz", col_types = "ccccdddd")

twas_pairs <- bind_rows(
  twas_dp |>
    summarise(has_coloc_hit = any(coloc_pp >= 0.8),
              .by = c(trait, gene_id)) |>
    mutate(phenos = "DP", .before = 1),
  twas_kp |>
    summarise(has_coloc_hit = any(coloc_pp >= 0.8),
              .by = c(trait, gene_id)) |>
    mutate(phenos = "KP", .before = 1),
) |>
  left_join(genes, by = "gene_id", relationship = "many-to-one") |>
  mutate(nearby_hit = has_nearby_hits(chrom, start, end, 0),
         .by = c(phenos, trait))

frac_pairs <- twas_pairs |>
  summarise(frac_nearby_hit = mean(nearby_hit),
            ci95_low = frac_nearby_hit - qnorm(0.975) * sd(nearby_hit) / sqrt(n()),
            ci95_hi = frac_nearby_hit + qnorm(0.975) * sd(nearby_hit) / sqrt(n()),
            .by = c(phenos, gene_biotype))

frac_pairs |>
  mutate(gene_biotype = str_replace(gene_biotype, "_", "-") |> fct_rev(),
         phenos = fct_rev(phenos)) |>
  ggplot(aes(x = phenos, y = frac_nearby_hit, ymin = ci95_low, ymax = ci95_hi, fill = phenos)) +
  facet_wrap(~gene_biotype) +
  geom_col(width = 0.5, show.legend = FALSE) +
  geom_errorbar(width = 0.25) +
  expand_limits(y = 1) +
  scale_y_continuous(expand = c(0, 0), labels = scales::label_percent()) +
  scale_fill_manual(values = c("#ad611f", "#13918d")) +
  theme_classic() +
  theme(
    axis.text = element_text(color = "black"),
  ) +
  xlab("xTWAS phenotypes") +
  ylab("Gene-trait pairs w/ nearby same-trait hit")

ggsave("figures/figure3/figure3g.png", width = 2.5, height = 2.9, device = png)

twas_pairs |>
  summarise(fraction = mean(nearby_hit),
            .by = c(phenos, gene_biotype))
