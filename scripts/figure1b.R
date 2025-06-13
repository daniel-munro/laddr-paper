## Example latent RNA phenotype and TWAS association

library(tidyverse)
library(patchwork)

bin_in_exon <- function(gene, rel_posns, anno) {
  gene_anno <- anno |>
    filter(gene_id == unique(gene))
  
  gene_info <- gene_anno |>
    filter(type == "gene") |>
    mutate(tss = if_else(strand == "+", start, end)) |>
    select(seqnames, tss, strand)
  
  bin_grange <- tibble(rel_pos = rel_posns) |>
    mutate(strand = gene_info$strand,
           pos = if_else(strand == "+",
                         gene_info$tss + rel_pos,
                         gene_info$tss - rel_pos)) |>
    with(
      GenomicRanges::GRanges(
        seqnames = gene_info$seqnames,
        ranges = IRanges::IRanges(start = pos, end = pos)
      )
    )
  
  exons_grange <- gene_anno |>
    filter(type == "exon") |>
    with(
      GenomicRanges::GRanges(
        seqnames = seqnames,
        ranges = IRanges::IRanges(start = start, end = end)
      )
    )
  
  GenomicRanges::countOverlaps(bin_grange, exons_grange, ignore.strand = TRUE) > 0
}

load_weights <- function(filename) {
  load(filename)
  i_best_model <- which.min(apply(cv.performance["pval", , drop=FALSE],2,min,na.rm=T))
  colnames(wgt.matrix)[i_best_model] <- "weight"
  wgt.matrix[, "weight", drop=FALSE] |>
    as_tibble(rownames = "variant_id") |>
    separate_wider_delim(variant_id, "_", names = c("chrom", "pos"), too_many = "drop") |>
    mutate(pos = as.integer(pos))
}

gene_id <- "ENSG00000170291"
PC <- "PC2"
trait <- "UKB_20022_Birth_weight"
window <- 0.3 # Mb up and downstream of TSS to show TWAS weights and GWAS sumstats

gene_names <- read_tsv("data/processed/protein_coding_genes.tsv", col_types = "cc----") |>
  deframe()
gene_name <- gene_names[gene_id]

## Overall coverage and phenotype-stratified coverage

anno <- rtracklayer::import("data/ref/Homo_sapiens.GRCh38.113.chr.gtf.gz") |>
  as_tibble()

covg <- read_tsv(
  str_glue("data/analyses/twas_example/model_info/inspect.{gene_id}.Geuvadis.tsv.gz"),
  col_types = cols(gene_id = "c", pos = "i", .default = "d")
) |>
  mutate(bin = seq_len(n()),
         in_exon = bin_in_exon(gene_id, pos, anno))

deciles <- covg |>
  select(bin,
         starts_with(str_glue("decile_{PC}_"))) |>
  pivot_longer(-bin, names_to = "decile", values_to = "median_coverage") |>
  mutate(decile = str_split_i(decile, "_", 3) |> fct_inorder(),
         decile_order = decile |>
           fct_relevel("1", "10", "2", "9", "3", "8", "4", "7", "5", "6"))

p1 <- covg |>
  filter(in_exon) |>
  ggplot(aes(x = bin, y = 1)) +
  annotate("segment", x = min(covg$bin[covg$in_exon]), xend = max(covg$bin[covg$in_exon]), y = 1, yend = 1) +
  geom_tile(fill = "black") +
  expand_limits(x = covg$bin) +
  theme_bw() +
  theme(
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    panel.grid = element_blank(),
    panel.border = element_blank(),
  ) +
  xlab("Exonic coverage bins (not to scale)") +
  ylab(NULL) +
  ggtitle(gene_name)

p2 <- covg |>
  ggplot(aes(x = bin, y = covg_50th, ymin = covg_25th, ymax = covg_75th)) +
  geom_ribbon(fill = "#ccc") +
  geom_line(aes(y = covg_25th), color = "#888") +
  geom_line(aes(y = covg_75th), color = "#888") +
  geom_line() +
  theme_bw() +
  theme(
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    panel.grid = element_blank(),
  ) +
  # xlab(str_glue("{gene_name} range (not to scale)")) +
  xlab(NULL) +
  ylab("Coverage") +
  ggtitle(str_glue("Median and interquartile {gene_name} coverage"))

p3 <- deciles |>
  ggplot(aes(x = bin, y = median_coverage, color = decile, group = decile_order)) +
  geom_line() +
  scale_color_manual(values = c("#208dff","#5a92ee","#7798dd","#8c9ecc","#9da4bb",
                                "#aaaaaa","#c9928a","#df776b","#f0554d","#fd1330"),
                     guide = guide_legend(override.aes = list(linewidth = 1))) +
  theme_bw() +
  theme(
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    legend.key.height = unit(8, "pt"),
    panel.grid = element_blank(),
  ) +
  # xlab(str_glue("{gene_name} range (not to scale)")) +
  xlab(NULL) +
  ylab("Median coverage for\nsamples in decile") +
  labs(color = "Phenotype\ndecile") +
  ggtitle(str_glue("{PC}-stratified {gene_name} coverage"))

p1 / p2 / p3 + plot_layout(heights = c(1, 10, 10))
ggsave("figures/figure1/figure1b_top.png", width = 6, height = 4, device = png)

## TWAS weights and GWAS sumstats

weights <- load_weights(str_glue("data/analyses/twas_example/weights/{gene_id}__{PC}.wgt.RDat"))

gwas <- read_tsv(
  str_glue("data/analyses/twas_example/gwas/{trait}-{gene_id}.tsv.gz"),
  col_types = cols(chromosome = "c", position = "i", zscore = "d", .default = "c")
)

gene_anno <- anno[anno$gene_id == gene_id, ]
tss <- gene_anno |>
  filter(type == "gene") |>
  with(if_else(strand == "+", start, end) / 1000000)

exons <- gene_anno |>
  filter(type == "exon") |>
  mutate(start = start / 1000000,
         end = end / 1000000)

p4 <- ggplot(exons, aes(xmin = start, xmax = end, ymin = 0, ymax = 2)) +
  annotate("segment", x = min(exons$start), xend = max(exons$end), y = 1, yend = 1) +
  geom_rect(fill = "black") +
  coord_cartesian(xlim = c(tss - window, tss + window), expand = 0) +
  theme_bw() +
  theme(
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    panel.grid = element_blank(),
    panel.border = element_blank(),
  ) +
  xlab(str_glue("{gene_name} location")) +
  ylab(NULL)

p5 <- weights |>
  mutate(pos = pos / 1000000) |>
  ggplot(aes(x = pos, y = weight)) +
  geom_vline(xintercept = tss, color = "#aaa") +
  geom_point(size = 1) +
  coord_cartesian(xlim = c(tss - window, tss + window)) +
  theme_bw() +
  theme(
    axis.text = element_text(color = "black"),
    panel.grid = element_blank(),
  ) +
  xlab(str_glue("{unique(weights$chrom)} position (Mb)")) +
  ylab("Weight") +
  ggtitle(str_glue("TWAS weights: {gene_name} {PC}"))

p6 <- gwas |>
  mutate(pos = position / 1000000) |>
  ggplot(aes(x = pos, y = zscore)) +
  geom_vline(xintercept = tss, color = "#aaa") +
  geom_point(size = 1) +
  coord_cartesian(xlim = c(tss - window, tss + window)) +
  theme_bw() +
  theme(
    axis.text = element_text(color = "black"),
    panel.grid = element_blank(),
  ) +
  xlab(str_glue("{unique(gwas$chromosome)} position (Mb)")) +
  ylab("Z-score") +
  ggtitle(str_glue("GWAS: {trait}"))

p4 / p5 / p6 + plot_layout(heights = c(1, 10, 10))
ggsave("figures/figure1/figure1b_bottom.png", width = 5, height = 4, device = png)
