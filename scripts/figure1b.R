## Example latent RNA phenotype and TWAS association

library(tidyverse)
library(patchwork)

bin_in_exon_iso <- function(bin_grange, seqnames, starts, ends) {
  exons_grange <- GenomicRanges::GRanges(
    seqnames = seqnames,
    ranges = IRanges::IRanges(start = starts, end = ends)
  )
  
  GenomicRanges::countOverlaps(bin_grange, exons_grange, ignore.strand = TRUE) > 0
}

bins_in_exons_iso <- function(gene, rel_posns, anno) {
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
  
  gene_anno |>
    filter(type == "exon") |>
    reframe(
      tibble(
        bin = seq_len(length(bin_grange)),
        in_exon = bin_in_exon_iso(bin_grange, seqnames, start, end),
      ),
      .by = transcript_id
    )
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

## ELP5
# gene_id <- "ENSG00000170291"
# PC <- "PC2"
# trait <- "UKB_20022_Birth_weight"
# window <- 0.3 # Mb up and downstream of TSS to show TWAS weights and GWAS sumstats
# tissue <- "Geuvadis"
## YPEL3
# gene_id <- "ENSG00000090238"
# PC <- "PC4"
# window <- 0.5 # Mb up and downstream of TSS to show TWAS weights and GWAS sumstats
# trait <- "Astle_et_al_2016_Eosinophil_counts"
# tissue <- "Geuvadis"
## IL7
# gene_id <- "ENSG00000104432"
# PC <- "PC2"
# window <- 0.5 # Mb up and downstream of TSS to show TWAS weights and GWAS sumstats
# trait <- "Astle_et_al_2016_Lymphocyte_counts"
# tissue <- "SPLEEN"
## H2AZ2
# gene_id <- "ENSG00000105968"
# PC <- "PC3"
# window <- 0.5 # Mb up and downstream of TSS to show TWAS weights and GWAS sumstats
# trait <- "Astle_et_al_2016_Red_blood_cell_count"
# tissue <- "WHLBLD"
## ZNF844
gene_id <- "ENSG00000223547"
PC <- "PC2"
window <- 0.5 # Mb up and downstream of TSS to show TWAS weights and GWAS sumstats
trait <- "UKB_50_Standing_height"
tissue <- "NERVET"

DP <-str_replace(PC, "PC", "DP")
gene_names <- read_tsv("data/processed/pcg_and_lncrna.tsv", col_types = "cc-----") |>
  deframe()
gene_name <- gene_names[gene_id]
trait_names <- read_tsv("data/pantry/geuvadis/twas/gwas_metadata.txt",
                        col_types = cols(Tag = "c", Phenotype = "c", .default = "-")) |>
  deframe()
trait_name <- trait_names[trait]

## Overall coverage and phenotype-stratified coverage

anno <- rtracklayer::import("data/ref/Homo_sapiens.GRCh38.113.chr.gtf.gz") |>
  as_tibble()

covg <- read_tsv(
  str_glue("data/analyses/twas_example/model_info/inspect.{gene_id}.{tissue}.tsv.gz"),
  col_types = cols(gene_id = "c", pos = "i", .default = "d")
) |>
  mutate(bin = seq_len(n()))

iso_bins <- bins_in_exons_iso(gene_id, covg$pos, anno)

iso_ranges <- iso_bins |>
  filter(in_exon) |>
  summarise(start = min(bin),
            end = max(bin),
            .by = transcript_id)

deciles <- covg |>
  select(bin,
         starts_with(str_glue("decile_{PC}_"))) |>
  pivot_longer(-bin, names_to = "decile", values_to = "median_coverage") |>
  mutate(decile = str_split_i(decile, "_", 3) |> fct_inorder(),
         decile_order = decile |>
           fct_relevel("1", "10", "2", "9", "3", "8", "4", "7", "5", "6"))

p1 <- iso_bins |>
  filter(in_exon) |>
  ggplot(aes(xmin = bin, xmax = bin + 1.2, y = transcript_id)) +
  geom_linerange(aes(xmin = start, xmax = end), data = iso_ranges, linewidth = 0.5) +
  geom_linerange(linewidth = 3) +
  expand_limits(x = covg$bin) +
  theme_bw() +
  theme(
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    panel.grid = element_blank(),
    panel.border = element_blank(),
  ) +
  xlab("Exons aligned to coverage bins (not to scale)") +
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
  ggtitle(str_glue("Median and interquartile {gene_name} coverage in {tissue}"))

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
  xlab("Bins along gene start to end") +
  ylab("Median coverage for\nsamples in decile") +
  labs(color = "Phenotype\ndecile") +
  ggtitle(str_glue("{DP}-stratified {gene_name} coverage in {tissue}"))

p1 / p2 / p3 + plot_layout(heights = c(1, 2, 2))
ggsave("figures/figure1/figure1b_top.png", width = 6, height = 4.5, device = png)

## TWAS weights and GWAS sumstats

weights <- load_weights(str_glue("data/analyses/twas_example/weights/{gene_id}__{PC}.{tissue}.wgt.RDat"))

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
  mutate(pos = pos / 1e6) |>
  ggplot(aes(x = pos, y = weight)) +
  geom_vline(xintercept = tss, color = "#aaa") +
  geom_point(size = 1) +
  coord_cartesian(xlim = c(tss - window, tss + window),
                  ylim = c(min(weights$weight) - 0.6, max(weights$weight) + 0.6),
                  expand = 0) +
  theme_bw() +
  theme(
    axis.text = element_text(color = "black"),
    panel.grid = element_blank(),
  ) +
  xlab(str_glue("{unique(weights$chrom)} position (Mb)")) +
  ylab("Weight") +
  ggtitle(str_glue("TWAS weights: {gene_name} {DP} in {tissue}"))

p6 <- gwas |>
  mutate(pos = position / 1000000) |>
  ggplot(aes(x = pos, y = zscore)) +
  geom_vline(xintercept = tss, color = "#aaa") +
  geom_point(size = 1) +
  coord_cartesian(xlim = c(tss - window, tss + window),
                  ylim = c(min(gwas$zscore) - 1, max(gwas$zscore) + 1),
                  expand = 0) +
  theme_bw() +
  theme(
    axis.text = element_text(color = "black"),
    panel.grid = element_blank(),
  ) +
  xlab(str_glue("{unique(gwas$chromosome)} position (Mb)")) +
  ylab("Z-score") +
  ggtitle(str_glue("GWAS: {trait_name}"))

p4 / p5 / p6 + plot_layout(heights = c(1, 10, 10))
ggsave("figures/figure1/figure1b_bottom.png", width = 5, height = 4.5, device = png)
