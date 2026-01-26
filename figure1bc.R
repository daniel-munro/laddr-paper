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

DDP <-str_replace(PC, "PC", "DDP")
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

deciles1 <- covg |>
  select(bin,
         starts_with(str_glue("decile_PC1_"))) |>
  pivot_longer(-bin, names_to = "decile", values_to = "median_coverage") |>
  mutate(decile = str_split_i(decile, "_", 3) |> fct_inorder(),
         decile_order = decile |>
           fct_relevel("1", "10", "2", "9", "3", "8", "4", "7", "5", "6"))

deciles2 <- covg |>
  select(bin,
         starts_with(str_glue("decile_PC2_"))) |>
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
  xlab("Isoforms scaled to align with coverage bins") +
  ylab(NULL) +
  ggtitle(gene_name)

p2 <- covg |>
  ggplot(aes(x = bin, y = covg_50th, ymin = covg_25th, ymax = covg_75th)) +
  geom_ribbon(fill = "#ccc") +
  geom_line(aes(y = covg_25th), color = "#888") +
  geom_line(aes(y = covg_75th), color = "#888") +
  geom_line() +
  scale_y_continuous(breaks = c(0, 20, 40)) +
  theme_bw() +
  theme(
    axis.text.x = element_blank(),
    axis.text.y = element_text(color = "black"),
    axis.ticks.x = element_blank(),
    panel.grid = element_blank(),
  ) +
  xlab(NULL) +
  ylab("Coverage") +
  ggtitle("RNA-seq coverage")

p3 <- deciles1 |>
  ggplot(aes(x = bin, y = median_coverage, color = decile, group = decile_order)) +
  geom_line(key_glyph = "rect") +
  scale_y_continuous(breaks = c(0, 20, 40)) +
  scale_color_manual(values = c("#208dff","#5a92ee","#7798dd","#8c9ecc","#9da4bb",
                                "#aaaaaa","#c9928a","#df776b","#f0554d","#fd1330")) +
  theme_bw() +
  theme(
    axis.text.x = element_blank(),
    axis.text.y = element_text(color = "black"),
    axis.ticks.x = element_blank(),
    legend.key.height = unit(6, "pt"),
    legend.key.width = unit(8, "pt"),
    panel.grid = element_blank(),
  ) +
  xlab(NULL) +
  ylab("Coverage") +
  labs(color = NULL) +
  ggtitle(str_glue("DDP1-stratified coverage"))

p4 <- deciles2 |>
  ggplot(aes(x = bin, y = median_coverage, color = decile, group = decile_order)) +
  geom_line(key_glyph = "rect") +
  scale_y_continuous(breaks = c(0, 20, 40)) +
  scale_color_manual(values = c("#208dff","#5a92ee","#7798dd","#8c9ecc","#9da4bb",
                                "#aaaaaa","#c9928a","#df776b","#f0554d","#fd1330")) +
  theme_bw() +
  theme(
    axis.text.x = element_blank(),
    axis.text.y = element_text(color = "black"),
    axis.ticks.x = element_blank(),
    legend.key.height = unit(6, "pt"),
    legend.key.width = unit(8, "pt"),
    panel.grid = element_blank(),
  ) +
  xlab("Bins along gene start to end") +
  ylab("Coverage") +
  labs(color = NULL) +
  ggtitle(str_glue("DDP2-stratified coverage"))

## TWAS weights and GWAS sumstats

xqtl <- read_tsv(
  str_glue("data/analyses/twas_example/cis_nominal/{gene_id}__{PC}.{tissue}.tsv.gz"),
  col_types = cols(variant_id = "c", pval_nominal = "d", .default = "-")
) |>
  separate_wider_delim(variant_id, "_", names = c("chrom", "pos", "ref", "alt", "b38")) |>
  mutate(pos = as.integer(pos))

gwas <- read_tsv(
  str_glue("data/analyses/twas_example/gwas/{trait}-{gene_id}.tsv.gz"),
  col_types = cols(chromosome = "c", position = "i", zscore = "d", pvalue = "d", .default = "c")
)

gene_anno <- anno[anno$gene_id == gene_id, ]
tss <- gene_anno |>
  filter(type == "gene") |>
  with(if_else(strand == "+", start, end) / 1000000)

exons <- gene_anno |>
  filter(type == "exon") |>
  mutate(start = start / 1000000,
         end = end / 1000000)

p5 <- ggplot(exons, aes(xmin = start, xmax = end + 2e-3, ymin = 0, ymax = 2)) +
  annotate("segment", x = min(exons$start), xend = max(exons$end), y = 1, yend = 1) +
  geom_rect(fill = "black") +
  coord_cartesian(xlim = c(tss - window, tss + window), expand = 0) +
  theme_bw() +
  theme(
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    # margins = margin_part(t = 0, r = 0, l = 0),
    margins = margin_part(t = 20, b = 0),
    panel.grid = element_blank(),
    panel.border = element_blank(),
  ) +
  xlab(str_glue("{gene_name} location")) +
  ylab(NULL)

p6 <- xqtl |>
  mutate(pos = pos / 1e6) |>
  ggplot(aes(x = pos, y = -log10(pval_nominal))) +
  geom_vline(xintercept = tss, color = "#aaa") +
  geom_point(size = 0.5) +
  coord_cartesian(xlim = c(tss - window, tss + window),
                  ylim = c(0, max(-log10(xqtl$pval_nominal)) + 2),
                  expand = 0) +
  scale_y_continuous(breaks = c(0, 5, 10)) +
  theme_bw() +
  theme(
    axis.text = element_text(color = "black"),
    # margins = margin_part(r = 0, l = 0),
    panel.grid = element_blank(),
  ) +
  xlab(str_glue("{unique(xqtl$chrom)} position (Mb)")) +
  ylab(expression(-log[10]*" P")) +
  ggtitle(str_glue("xQTL: {gene_name} {DDP}"))

p7 <- gwas |>
  mutate(pos = position / 1000000) |>
  ggplot(aes(x = pos, y = -log10(pvalue))) +
  geom_vline(xintercept = tss, color = "#aaa") +
  geom_point(size = 0.5) +
  coord_cartesian(xlim = c(tss - window, tss + window),
                  ylim = c(0, max(-log10(gwas$pvalue)) + 2),
                  expand = 0) +
  scale_y_continuous(breaks = c(0, 10, 20)) +
  theme_bw() +
  theme(
    axis.text = element_text(color = "black"),
    # margins = margin_part(t = 20, r = 0, b = 0, l = 0),
    panel.grid = element_blank(),
  ) +
  xlab(str_glue("{unique(gwas$chromosome)} position (Mb)")) +
  ylab(expression(-log[10]*" P")) +
  ggtitle(str_glue("GWAS: {trait_name}"))

p1 / p2 / p3 / p4 / p5 / p6 / p7 + plot_layout(heights = c(1, 2, 2, 2, 0.3, 2, 2))
ggsave("figures/figure1/figure1bc.png", width = 5, height = 10, device = png)
