# Calculate values in the paper text

library(tidyverse)

tissues_gtex <- read_lines("data/info/tissues.gtex.txt")

###

# "We trained RNA phenotype gene models using 28,637 samples from 54 human tissues (GTEx) and 33 human cancer types (TCGA)."

samples <- read_tsv(
  "data/info/coverage_manifest-gtextcga-full.tsv",
  col_names = c("tissue", "sample", "bigwig"),
  col_types = "ccc"
) |>
  mutate(dataset = str_split_i(bigwig, "/", 1))

nrow(samples)
samples |>
  distinct(dataset, tissue) |>
  count(dataset)

###

phenos <- tibble(tissue = tissues_gtex) |>
  reframe(
    read_tsv(
      str_glue("data/phenos/gtextcga-full/{tissue}-latent.phenotype_groups.txt.gz"),
      col_names = c("phenotype_id", "gene_id"),
      col_types = "cc"
    ),
    .by = tissue
  )

# "producing an average of 191,801 (SD 27.1) phenotypes per tissue across an average of 19,410 (SD 10.6) protein-coding genes."

phenos |>
  count(tissue) |>
  summarise(n_phenos_mean = mean(n),
            n_phenos_sd = sd(n))
phenos |>
  distinct(tissue, gene_id) |>
  count(tissue) |>
  summarise(n_phenos_mean = mean(n),
            n_phenos_sd = sd(n))

###

qtls_gtextcga_full <- read_tsv("data/processed/gtextcga-full.qtls.tsv.gz", col_types = "cciccd")

# "We found 4,265 to 48,100 conditionally independent xQTLs per tissue for 3,322 to 16,084 genes."

qtls_gtextcga_full |>
  count(tissue, sort = TRUE) |>
  slice(1:3, 47:n())

qtls_gtextcga_full |>
  distinct(tissue, gene_id) |>
  count(tissue, sort = TRUE) |>
  slice(1:3, 47:n())

# "finding that on average, 6.5% more cis-QTLs were found for latent RNA phenotypes than for explicit RNA phenotypes per tissue"

qtls_pantry <- read_tsv("data/pantry/processed/gtex.sep.qtls.tsv.gz",
                        col_types = "ccicccd")

full_join(
  qtls_gtextcga_full |>
    count(tissue),
  qtls_pantry |>
    count(tissue),
  by = "tissue",
  relationship = "one-to-one"
) |>
  mutate(pct_more = (n.x / n.y - 1) * 100) |>
  summarise(ave_pct_more = mean(pct_more))

# "with PC1 phenotypes producing 5.6 times as many xQTLs as PC8 phenotypes and 13.5 times as many xQTLs as PC16 phenotypes"

qtls_gtextcga_full |>
  separate_wider_delim(phenotype_id, "__", names = c("gene", "PC")) |>
  count(PC) |>
  summarise(PC1_PC8 = n[PC == "PC1"] / n[PC == "PC8"],
            PC1_PC16 = n[PC == "PC1"] / n[PC == "PC16"])

###

# "and examined the maximum Pearson r2 between each residual latent phenotype and an explicit phenotype of the same gene. For PC1 latent phenotypes, residualization reduced these values from mean 0.49 (SD 0.24) to mean 0.089 (SD 0.088) (Figure XXa)."

corrs_max <- read_tsv("data/processed/latent_explicit_corrs.tsv.gz", col_types = "ccccd") |>
  mutate(r2_max = r^2)

corrs_max |>
  filter(PC == "PC1") |>
  summarise(r2_max_mean = mean(r2_max),
            r2_max_sd = sd(r2_max),
            .by = latent)

# "For all latent phenotypes, residualization reduced the max r2 values from mean 0.089 (SD 0.15) to mean 0.036 (SD 0.053)."

corrs_max |>
  summarise(r2_max_mean = mean(r2_max),
            r2_max_sd = sd(r2_max),
            .by = latent)

###

# "Explicit modalities resulted in 24,170 total xQTLs, while full latent phenotypes resulted in 27,521 xQTLs.
# Adding residual latent phenotypes to explicit phenotypes resulted in a net increase of 7,080 xQTLs for a total of 31,250, 20,140 of which were called for explicit phenotypes and 11,110 of which were called for residual latent phenotypes"

qtls_geuvadis <- read_tsv("data/processed/geuvadis.qtls.tsv.gz", col_types = "ccccdci")

qtls_geuvadis |>
  count(version)

qtls_geuvadis |>
  filter(version == "residual-cross_latent") |>
  mutate(is_latent = modality == "latent_residual") |>
  count(is_latent)

# "In terms of unique genes represented by the xQTLs, addition of the latent residual phenotypes increased the number of xGenes from 13,926 to 16,176"

qtls_geuvadis |>
  distinct(version, gene_id) |>
  count(version)

###

# "Applying these to GWAS data for a collection of 114 complex traits resulted in a median of 21,178 significant TWAS associations per tissue, including a total of 18,867 unique gene-trait pairs, 7,658 of which include hits with strong evidence of colocalization at the level of shared causal variant"
# TODO update to use all GTEX tissues

twas_gtex5 <- read_tsv("data/processed/gtex5.twas_hits.tsv.gz", col_types = "ccccdddd")
twas_pantry <- read_tsv("data/pantry/processed/gtex.twas_hits.tsv.gz", col_types = "cccc----d----d")

twas_gtex5 |>
  count(tissue) |>
  summarise(n_median = median(n))

twas_gtex5 |>
  distinct(gene_id, trait) |>
  count()

twas_gtex5 |>
  filter(coloc_pp > 0.8) |>
  distinct(gene_id, trait) |>
  count()

# "Compared to six-modality xTWAS and using the same p-value threshold, latent phenotypes resulted in an average of 54% more associations per tissue, 35% more unique gene-trait pairs overall with associations, and 31% more unique gene-trait pairs overall with strongly colocalizing associations."

left_join(
  twas_gtex5 |>
    filter(twas_p < 5e-8 / 6) |>
    count(tissue, name = "n_latent"),
  twas_pantry |>
    count(tissue, name = "n_explicit"),
  by = "tissue"
) |>
  mutate(percent_inc = (n_latent / n_explicit - 1) * 100) |>
  summarise(mean_percent_inc = mean(percent_inc))

bind_cols(
  twas_gtex5 |>
    filter(twas_p < 5e-8 / 6,
           tissue %in% twas_pantry$tissue) |>
    distinct(gene_id, trait) |>
    count(name = "n_latent"),
  twas_pantry |>
    filter(tissue %in% twas_gtex5$tissue) |>
    distinct(gene_id, trait) |>
    count(name = "n_explicit"),
) |>
  mutate(percent_inc = (n_latent / n_explicit - 1) * 100)

bind_cols(
  twas_gtex5 |>
    filter(twas_p < 5e-8 / 6,
           tissue %in% twas_pantry$tissue,
           coloc_pp > 0.8) |>
    distinct(gene_id, trait) |>
    count(name = "n_latent"),
  twas_pantry |>
    filter(tissue %in% twas_gtex5$tissue,
           COLOC.PP4 > 0.8) |>
    distinct(gene_id, trait) |>
    count(name = "n_explicit"),
) |>
  mutate(percent_inc = (n_latent / n_explicit - 1) * 100)

###

qtls_gtex5_full <- read_tsv("data/processed/gtex5-full.qtls.tsv.gz", col_types = "cciccd")

qtls_gtex_full <- read_tsv("data/processed/gtex-full.qtls.tsv.gz", col_types = "cciccd")

qtl_counts <- full_join(
  qtls_gtex5_full |>
    count(tissue, name = "n_gtex5"),
  qtls_gtex_full |>
    count(tissue, name = "n_gtex"),
  by = "tissue",
  relationship = "one-to-one"
) |>
  full_join(
    qtls_gtextcga_full |>
      count(tissue, name = "n_gtextcga"),
    by = "tissue",
    relationship = "one-to-one"
  ) |>
  mutate(
    gtex5_gtex_pct_inc = (n_gtex / n_gtex5 - 1) * 100,
    gtex_gtextcga_pct_inc = (n_gtextcga / n_gtex - 1) * 100,
  )

# Surprisingly, using 54 tissues to fit models increased the number of independent xQTLs per tissue by only 2.3% on average (Figure 4a).
# Using 54 tissues plus 33 cancer datasets increased xQTLs by 3.0% on average compared to 54 tissues only (Figure 4b).

qtl_counts |>
  summarise(
    gtex5_gtex_pct_inc_mean = mean(gtex5_gtex_pct_inc),
    gtex_gtextcga_pct_inc_mean = mean(gtex_gtextcga_pct_inc),
  )
