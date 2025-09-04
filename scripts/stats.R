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

# "producing an average of 417,593 (SD 9,274) phenotypes per tissue across an average of 40,147 (SD 1,351) protein-coding genes."

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

# "We found 3,755 to 64,106 conditionally independent xQTLs per tissue for 3,419 to 29,007 genes."

qtls_gtextcga_full |>
  count(tissue, sort = TRUE) |>
  slice(1:3, 47:n())

qtls_gtextcga_full |>
  distinct(tissue, gene_id) |>
  count(tissue, sort = TRUE) |>
  slice(1:3, 47:n())

# "finding that on average, 95% more independent cis-QTLs were found for latent RNA phenotypes than for explicit RNA phenotypes per tissue"

qtls_pantry <- read_tsv("data/processed/gtex-pantry.qtls.tsv.gz", col_types = "ccicccd")

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

# "with PC1 phenotypes producing 7.4 times as many xQTLs as PC8 phenotypes and 18.6 times as many xQTLs as PC16 phenotypes"

qtls_gtextcga_full |>
  separate_wider_delim(phenotype_id, "__", names = c("gene", "PC")) |>
  count(PC) |>
  summarise(PC1_PC8 = n[PC == "PC1"] / n[PC == "PC8"],
            PC1_PC16 = n[PC == "PC1"] / n[PC == "PC16"])

###

# "Compared to six-modality xTWAS at the same p-value threshold of 5⨉10-8 and using the same Geuvadis dataset for transcriptomic models, latent phenotypes resulted in nearly the same number of total associations (24,697 vs. 24,644 for six-modality), 33% more unique gene-trait pairs with associations, and 37% more unique gene-trait pairs with strong evidence of colocalization at the level of shared causal variant."

twas_geuv <- read_tsv("data/processed/geuvadis-full.twas_hits.tsv.gz", col_types = "cccdddd")
twas_geuv_pantry <- read_tsv("data/processed/geuvadis-pantry.twas_hits.tsv.gz", col_types = "ccccdddd")
twas_geuv_both <- bind_rows(
  twas_geuv |> mutate(type = "latent"),
  twas_geuv_pantry |> mutate(type = "explicit"),
)

bind_cols(
  twas_geuv |>
    count(name = "n_latent"),
  twas_geuv_pantry |>
    count(name = "n_explicit"),
) |>
  mutate(percent_inc = (n_latent / n_explicit - 1) * 100)

bind_cols(
  twas_geuv |>
    distinct(gene_id, trait) |>
    count(name = "n_latent"),
  twas_geuv_pantry |>
    distinct(gene_id, trait) |>
    count(name = "n_explicit"),
) |>
  mutate(percent_inc = (n_latent / n_explicit - 1) * 100)

bind_cols(
  twas_geuv |>
    filter(coloc_pp > 0.8) |>
    distinct(gene_id, trait) |>
    count(name = "n_latent"),
  twas_geuv_pantry |>
    filter(coloc_pp > 0.8) |>
    distinct(gene_id, trait) |>
    count(name = "n_explicit"),
) |>
  mutate(percent_inc = (n_latent / n_explicit - 1) * 100)

# "Applying xTWAS to latent RNA phenotypes for 49 GTEx tissues and the same 114 traits resulted in a median of 26,501 significant TWAS associations per GTEx tissue, including a total of 64,314 unique gene-trait pairs, 28,116 of which include strongly colocalizing associations"

twas_gtex <- read_tsv("data/processed/gtextcga-full.twas_hits.tsv.gz", col_types = "ccccdddd")

twas_gtex |>
  count(tissue) |>
  summarise(n_median = median(n))

twas_gtex |>
  distinct(gene_id, trait) |>
  count()

twas_gtex |>
  filter(coloc_pp > 0.8) |>
  distinct(gene_id, trait) |>
  count()

###

# "Despite long non-coding RNAs (lncRNAs) having higher transcript complexity (i.e., splice variants per exon) than protein-coding mRNAs, xQTL mapping across six explicit modalities produced an average of 0.16 independent xQTLs per gene per tissue for lncRNAs, which is only 29% of the rate for protein-coding genes."

genes <- read_tsv("data/processed/pcg_and_lncrna.tsv", col_types = "c-c----")

genes |>
  left_join(count(qtls_pantry, gene_id), by = "gene_id", relationship = "one-to-one") |>
  replace_na(list(n = 0)) |>
  mutate(n_per_tissue = n / n_distinct(qtls_pantry$tissue)) |>
  summarise(mean_n_per_tissue = mean(n_per_tissue),
            .by = gene_biotype)

# "it is similar to the rate of total annotated exons per gene (i.e., the product of the rate of exons per isoform and the rate of isoforms per gene), for which lncRNAs have 28.9% the rate of protein-coding genes."

anno <- rtracklayer::import("data/ref/Homo_sapiens.GRCh38.113.chr.gtf.gz") |>
  as_tibble() |>
  filter(gene_id %in% genes$gene_id)

anno |>
  filter(type == "exon") |>
  count(gene_id, gene_biotype) |>
  summarise(mean_exons_per_gene = mean(n),
            .by = gene_biotype)

# "While latent RNA phenotype xQTL mapping also resulted in a lower rate of independent xQTLs for lncRNAs than for protein-coding genes, there were 0.39 xQTLs per gene per tissue, which was 44% of the rate for protein-coding genes."

genes |>
  left_join(count(qtls_gtextcga_full, gene_id), by = "gene_id", relationship = "one-to-one") |>
  replace_na(list(n = 0)) |>
  mutate(n_per_tissue = n / n_distinct(qtls_gtextcga_full$tissue)) |>
  summarise(mean_n_per_tissue = mean(n_per_tissue),
            .by = gene_biotype)

###

# "35.6% of latent phenotype TWAS gene-trait pairs involving protein-coding genes had a nearby gene associated with the same trait, compared to 31.4% of explicit phenotype TWAS pairs."

# "For gene-trait pairs involving lncRNAs, the percentages were higher, at 59.9% for latent phenotype TWAS and 57.7% for explicit phenotype TWAS."

# See figureS1.R

###

# "and examined the maximum Pearson r2 between each residual latent phenotype and an explicit phenotype of the same gene. For PC1 latent phenotypes, residualization reduced these values from mean 0.40 (SD 0.27) to mean 0.068 (SD 0.093) (Figure 4a)."

corrs_max <- read_tsv("data/processed/latent_explicit_corrs.tsv.gz", col_types = "ccccd") |>
  mutate(r2_max = r^2)

corrs_max |>
  filter(PC == "PC1") |>
  summarise(r2_max_mean = mean(r2_max),
            r2_max_sd = sd(r2_max),
            .by = latent)

# "For all latent phenotypes, residualization reduced the max r2 values from mean 0.072 (SD 0.14) to mean 0.030 (SD 0.055)."

corrs_max |>
  summarise(r2_max_mean = mean(r2_max),
            r2_max_sd = sd(r2_max),
            .by = latent)

###

# "Explicit modalities resulted in 22,470 total xQTLs, while full latent phenotypes resulted in 30,158 xQTLs.
# Adding residual latent phenotypes to explicit phenotypes resulted in a net increase of 9,408 xQTLs (42%) for a total of 31,878, 16,544 of which were called for explicit phenotypes and 15,334 of which were called for residual latent phenotypes"

qtls_geuvadis <- read_tsv("data/processed/geuvadis.qtls.tsv.gz", col_types = "ccccdci")

qtls_geuvadis |>
  count(version)

qtls_geuvadis |>
  filter(version == "residual-cross_latent") |>
  mutate(is_latent = modality == "latent_residual") |>
  count(is_latent)

# "In terms of unique genes represented by the xQTLs, addition of the latent residual phenotypes increased the number of xGenes from 12,872 to 16,950"

qtls_geuvadis |>
  distinct(version, gene_id) |>
  count(version)

###

#  In every case, holding out one modality increased the number of latent residual phenotypes, resulting in a similar total number of independent cis-QTLs

qtls_held_out <- read_tsv(
  "data/processed/held_out-geuvadis.qtls.tsv.gz", col_types = "ccccdci"
)

qtls_held_out |>
  filter(modality == "latent_residual") |>
  count(held_out, sort = TRUE)

###

twas_pantry <- read_tsv("data/processed/geuvadis-pantry.twas_hits.tsv.gz", col_types = "ccc---d-")

twas_resid <- read_tsv("data/processed/geuvadis-residual.twas_hits.tsv.gz", col_types = "cc---d-") |>
  mutate(modality = "latent", .before = 1)

# "This resulted in a 37% increase in unique gene-trait association pairs"

bind_rows(twas_pantry, twas_resid) |>
  summarise(latent_only = all(modality == "latent"),
            .by = c(trait, gene_id)) |>
  with(mean(latent_only))

# "and residual latent phenotypes had the strongest associations for 7,655 (53%) of all unique pairs"

bind_rows(twas_pantry, twas_resid) |>
  slice_min(twas_p, n = 1, with_ties = FALSE, by = c(trait, gene_id)) |>
  summarise(n_latent = sum(modality == "latent"),
            frac_latent = mean(modality == "latent"))

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

# Surprisingly, using 54 tissues to fit models increased the number of independent xQTLs per tissue by only 4.2% on average (Figure 4a).
# Using 54 tissues plus 33 cancer datasets increased xQTLs by 1.7% on average compared to 54 tissues only (Figure 4b).

qtl_counts |>
  summarise(
    gtex5_gtex_pct_inc_mean = mean(gtex5_gtex_pct_inc),
    gtex_gtextcga_pct_inc_mean = mean(gtex_gtextcga_pct_inc),
  )

###

qtls_prune <- read_tsv("data/processed/prune-BRNCTXB.qtls.tsv.gz", col_types = "cicciccd") |>
  filter(map_group %in% c("latent", "pantry"))

# "Latent RNA phenotypes were less affected by gene annotation sparsity than were explicit phenotypes, with independent xQTLs dropping X% between 100% and 0% inclusion of non-canonical isoforms, compared to X% for explicit phenotypes"

qtls_prune |>
  count(map_group, pruning) |>
  summarise((n[pruning == 0] / n[pruning == 100] - 1) * 100,
            .by = map_group)

###

qtl_counts_seqsim <- read_tsv("data/processed/seqsim.qtls.tsv.gz", col_types = "cccdci") |>
  count(reads) |>
  deframe()

# "When we used each of these variations as input for latent RNA phenotyping and mapped cis-QTLs, we observed a relatively small drop (-6.4%) in discoveries from the 50 bp truncated reads, and larger drops from the single-end (-15%) and lowered sequencing depth (-49% for 50% of reads, -74% for 25% of reads) simulations"

(qtl_counts_seqsim["pe-50bp-100pct"] / qtl_counts_seqsim["pe-75bp-100pct"] - 1) * 100
(qtl_counts_seqsim["se-75bp-100pct"] / qtl_counts_seqsim["pe-75bp-100pct"] - 1) * 100
(qtl_counts_seqsim["pe-75bp-50pct"] / qtl_counts_seqsim["pe-75bp-100pct"] - 1) * 100
(qtl_counts_seqsim["pe-75bp-25pct"] / qtl_counts_seqsim["pe-75bp-100pct"] - 1) * 100
