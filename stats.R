# Calculate values in the paper text

library(tidyverse)

tissues_gtex <- read_lines("data/info/tissues.gtex.txt")

###

qtls_gtex_ddp <- read_tsv("data/processed/gtextcga-full.qtls.tsv.gz", col_types = "cciccd")

qtls_gtex_kdp <- read_tsv("data/processed/gtex-pantry.qtls.tsv.gz", col_types = "ccicccd")

qtls_gtex_rddp <- read_tsv("data/processed/gtex-residual-cross.qtls.tsv.gz", col_types = "ccicccd")

qtl_counts_gtex <- full_join(
  qtls_gtex_ddp |>
    count(tissue, name = "n_ddp"),
  qtls_gtex_kdp |>
    count(tissue, name = "n_kdp"),
  by = "tissue",
  relationship = "one-to-one"
) |>
  full_join(
    qtls_gtex_rddp |>
      count(tissue, name = "n_rddp"),
    by = "tissue",
    relationship = "one-to-one"
  )

# "Applied to GTEx, LaDDR identified 95% more independent xQTLs per tissue on average than the six transcriptional regulation modes implemented in Pantry."

qtl_counts_gtex |>
  mutate(pct_more = ((n_ddp - n_kdp) / n_kdp) * 100) |>
  summarise(ave_pct_more = mean(pct_more))

# "Residualizing known modalities prior to LaDDR increased discovery by an additional 81% per tissue on average"

qtl_counts_gtex |>
  mutate(pct_more = ((n_rddp - n_ddp) / n_kdp) * 100) |>
  summarise(ave_pct_more = mean(pct_more))

###

twas_gtex_ddp <- read_tsv("data/processed/gtextcga-full.twas_hits.tsv.gz", col_types = "ccccdddd")

twas_gtex_kdp <- read_tsv("data/processed/gtex-pantry.twas_hits.tsv.gz", col_types = "cccc---dd")

# "LaDDR uncovered 11,790 unique gene–trait pairs with significant associations per tissue on average, versus 8,579 from knowledge-driven phenotypes."

twas_gtex_ddp |>
  distinct(tissue, gene_id, trait) |>
  count(tissue) |>
  summarise(mean_pairs = mean(n))

twas_gtex_kdp |>
  distinct(tissue, gene_id, trait) |>
  count(tissue) |>
  summarise(mean_pairs = mean(n))

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

phenos_gtex_ddp <- tibble(tissue = tissues_gtex) |>
  reframe(
    read_tsv(
      str_glue("data/phenos/gtextcga-full/{tissue}-latent.phenotype_groups.txt.gz"),
      col_names = c("phenotype_id", "gene_id"),
      col_types = "cc"
    ),
    .by = tissue
  )

# "producing an average of 417,593 (SD 9,274) phenotypes per tissue across an average of 40,147 (SD 1,351) protein-coding genes."

phenos_gtex_ddp |>
  count(tissue) |>
  summarise(n_phenos_mean = mean(n),
            n_phenos_sd = sd(n))
phenos_gtex_ddp |>
  distinct(tissue, gene_id) |>
  count(tissue) |>
  summarise(n_phenos_mean = mean(n),
            n_phenos_sd = sd(n))

###

# "We found 3,755 to 64,106 conditionally independent xQTLs per tissue for 3,419 to 29,007 genes."

qtls_gtex_ddp |>
  count(tissue, sort = TRUE) |>
  slice(1:3, 47:n())

qtls_gtex_ddp |>
  distinct(tissue, gene_id) |>
  count(tissue, sort = TRUE) |>
  slice(1:3, 47:n())

# "finding that on average, 95% more independent cis-QTLs were found for latent RNA phenotypes than for explicit RNA phenotypes per tissue"

full_join(
  qtls_gtex_ddp |>
    count(tissue),
  qtls_gtex_kdp |>
    count(tissue),
  by = "tissue",
  relationship = "one-to-one"
) |>
  mutate(pct_more = (n.x / n.y - 1) * 100) |>
  summarise(ave_pct_more = mean(pct_more))

# "with PC1 phenotypes producing 7.4 times as many xQTLs as PC8 phenotypes and 18.6 times as many xQTLs as PC16 phenotypes"

qtls_gtex_ddp |>
  separate_wider_delim(phenotype_id, "__", names = c("gene", "PC")) |>
  count(PC) |>
  summarise(PC1_PC8 = n[PC == "PC1"] / n[PC == "PC8"],
            PC1_PC16 = n[PC == "PC1"] / n[PC == "PC16"])

###

# "Compared to six-modality xTWAS at the same p-value threshold of 5⨉10-8 and using the same Geuvadis dataset for transcriptomic models, latent phenotypes resulted in nearly the same number of total associations (24,697 vs. 24,644 for six-modality), 33% more unique gene-trait pairs with associations, and 37% more unique gene-trait pairs with strong evidence of colocalization at the level of shared causal variant."

twas_geuv_ddp <- read_tsv("data/processed/geuvadis-full.twas_hits.tsv.gz", col_types = "cccdddd")
twas_geuv_kdp <- read_tsv("data/processed/geuvadis-pantry.twas_hits.tsv.gz", col_types = "ccccdddd")

bind_cols(
  twas_geuv_ddp |>
    count(name = "n_latent"),
  twas_geuv_kdp |>
    count(name = "n_explicit"),
) |>
  mutate(percent_inc = (n_latent / n_explicit - 1) * 100)

bind_cols(
  twas_geuv_ddp |>
    distinct(gene_id, trait) |>
    count(name = "n_latent"),
  twas_geuv_kdp |>
    distinct(gene_id, trait) |>
    count(name = "n_explicit"),
) |>
  mutate(percent_inc = (n_latent / n_explicit - 1) * 100)

bind_cols(
  twas_geuv_ddp |>
    filter(coloc_pp > 0.8) |>
    distinct(gene_id, trait) |>
    count(name = "n_latent"),
  twas_geuv_kdp |>
    filter(coloc_pp > 0.8) |>
    distinct(gene_id, trait) |>
    count(name = "n_explicit"),
) |>
  mutate(percent_inc = (n_latent / n_explicit - 1) * 100)

# "Applying xTWAS to latent RNA phenotypes for 49 GTEx tissues and the same 114 traits resulted in a median of 26,501 significant TWAS associations per GTEx tissue, including a total of 64,314 unique gene-trait pairs, 28,116 of which include strongly colocalizing associations"

twas_gtex_ddp |>
  count(tissue) |>
  summarise(n_median = median(n))

twas_gtex_ddp |>
  distinct(gene_id, trait) |>
  count()

twas_gtex_ddp |>
  filter(coloc_pp > 0.8) |>
  distinct(gene_id, trait) |>
  count()

###

# "Despite long non-coding RNAs (lncRNAs) having higher transcript complexity (i.e., splice variants per exon) than protein-coding mRNAs, xQTL mapping across six explicit modalities produced an average of 0.16 independent xQTLs per gene per tissue for lncRNAs, which is only 29% of the rate for protein-coding genes."

genes <- read_tsv("data/processed/pcg_and_lncrna.tsv", col_types = "ccc----")

genes |>
  left_join(count(qtls_gtex_kdp, gene_id), by = "gene_id", relationship = "one-to-one") |>
  replace_na(list(n = 0)) |>
  mutate(n_per_tissue = n / n_distinct(qtls_gtex_kdp$tissue)) |>
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
  left_join(count(qtls_gtex_ddp, gene_id), by = "gene_id", relationship = "one-to-one") |>
  replace_na(list(n = 0)) |>
  mutate(n_per_tissue = n / n_distinct(qtls_gtex_ddp$tissue)) |>
  summarise(mean_n_per_tissue = mean(n_per_tissue),
            .by = gene_biotype)

###

# "Expression had the lowest ratio of colocalizing TWAS hits to xQTLs (), and rDDPs had the second-lowest (, Figure 5d-f)."

twas_gtex_rddp <- read_tsv("data/processed/gtex-residual.twas_hits.tsv.gz", col_types = "ccc---dd") |>
  mutate(modality = "latent_residual", .before = 2)

coloc <- bind_rows(twas_gtex_kdp, twas_gtex_rddp) |>
  summarise(
    coloc_n = sum(coloc_pp > 0.8),
    .by = c(tissue, modality)
  )

coloc_qtl_ratio <- qtls_gtex_rddp |>
  count(tissue, modality, name = "n_qtls") |>
  full_join(coloc, by = c("tissue", "modality"), relationship = "one-to-one") |>
  mutate(ratio = coloc_n / n_qtls)

coloc_qtl_ratio |>
  summarise(mean_ratio = mean(ratio), .by = modality) |>
  arrange(mean_ratio)

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

qtls_geuv <- read_tsv("data/processed/geuvadis.qtls.tsv.gz", col_types = "ccccdci")

qtls_geuv |>
  count(version)

qtls_geuv |>
  filter(version == "residual-cross_latent") |>
  mutate(is_latent = modality == "latent_residual") |>
  count(is_latent)

# "In terms of unique genes represented by the xQTLs, addition of the latent residual phenotypes increased the number of xGenes from 12,872 to 16,950"

qtls_geuv |>
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

twas_geuv_rddp <- read_tsv("data/processed/geuvadis-residual.twas_hits.tsv.gz", col_types = "cc---d-") |>
  mutate(modality = "latent", .before = 1)

# "This resulted in a 58% increase in unique gene-trait association pairs"

bind_rows(twas_geuv_kdp, twas_geuv_rddp) |>
  summarise(latent_only = all(modality == "latent"),
            .by = c(trait, gene_id)) |>
  with(sum(latent_only) / sum(!latent_only))

# "and residual latent phenotypes had the strongest associations for 7,655 (53%) of all unique pairs"

bind_rows(twas_geuv_kdp, twas_geuv_rddp) |>
  slice_min(twas_p, n = 1, with_ties = FALSE, by = c(trait, gene_id)) |>
  summarise(n_latent = sum(modality == "latent"),
            frac_latent = mean(modality == "latent"))

###

# "Similarly, across GTEx tissues there was a 59% increase in unique gene-trait association pairs on average"

bind_rows(twas_gtex_kdp, twas_gtex_rddp) |>
  summarise(latent_only = all(modality == "latent_residual"),
            .by = c(tissue, trait, gene_id)) |>
  summarise(n_not_latent_only = sum(!latent_only),
            n_latent_only = sum(latent_only),
            .by = c(tissue)) |>
  mutate(percent_inc = n_latent_only / n_not_latent_only) |>
  with(mean(percent_inc))

# "and residual latent phenotypes had the strongest associations for 6,926 (51%) of all unique pairs on average"

bind_rows(twas_gtex_kdp, twas_gtex_rddp) |>
  slice_min(twas_p, n = 1, with_ties = FALSE, by = c(tissue, trait, gene_id)) |>
  summarise(n_latent = sum(modality == "latent_residual"),
            frac_latent = mean(modality == "latent_residual"),
            .by = tissue) |>
  summarise(mean_n_latent = mean(n_latent),
            mean_frac_latent = mean(frac_latent))

###

twas_examples <- tribble(
  ~tissue,   ~gene_id,          ~trait,
  "ADRNLG",  "ENSG00000134480", "UKB_1180_Morning_or_evening_person_chronotype",
  "HRTLV",   "ENSG00000178882", "GLGC_Mc_TG",
  "PNCREAS", "ENSG00000136267", "MAGIC_FastingGlucose",
  "BRNACC",  "ENSG00000239268", "UKB_1200_Sleeplessness_or_insomnia",
)

twas_example_hits <- twas_gtex_ddp |>
  semi_join(twas_examples, by = c("tissue", "gene_id", "trait")) |>
  left_join(select(genes, gene_id, gene_name), by = "gene_id", relationship = "many-to-one")

# Confirm they are present in rDDP hits but not in KDP hits

semi_join(twas_gtex_rddp, twas_examples, by = c("tissue", "gene_id", "trait"))
semi_join(twas_gtex_kdp, twas_examples, by = c("gene_id", "trait"))

# The “Morning/evening person (chronotype)” trait had data-driven xTWAS hits for CCNH (Cyclin H) in adrenal gland tissue.

filter(twas_example_hits, gene_name == "CCNH")

# The “Triglycerides” trait had data-driven xTWAS hits for RFLNA (Refilin A) in heart - left ventricle.

filter(twas_example_hits, gene_name == "RFLNA")

# The “Fasting glucose” trait had data-driven xTWAS hits for DGKB (Diacylglycerol kinase beta) in pancreas.

filter(twas_example_hits, gene_name == "DGKB")

# Finally, the “Sleeplessness / insomnia” trait had data-driven xTWAS hits for LINC03051 (Long intergenic non-protein coding RNA 3051) in .

filter(twas_example_hits, gene_name == "LINC03051")

###

qtls_gtex5_full <- read_tsv("data/processed/gtex5-full.qtls.tsv.gz", col_types = "cciccd")

qtls_gtex_full <- read_tsv("data/processed/gtex-full.qtls.tsv.gz", col_types = "cciccd")

qtl_counts_models <- full_join(
  qtls_gtex5_full |>
    count(tissue, name = "n_gtex5"),
  qtls_gtex_full |>
    count(tissue, name = "n_gtex"),
  by = "tissue",
  relationship = "one-to-one"
) |>
  full_join(
    qtls_gtex_ddp |>
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

qtl_counts_models |>
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

# "When we used each of these variations as input for latent RNA phenotyping and mapped cis-QTLs, we observed a relatively small drop (-2.7%) in discoveries from the 50 bp truncated reads, and larger drops from the single-end (-7.2%) and lowered sequencing depth (-41% for 50% of reads, -66% for 25% of reads) simulations"

(qtl_counts_seqsim["pe-50bp-100pct"] / qtl_counts_seqsim["pe-75bp-100pct"] - 1) * 100
(qtl_counts_seqsim["se-75bp-100pct"] / qtl_counts_seqsim["pe-75bp-100pct"] - 1) * 100
(qtl_counts_seqsim["pe-75bp-50pct"] / qtl_counts_seqsim["pe-75bp-100pct"] - 1) * 100
(qtl_counts_seqsim["pe-75bp-25pct"] / qtl_counts_seqsim["pe-75bp-100pct"] - 1) * 100

###

# "The original samples had an average of 30.1 million reads (SD 8.6 million)."

read_tsv("data/seqsim/read_counts.tsv", col_types = "ci") |>
  summarise(reads_mean = mean(reads),
            reads_sd = sd(reads))
