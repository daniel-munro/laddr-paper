## Replication of GTEx LCL DDP xQTLs in Geuvadis

library(tidyverse)

qtlrep <- read_tsv(
  "data/qtl/gtextcga-full_LCL_ddp_in_geuvadis.tsv.gz",
  col_types = "ccdddidddd"
) |>
  mutate(
    geuvadis_signif = geuvadis_pval_nominal < geuvadis_pval_nominal_threshold,
    gtex_slope = pmax(-4, pmin(gtex_slope, 4))
  ) |>
  # arrange(geuvadis_signif)
  arrange(desc(gtex_pval_nominal)) |>
  mutate(gtex_pval_nominal = pmax(1e-40, gtex_pval_nominal))

qtlrep |>
  ggplot(aes(x = gtex_slope, y = geuvadis_slope, color = -log10(gtex_pval_nominal), alpha = geuvadis_signif)) +
  geom_hline(yintercept = 0, linewidth = 0.3, color = "#cccccc") +
  geom_vline(xintercept = 0, linewidth = 0.3, color = "#cccccc") +
  geom_point(size = 0.5) +
  scale_x_continuous(expand = 0) +
  scale_y_continuous(limits = c(-4, 4), expand = 0) +
  scale_color_viridis_c(breaks = c(10, 20, 30, 40), labels = c(10, 20, 30, "40+")) +
  scale_alpha_manual(values = c(0.2, 1)) +
  coord_fixed() +
  theme_classic() +
  guides(alpha = "none") +
  xlab("xQTL slope, GTEx LCL") +
  ylab("Slope for same phenotype-variant, Geuvadis LCL") +
  labs(color = expression(-log[10](italic(p)[plain("GTEx")])))

ggsave("figures/figureS2.png", width = 5, height = 4, device = png)
