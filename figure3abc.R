## Latent xQTLs and TWAS

library(tidyverse)

#############
## Panel a ## Full latent vs. explicit xQTLs bar plot
#############

qtls_gtextcga_full <- read_tsv("data/processed/gtextcga-full.qtls.tsv.gz", col_types = "cciccd")

qtls_pantry <- read_tsv("data/processed/gtex-pantry.qtls.tsv.gz", col_types = "ccicccd")

qtls <- bind_rows(
  qtls_gtextcga_full |>
    mutate(type = "Data-driven"),
  qtls_pantry |>
    mutate(type = "Knowledge-driven")
) |>
  count(type, tissue) |>
  arrange(desc(type), desc(n)) |>
  mutate(tissue = fct_inorder(tissue),
         type = fct_inorder(type))

ggplot(qtls, aes(x = tissue, y = n / 1000, fill = type)) +
  geom_col(position = "dodge") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.04))) +
  scale_fill_manual(values = c("#ad611f", "#13918d")) +
  # expand_limits(y = max(qtls$n * 1.03 / 1000)) +
  theme_bw() +
  theme(
    axis.text.x = element_blank(),
    axis.text.y = element_text(color = "black"),
    axis.ticks.x = element_blank(),
    legend.key.size = unit(12, "pt"),
    legend.position = "inside",
    legend.position.inside = c(0.75, 0.8),
    panel.grid = element_blank(),
  ) +
  xlab("GTEx tissues") +
  ylab("xQTLs (×1000)") +
  labs(fill = NULL)

ggsave("figures/figure3/figure3a.png", width = 3.5, height = 3, device = png)

#############
## Panel b ## Full latent vs. explicit xQTLs scatter plot
#############

gtex_colors <- read_tsv(
  "data/pantry/gtex/tissueInfo.tsv",
  col_types = cols(tissueSiteDetailAbbr = "c", colorHex = "c", .default = "-")
) |>
  mutate(colorHex = str_c("#", colorHex)) |>
  deframe()

qtls_wide <- qtls |>
  pivot_wider(id_cols = tissue, names_from = "type", values_from = "n")

ggplot(qtls_wide, aes(x = `Knowledge-driven` / 1000, y = `Data-driven` / 1000, color = tissue)) +
  geom_abline(slope = 1, intercept = 0, linetype = 3) +
  geom_point(show.legend = FALSE) +
  expand_limits(
    x = c(0, max(qtls_wide$`Knowledge-driven`) * 1.04 / 1000),
    y = c(0, max(qtls_wide$`Data-driven`) * 1.04 / 1000),
  ) +
  scale_y_continuous(breaks = c(0, 10, 20, 30, 40, 50, 60)) +
  scale_color_manual(values = gtex_colors) +
  coord_fixed(expand = 0) +
  theme_bw() +
  theme(
    axis.text = element_text(color = "black"),
    panel.grid = element_blank(),
  ) +
  xlab("KDP xQTLs (×1000)") +
  ylab("DDP xQTLs (×1000)") +
  labs(fill = NULL)

ggsave("figures/figure3/figure3b.png", width = 2.5, height = 3.5, device = png)

#############
## Panel c ## xQTLs by PC number
#############

qtls_pc_count <- qtls_gtextcga_full |>
  mutate(PC = str_split_i(phenotype_id, "__PC", 2) |>
           as.integer()) |>
  count(PC) |>
  mutate(mean_per_tissue = n / n_distinct(qtls_gtextcga_full$tissue))

qtls_pc_count |>
  ggplot(aes(x = PC, y = mean_per_tissue / 1000)) +
  geom_col(width = 0.5, fill = "black") +
  scale_x_continuous(expand = c(0.03, 0), breaks = 1:16) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.01)), breaks = c(0, 2, 4, 6, 8)) +
  # expand_limits(y = 1.01 * max(qtls_pc_count$mean_per_tissue) / 1000) +
  theme_classic() +
  theme(
    axis.text = element_text(color = "black"),
  ) +
  xlab("DDP rank in gene") +
  ylab("Mean xQTLs per tissue (×1000)")

ggsave("figures/figure3/figure3c.png", width = 2.75, height = 3.5, device = png)
