## Pruned annotation xQTLs

modalities <- c(
  latent = "Data-driven",
  alt_TSS = "Alt. TSS",
  alt_polyA = "Alt. polyA",
  isoforms = "Isoform ratio",
  stability = "RNA stability",
  splicing = "Intron excision",
  expression = "Expression"
)

modality_colors <- c(
  `Data-driven` = "#13918d",
  `Alt. TSS` = "#896090",
  `Alt. polyA` = "#d97f26",
  `Isoform ratio` = "#6a90cd",
  `RNA stability` = "#ddb23c",
  `Intron excision` = "#59a257",
  Expression = "#bf4042"
)

map_groups <- c(
  latent = "Data-driven",
  pantry = "Knowledge-driven"
)

qtls_prune <- read_tsv("data/processed/prune-BRNCTXB.qtls.tsv.gz", col_types = "cicciccd") |>
  filter(map_group %in% names(map_groups)) |>
  mutate(modality = factor(modalities[modality], levels = names(modality_colors)),
         map_group = factor(map_groups[map_group], levels = map_groups),
         pruning = fct_reorder(as.character(pruning), pruning))

n_qtls_100p <- qtls_prune |>
  filter(pruning == 100) |>
  count(map_group)

qtls_prune |>
  count(map_group, pruning, modality) |>
  filter(map_group %in% c("Data-driven", "Knowledge-driven")) |>
  ggplot(aes(x = pruning, y = n / 1000, fill = modality)) +
  # facet_wrap(~map_group, scales = "free_y") +
  facet_wrap(~map_group) +
  geom_hline(mapping = aes(yintercept = n / 1000), data = n_qtls_100p,
             linetype = 2) +
  geom_col(width = 0.8) +
  # scale_y_continuous(breaks = c(0, 5, 10, 15, 20), expand = expansion(mult = c(0, 0.03))) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.03))) +
  scale_fill_manual(values = modality_colors) +
  theme_bw() +
  theme(
    axis.text = element_text(color = "black"),
    legend.key.size = unit(12, "pt"),
    panel.grid = element_blank(),
  ) +
  xlab("% of non-canonical isoforms kept") +
  ylab("xQTLs (×1000)") +
  labs(fill = "Modality")

ggsave("figures/figure5/figure5c.png", width = 5, height = 3, device = png)
