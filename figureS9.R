## Read length simulation QTLs

library(tidyverse)

read_types <- c(
  `pe-75bp-100pct` = "75 bp P.E.",
  `pe-50bp-100pct` = "50 bp P.E.",
  `se-75bp-100pct` = "75 bp S.E. R1&2",
  `se-50bp-100pct` = "50 bp S.E. R1&2",
  `se1-75bp-100pct` = "75 bp S.E. R1",
  `se1-50bp-100pct` = "50 bp S.E. R1",
  `se2-75bp-100pct` = "75 bp S.E. R2",
  `se2-50bp-100pct` = "50 bp S.E. R2",
  `pe-75bp-50pct` = "75 bp P.E. 1/2 depth",
  `pe-75bp-25pct` = "75 bp P.E. 1/4 depth"
)

qtls <- read_tsv("data/processed/seqsim.qtls.tsv.gz", col_types = "cccdci") |>
  mutate(reads = factor(read_types[reads], levels = read_types))

qtls |>
  count(reads) |>
  mutate(rel_amount = n / n[reads == "75 bp P.E."]) |>
  ggplot(aes(x = reads, y = rel_amount)) +
  geom_col(width = 0.5, fill = "black") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.04)), labels = scales::label_percent()) +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(hjust = 1, vjust = 1, angle = 45, color = "black"),
    axis.text.y = element_text(color = "black"),
  ) +
  xlab("Simulated read variations") +
  ylab("# xQTLs relative to # for 75 bp P.E.")

ggsave("figures/figureS9.png", width = 3.5, height = 3.6, device = png)
