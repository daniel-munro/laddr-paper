## Create simple plots for use in method diagram

library(tidyverse)

#########################
## base-level coverage ##
#########################

bw <- rtracklayer::import.bw("~/GitHub/LaDDR/examples/data_input/covg_bigwig/dset1/HG00108.bw")

covg <- as_tibble(bw) |>
    filter(start > 16790,
           end < 17900) |>
    mutate(score = if_else(start > 17290 & start < 17330, score / 10, score),
           # group = round(start, -1)) |>
           group = floor(start / 4)) |>
    summarise(start = min(start),
              end = max(end),
              score = mean(score),
              .by = group)
covg |>
    # filter(start > 17100,
    #        end < 18500) |>
    ggplot(aes(xmin = start, xmax = end + 1, ymin = 0, ymax = score)) +
    geom_rect(fill = "black") +
    theme_classic() +
    theme(
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        plot.margin = margin(0, 0, 0, 0),
    )
ggsave("figures/figure1/figure1a_part1.png", width = 3, height = 0.7)

##########
## Bins ##
##########

bins <- read_tsv(
    "data/adaptive/bins_variance_diff/384.bed.gz",
    col_types = "ciiccc",
    col_names = c("chrom", "start", "end", "name", "type", "strand")
) |>
    separate(name, c("gene_id", "start_gene", "end_gene"), sep = "_", convert = TRUE) |>
    filter(gene_id == "ENSG00000131591") |> #unique(gene_id)[6]) |>
    filter(start > 1090900,
           end < 1111500) |>
    filter(start > lag(start + 24)) |>
    select(start, end)
bins |>
    bind_rows(tibble(start = min(bins$start),
                     end = max(bins$end))) |>
    ggplot(aes(xmin = start, xmax = end, ymin = 0.6, ymax = 1.4)) +
    geom_rect(color = "black", fill = NA, linewidth = 0.2) +
    theme_classic() +
    theme(
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        plot.margin = margin(0, 0, 0, 0),
    )
ggsave("figures/figure1/figure1a_part2.png", width = 4, height = 0.2)

#####################
## Binned coverage ##
#####################

covg |>
    mutate(group = round(cumsum(score) / 80)) |>
    summarise(score = mean(score),
              .by = group) |>
    mutate(group = rank(group)) |>
    # ggplot(aes(xmin = group, xmax = group + 0.6, ymin = 0, ymax = score)) +
    # geom_rect(fill = "black") +
    ggplot(aes(x = group, y = score)) +
    geom_col(fill = "black", width = 0.6) +
    theme_classic() +
    theme(
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        plot.margin = margin(0, 0, 0, 0),
    ) +
    xlab(NULL) +
    ylab(NULL)
ggsave("figures/figure1/figure1a_part3.png", width = 3, height = 0.7)
