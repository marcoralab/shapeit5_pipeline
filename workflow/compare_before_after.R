library(readr)
suppressPackageStartupMessages(library(dplyr))
library(ggplot2)

before <-
  "ADSP.WGS.SNPs.compact_filter.combined.chr22.withac.bcf.stats.af" |>
  read_tsv(col_types = cols(.default = "i", af = "d"))

after <-
  "ADSP.WGS.SNPs.compact_filter.combined.chr22.phased_rare.bcf.stats.af" |>
  read_tsv(col_types = cols(.default = "i", af = "d"))

binned_sums <- function(x) {
  x |>
    mutate(bin = cut(af, breaks = c(-1, 0, 1e-5, 1e-4, 1e-3, (1:100) / 100)),
           totvars = snp + indel) |>
    group_by(bin) |>
    summarise(af = weighted.mean(af, totvars),
              across(where(is.integer), sum))
}

before |>
  anti_join(after, by = "af")

after |>
  anti_join(before, by = "af")


after_binned <- after |>
  binned_sums() |>
  mutate(stage = "after")

before_binned <- before |>
  binned_sums() |>
  mutate(stage = "before")

binned <- bind_rows(after_binned, before_binned)

plt <- binned |>
  ggplot(aes(x = af, color = stage)) +
    geom_line(aes(y = snp), linetype = "dashed") +
    geom_line(aes(y = indel), linetype = "dotted") +
    geom_line(aes(y = totvars)) +
    scale_x_log10() +
    theme_bw() +
    theme(legend.position = "bottom")

ggsave("compare_before_after.png", plot = plt, width = 8, height = 6, dpi = 300)

plt2 <- binned |>
  group_by(stage) |>
  mutate(bin_number = seq_along(bin)) |>
  ungroup() |>
  filter(af < 3e-2) |>
  ggplot(aes(x = bin_number, color = stage)) +
    geom_line(aes(y = snp), linetype = "dashed") +
    geom_line(aes(y = indel), linetype = "dotted") +
    geom_line(aes(y = totvars)) +
    theme_bw() +
    theme(legend.position = "bottom")

ggsave("compare_before_after_zoom.png", plot = plt2, width = 8, height = 6, dpi = 300)
