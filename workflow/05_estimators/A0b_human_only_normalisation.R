# A0b: human-only normalisation of the legacy UBCS signal.
# R = (UBCS_H at fusion flank - interior baseline) / (mean UBCS_H at control ends - interior baseline),
# i.e. fraction of "telomere-like" signal retained at the fusion flank on the human branch,
# without dividing by chimp. Same for the chimp branch as a control (should be ~0 at 2a flank in Pan,
# which is not subtelomeric there, and ~1 at 2b flank, which is).
suppressPackageStartupMessages(library(dplyr))
src <- readLines("scripts/A0_reproduce_2022.R")
eval(parse(text = src[1:grep("^LEGACY_CTRL", src)]))
d <- read_ubcs("legacy/ubcs_2021/wyniki_2021_03_07_v4.tsv")
FUSION <- 113.57e6
ctrl_chr <- paste0("chr", c(1, 3:12, 16, 17))
rows <- list()
for (b in c(300, 20, 2)) {
  dt <- select_set(d, b, if (b == 2) "query" else NA)
  for (D in c(5, 10, 15) * 1e6) for (who in c("target", "query")) {
    x <- dt[dt$derived_in_type == who, ]
    ends <- x[x$chrom %in% ctrl_chr, ] %>% filter(begin < D | begin + D >= last_begin) %>%
      mutate(tel = paste0(chrom, ifelse(begin < D, "p", "q"))) %>% group_by(tel) %>%
      summarise(s = sum(UBCS), .groups = "drop")
    mid <- x[x$chrom %in% ctrl_chr & x$begin > 40e6 & x$begin + 40e6 < x$last_begin, ] %>%
      group_by(chrom) %>% summarise(s = mean(UBCS) * D / 1e6, .groups = "drop")
    c2 <- x[x$chrom == "chr2", ]
    f2b <- sum(c2$UBCS[c2$begin >= floor(FUSION / 1e6) * 1e6 + 1e6 & c2$begin < FUSION + D])
    f2a <- sum(c2$UBCS[c2$begin < floor(FUSION / 1e6) * 1e6 & c2$begin >= FUSION - D])
    base <- mean(mid$s); tel <- mean(ends$s)
    rows[[length(rows) + 1]] <- data.frame(bins = b, D_Mb = D / 1e6,
      branch = ifelse(who == "target", "human", "chimp"), ends_mean = tel, ends_sd = sd(ends$s),
      interior = base, flank2a = f2a, flank2b = f2b,
      R2a = (f2a - base) / (tel - base), R2b = (f2b - base) / (tel - base),
      q2b_among_ends = mean(ends$s <= f2b))
  }
}
out <- bind_rows(rows)
dir.create("results/A0", showWarnings = FALSE, recursive = TRUE)
write.table(out, "results/A0/A0b_human_only_normalisation.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
print(out %>% mutate(across(where(is.numeric), ~ round(.x, 2))), row.names = FALSE)
