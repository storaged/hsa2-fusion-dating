# A0: Reproduce the published HSA2 fusion-time estimate (Poszewiecka et al. 2022)
# from the legacy UBCS tables (ubcs_2021/wyniki_2021_03_07_v4.tsv), following
# legacy/ubcs_2021/fuzja_KG.R, then run sensitivity analyses on choices that the
# original code fixed implicitly.
#
# Legacy logic (fuzja_KG.R):
#   ratio_2b  = sum(UBCS human-derived) / sum(UBCS chimp-derived) over chr2:[114, 129) Mb
#   ratio_ctl = mean over control telomeres of the same ratio over the terminal D Mb
#   T         = 6 * (1 - ratio_2b / ratio_ctl)
#   CI        = bootstrap over control telomeres (7 draws) and D in 10..15 Mb (7 draws),
#               median over D, 200 replicates. The chr2 ratio is NOT resampled.

suppressPackageStartupMessages({ library(dplyr); library(tidyr) })
set.seed(20260927)

args <- commandArgs(trailingOnly = TRUE)
legacy_file <- if (length(args) >= 1) args[1] else "legacy/ubcs_2021/wyniki_2021_03_07_v4.tsv"
out_dir     <- if (length(args) >= 2) args[2] else "results/A0"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

FUSION_SITE <- 114e6
T_SPLIT     <- 6

read_ubcs <- function(f) {
  d <- read.table(f, header = FALSE, stringsAsFactors = FALSE)
  colnames(d) <- c("chrom", "begin", "end", "window_size", "number_of_bins", "derived_in_type",
                   "target", "query", "outgroup", "ffrom", "with_snps", "p",
                   "expected_BCS", "actual_BCS", "UBCS")
  d$with_snps <- d$with_snps == "True"
  d <- distinct(d)
  last <- d %>% group_by(chrom) %>% summarise(last_begin = max(begin), .groups = "drop")
  d <- inner_join(d, last, by = "chrom")
  d[d$begin != d$last_begin, ]   # drop the partial last Mb, as in utils.R
}

select_set <- function(d, bins, ffrom = NA, query = "panTro6") {
  s <- d[d$target == "hg38" & d$query == query & d$with_snps & d$window_size == 300 &
         d$number_of_bins == bins, ]
  if (!is.na(ffrom)) s <- s[s$ffrom == ffrom, ]
  s
}

# control telomeres: legacy default {1..12, 16, 17}, both arms
LEGACY_CTRL <- paste0("chr", c(1:12, 16, 17))

telomere_ratios <- function(dt, D, ctrl = LEGACY_CTRL) {
  st <- dt[dt$chrom %in% ctrl, ] %>%
    filter(begin < D | begin + D >= last_begin) %>%
    mutate(tel = paste0(chrom, "_", ifelse(begin < D, "p", "q"))) %>%
    group_by(tel) %>%
    summarise(r = sum(UBCS[derived_in_type == "target"]) / sum(UBCS[derived_in_type == "query"]),
              .groups = "drop")
  setNames(st$r, st$tel)
}

fusion_ratio <- function(dt, side = c("2b", "2a"), D = 15e6) {
  side <- match.arg(side)
  c2 <- dt[dt$chrom == "chr2", ]
  w <- if (side == "2b") c2$begin >= FUSION_SITE & c2$begin < FUSION_SITE + D
       else              c2$begin <  FUSION_SITE & c2$begin >= FUSION_SITE - D
  x <- c2[w, ]
  sum(x$UBCS[x$derived_in_type == "target"]) / sum(x$UBCS[x$derived_in_type == "query"])
}

# --- 1. exact legacy procedure -----------------------------------------------
legacy_estimate <- function(dt, nboot = 200, agg = mean, Ds = 10:15, D_fusion = 15e6,
                            side = "2b", resample_fusion = FALSE, ctrl = LEGACY_CTRL) {
  r2 <- fusion_ratio(dt, side, D_fusion)
  point <- median(sapply(Ds, function(D) T_SPLIT * (1 - r2 / agg(telomere_ratios(dt, D * 1e6, ctrl)))))
  boot <- replicate(nboot, {
    Dsel <- sample(Ds, length(Ds) + 1, replace = TRUE)
    median(sapply(Dsel, function(D) {
      rt <- telomere_ratios(dt, D * 1e6, ctrl)
      rt <- sample(rt, 7, replace = TRUE)
      r2b <- r2
      if (resample_fusion) {   # resample 1-Mb bins of the fusion flank too
        c2 <- dt[dt$chrom == "chr2", ]
        w  <- if (side == "2b") c2$begin >= FUSION_SITE & c2$begin < FUSION_SITE + D_fusion
              else              c2$begin <  FUSION_SITE & c2$begin >= FUSION_SITE - D_fusion
        x <- c2[w, ]; b <- sample(unique(x$begin), replace = TRUE)
        tt <- sapply(b, function(bb) x$UBCS[x$begin == bb & x$derived_in_type == "target"])
        qq <- sapply(b, function(bb) x$UBCS[x$begin == bb & x$derived_in_type == "query"])
        r2b <- sum(unlist(tt)) / sum(unlist(qq))
      }
      T_SPLIT * (1 - r2b / agg(rt))
    }))
  })
  c(point = point, lo = unname(quantile(boot, 0.025)), med = median(boot),
    hi = unname(quantile(boot, 0.975)), ratio_fusion = r2)
}

d <- read_ubcs(legacy_file)
cat("Rows:", nrow(d), " queries:", paste(unique(d$query), collapse = ","), "\n")

sets <- list(bins300 = select_set(d, 300), bins20 = select_set(d, 20), bins2_q = select_set(d, 2, "query"))
res <- list()
for (nm in names(sets)) {
  dt <- sets[[nm]]
  if (nrow(dt) == 0) { cat("empty set", nm, "\n"); next }
  res[[paste0(nm, "|legacy")]]             <- legacy_estimate(dt)
  res[[paste0(nm, "|median_ctrl")]]        <- legacy_estimate(dt, agg = median)
  res[[paste0(nm, "|side2a")]]             <- legacy_estimate(dt, side = "2a")
  res[[paste0(nm, "|fusionD10")]]          <- legacy_estimate(dt, D_fusion = 10e6)
  res[[paste0(nm, "|fusionD20")]]          <- legacy_estimate(dt, D_fusion = 20e6)
  res[[paste0(nm, "|resample_fusion")]]    <- legacy_estimate(dt, resample_fusion = TRUE)
  res[[paste0(nm, "|all_noninv_ctrl")]]    <- legacy_estimate(dt, ctrl = paste0("chr", c(1, 3:12, 16, 17, 19, 20)))
}
out <- do.call(rbind, res)
out <- data.frame(set = sub("\\|.*", "", rownames(out)), variant = sub(".*\\|", "", rownames(out)),
                  round(out, 3), row.names = NULL)
print(out)
write.table(out, file.path(out_dir, "A0_reproduction_sensitivity.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

# per-telomere control ratios (for the plan/figures)
tr <- telomere_ratios(sets$bins300, 15e6, paste0("chr", 1:22))
write.table(data.frame(telomere = names(tr), ratio_human_over_chimp = tr),
            file.path(out_dir, "A0_control_telomere_ratios_D15.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
cat("fusion 2b ratio:", fusion_ratio(sets$bins300, "2b"), " 2a ratio:", fusion_ratio(sets$bins300, "2a"), "\n")
