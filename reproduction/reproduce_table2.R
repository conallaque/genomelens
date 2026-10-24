#!/usr/bin/env Rscript
# Recompute, from the printed totals of Kastrinos et al. 2021 (Gastroenterology, PMID 33839100), Table 2, the arithmetic its conclusions rest on,
# in base R only. The table compares colonoscopy schedules among carriers who are already known (a cohort of healthy 25-year-olds followed to 75,
# US health-care system perspective, 3% discounting as the paper states it, 2019 US dollars, perfect adherence). It has no "not identified" arm, so nothing here is the
# value of identifying a carrier. This reproduces the ARITHMETIC OF THE PRINTED TOTALS; it does not rebuild the paper's Markov model.
# 1. Every printed strategy is classified: on the efficiency frontier, or dominated.
# 2. Each printed ICER is set beside the one recomputed from the printed totals, allowing for the rounding of QALYs to three decimals.
# 3. The optimal strategy at 100,000 USD per QALY (the highest QALYs with an ICER no greater than the threshold) and its net monetary benefit,
#    NMB = threshold x incremental QALYs - incremental cost, against the next efficient strategy and against the least intensive strategy listed.
# 4. A check of the paper's statement that costs and QALYs are discounted at 3% a year: the most 3% discounting can give a cohort followed for 50
#    annual cycles (ages 25 to 75), at full health with no deaths, is compared with the printed QALYs.
# Usage: Rscript reproduction/reproduce_table2.R   (the output is compared with EXPECTED_OUTPUT.txt by tests/public).
args <- commandArgs(trailingOnly = FALSE)
here <- dirname(sub("^--file=", "", args[grep("^--file=", args)]))
if (length(here) == 0 || here == "") here <- "."
wtp <- 100000
if (!file.exists(file.path(here, "kastrinos2021_table2.tsv"))) {
  cat("kastrinos2021_table2.tsv is not shipped with this repository: the publisher's terms for the NIH author manuscript do not grant redistribution.\nSee reproduction/README.md for how to supply it (five columns, 20 printed strategies). Nothing was computed.\n", file = stderr())
  quit(status = 2)
}
t2 <- read.delim(file.path(here, "kastrinos2021_table2.tsv"), comment.char = "#", stringsAsFactors = FALSE, check.names = FALSE, quote = "")
money <- function(s) as.numeric(gsub("[$,]", "", s))
f2 <- function(x, nd = 2) formatC(x, format = "f", digits = nd, big.mark = ",")
sgn <- function(x, nd = 2) paste0(ifelse(x < 0, "-", "+"), f2(abs(x), nd))
classify <- function(q, cst) {
  n <- length(q)
  strong <- sapply(seq_len(n), function(i) any(sapply(setdiff(seq_len(n), i), function(j)
    cst[j] <= cst[i] && q[j] >= q[i] && (cst[j] < cst[i] || q[j] > q[i]))))
  alive <- order(cst); alive <- alive[!strong[alive]]
  ext <- integer(0)
  repeat {
    drop <- 0
    if (length(alive) >= 3) for (k in 2:(length(alive) - 1)) {
      a <- alive[k - 1]; b <- alive[k]; c <- alive[k + 1]
      if ((cst[b] - cst[a]) / (q[b] - q[a]) > (cst[c] - cst[b]) / (q[c] - q[b])) { drop <- k; break }
    }
    if (drop == 0) break
    ext <- c(ext, alive[drop]); alive <- alive[-drop]
  }
  list(status = ifelse(strong, "dominated", ifelse(seq_len(n) %in% ext, "extended-dominated", "frontier")), front = alive)
}
genes <- unique(t2$gene)
res <- lapply(genes, function(g) {
  d <- t2[t2$gene == g, ]
  q <- as.numeric(d$QALYs); cst <- money(d$Cost)
  c(list(d = d, q = q, cst = cst), classify(q, cst))
})
names(res) <- genes
cat("1. Classification of every printed strategy\n")
for (g in genes) {
  r <- res[[g]]
  for (i in seq_len(nrow(r$d))) cat(sprintf("  %s %s: %s (printed %s)\n", g, r$d$strategy_as_printed[i], r$status[i], r$d$ICER[i]))
}
cat("2. Printed ICERs against the previous strategy on the frontier, recomputed from the printed totals\n")
icers <- list()
for (g in genes) {
  r <- res[[g]]; front <- r$front
  icers[[g]] <- numeric(0)
  if (length(front) > 1) for (k in 2:length(front)) {
    a <- front[k - 1]; b <- front[k]
    dq <- r$q[b] - r$q[a]; dc <- r$cst[b] - r$cst[a]
    lo <- dc / (dq + 0.001); hi <- dc / (dq - 0.001)
    printed <- money(r$d$ICER[b])
    ok <- (lo - 1 <= printed) && (printed <= hi + 1)
    icers[[g]] <- c(icers[[g]], dc / dq)
    line <- sprintf("  %s %s vs %s: printed %s, recomputed %s, ", g, r$d$strategy_as_printed[b], r$d$strategy_as_printed[a], f2(printed), f2(dc / dq))
    line <- paste0(line, if (ok) "reproduces" else sprintf("does not reproduce (QALYs rounded to three decimals allow %s to %s; the printed ratio implies an increment of %.4f QALYs, the printed totals give %.4f)",
                                                       f2(lo), f2(hi), dc / printed, dq))
    cat(line, "\n", sep = "")
  }
}
cat("3. The optimal strategy at 100,000 USD per QALY, and its net monetary benefit against two comparators\n")
for (g in genes) {
  r <- res[[g]]; front <- r$front
  elig <- c(front[1], if (length(front) > 1) front[2:length(front)][icers[[g]] <= wtp])
  opt <- elig[which.max(r$q[elig])]
  cat(sprintf("  %s optimal: %s\n", g, r$d$strategy_as_printed[opt]))
  k <- which(front == opt)
  if (k > 1) {
    p <- front[k - 1]; dq <- r$q[opt] - r$q[p]; dc <- r$cst[opt] - r$cst[p]
    cat(sprintf("    against the next efficient strategy (%s): QALYs %s, cost %s, NMB %s\n", r$d$strategy_as_printed[p], sgn(dq, 3), sgn(dc), f2(wtp * dq - dc)))
  } else cat("    against the next efficient strategy: none (the optimal strategy is the reference of the frontier)\n")
  dq <- r$q[opt] - r$q[1]; dc <- r$cst[opt] - r$cst[1]
  cat(sprintf("    against the least intensive strategy listed (%s): QALYs %s, cost %s, NMB %s\n", r$d$strategy_as_printed[1], sgn(dq, 3), sgn(dc), f2(wtp * dq - dc)))
}
cat("4. A check of the stated discounting (the source states 3% a year, ages 25 to 75)\n")
cap <- sum(1.03 ^ -(0:49)); cap_end <- sum(1.03 ^ -(1:50))
cat(sprintf("  The most 3%% discounting can give over 50 annual cycles, at full health with no deaths, is %.2f QALYs (%.2f if the first cycle is also discounted)\n", cap, cap_end))
allover <- TRUE
for (g in genes) {
  qs <- res[[g]]$q
  allover <- allover && (min(qs) > cap)
  cat(sprintf("  %s: printed QALYs %.3f to %.3f, %s %.2f\n", g, min(qs), max(qs), if (min(qs) > cap) "all above" else "not all above", cap))
}
if (allover) cat("  The printed QALYs therefore cannot be 3% discounted over that horizon: the discounting convention behind the printed totals is unresolved, and every figure above is arithmetic on the totals as printed.\n")
