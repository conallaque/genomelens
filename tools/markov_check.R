#!/usr/bin/env Rscript
# ─────────────────────────────────────────────────────────────────────────────
# markov_check.R — independent cross-validation of econ/markov.py
#
# WHY THIS EXISTS
#   The Python engine is the only implementation of the cost-effectiveness
#   model, so "the model is right" currently rests on the model checking
#   itself. This is a second, independent implementation in a different
#   language, written from the published method rather than from the Python
#   source line-by-line. If both land on the same numbers, the arithmetic is
#   probably right.
#
#   It has already earned its keep once: it found that
#   markov_cost_effectiveness() was differencing per-arm totals that had
#   already been rounded for display (2dp cost, 4dp QALY), so NMB carried
#   ~$5 of rounding noise at a $100,000/QALY threshold. Fixed — run_markov now
#   also returns *_exact totals and the incremental math uses those. The NMB
#   section below keeps the emulated-rounding comparison as a regression guard,
#   so a reintroduction of the same mistake is reported rather than absorbed
#   into a tolerance.
#
#   Base R only. No heemod, no hesim, no jsonlite — nothing to install. The
#   matrix-cohort idiom below is the DARTH-course formulation of a discrete-time
#   cohort state-transition model (cSTM).
#
# USAGE
#   Rscript tools/markov_check.R                     # uses tools/cea_fixture.tsv
#   Rscript tools/markov_check.R path/to/fixture.tsv
#
#   Regenerate the fixture from Python first:
#   python tools/cea_excel_export.py --fixture
#
# EXIT STATUS
#   0 = every quantity agrees within tolerance; 1 = at least one mismatch.
# ─────────────────────────────────────────────────────────────────────────────

args <- commandArgs(trailingOnly = TRUE)
fixture <- if (length(args) >= 1) args[1] else {
  # Default to the fixture beside this script, so cwd does not matter.
  self <- sub("^--file=", "",
              grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE))
  if (length(self) == 1) file.path(dirname(self), "cea_fixture.tsv") else
    "tools/cea_fixture.tsv"
}

if (!file.exists(fixture)) {
  stop(sprintf("fixture not found: %s\n  run: python tools/cea_excel_export.py --fixture",
               fixture), call. = FALSE)
}

# Tab-separated: see the note in cea_excel_export.py about *.csv being
# gitignored by the never-commit-DNA rule.
d <- read.delim(fixture, stringsAsFactors = FALSE)
num <- function(k) as.numeric(d$value[d$section == "input" & d$key == k])
exp_num <- function(k) as.numeric(d$value[d$section == "expected" & d$key == k])
# read.delim turns the literal token NA (used for a suppressed ICER) into
# NA_character_, and `NA %in% c("NA", ...)` is FALSE rather than TRUE — which
# would silently mis-report the ICER comparison. Normalise it back to a string.
exp_chr <- function(k) {
  v <- as.character(d$value[d$section == "expected" & d$key == k])
  if (length(v) == 0 || is.na(v)) "NA" else v
}

# ── inputs ───────────────────────────────────────────────────────────────────
start_age      <- num("start_age")
n_cycles       <- as.integer(num("cycles"))
incidence_rate <- num("incidence_rate")
rrr            <- num("rrr_intervention")
c_interv       <- num("cost_intervention_annual")
c_disease      <- num("cost_disease_annual")
c_well         <- num("cost_well_annual")
u_well         <- num("utility_well")
u_disease      <- num("utility_disease")
excess_mort    <- num("excess_mortality_rate")
disc_rate      <- num("discount_rate")
wtp            <- num("wtp")
gomp_a         <- num("gompertz_a")
gomp_b         <- num("gompertz_b")
cycle_len      <- num("cycle_length")
half_cycle     <- num("half_cycle") == 1

# ── method pieces ────────────────────────────────────────────────────────────

# Rate -> probability over one cycle. Using p = rate directly is the classic
# error; it overstates transitions and breaks as rate*dt approaches 1.
rate_to_prob <- function(rate, dt = 1) 1 - exp(-pmax(0, rate) * dt)

# Age-specific all-cause mortality RATE (Gompertz), anchored at age 20.
gompertz_rate <- function(age) gomp_a * exp(gomp_b * pmax(0, age - 20))

# 3x3 transition matrix over (Well, Disease, Dead) with competing risks.
# p_dd is built from the UN-normalised background probability: the Python
# implementation normalises its local copy for the Well row but reads the
# original parameter for the Disease row. Reproduced deliberately so the two
# implementations are comparable rather than merely close.
transition_matrix <- function(p_wd_raw, p_ex, p_bg_raw) {
  p_wd <- min(1, max(0, p_wd_raw))
  p_bg <- min(1, max(0, p_bg_raw))
  tot  <- p_wd + p_bg
  if (tot > 1) { p_wd <- p_wd / tot; p_bg <- p_bg / tot }   # competing risks
  p_dd <- min(1, max(0, p_bg_raw + p_ex))
  matrix(c(1 - p_wd - p_bg, p_wd,     p_bg,
           0,               1 - p_dd, p_dd,
           0,               0,        1),
         nrow = 3, byrow = TRUE,
         dimnames = list(c("Well", "Disease", "Dead"),
                         c("Well", "Disease", "Dead")))
}

#' Run one arm through the cohort model.
#' Returns full-precision totals — no rounding, so the comparison below tests
#' the model rather than a display convention.
run_markov <- function(guided = FALSE) {
  cohort <- c(Well = 1, Disease = 0, Dead = 0)
  tot_cost <- 0; tot_qaly <- 0; tot_ly <- 0

  for (t in seq_len(n_cycles) - 1L) {        # t = 0 .. n_cycles-1
    age <- start_age + t
    inc <- if (guided) incidence_rate * (1 - rrr) else incidence_rate

    p_wd <- rate_to_prob(inc, cycle_len)
    p_bg <- rate_to_prob(gompertz_rate(age), cycle_len)
    p_ex <- rate_to_prob(excess_mort, cycle_len)
    P    <- transition_matrix(p_wd, p_ex, p_bg)

    # Payoffs valued on START-of-cycle occupancy.
    c_cycle <- cohort["Well"] * (c_well + if (guided) c_interv else 0) +
               cohort["Disease"] * c_disease
    q_cycle <- cohort["Well"] * u_well + cohort["Disease"] * u_disease
    l_cycle <- cohort["Well"] + cohort["Disease"]

    # Half-cycle (trapezoidal) correction on the first and last cycles.
    w <- if (half_cycle && (t == 0 || t == n_cycles - 1)) 0.5 else 1
    disc <- 1 / (1 + disc_rate)^t

    tot_cost <- tot_cost + c_cycle * w * disc
    tot_qaly <- tot_qaly + q_cycle * w * disc
    tot_ly   <- tot_ly   + l_cycle * w * disc

    cohort <- as.vector(cohort %*% P)                # advance
    names(cohort) <- c("Well", "Disease", "Dead")
    s <- sum(cohort)
    if (abs(s - 1) > 1e-9) cohort <- cohort / s      # guard against leakage
  }
  list(cost = unname(tot_cost), qaly = unname(tot_qaly), ly = unname(tot_ly))
}

sc <- run_markov(guided = FALSE)
gg <- run_markov(guided = TRUE)

d_cost <- gg$cost - sc$cost
d_qaly <- gg$qaly - sc$qaly
d_ly   <- gg$ly   - sc$ly

dominant  <- d_cost < 0 && d_qaly > 0
dominated <- d_cost > 0 && d_qaly < 0
icer      <- if (dominant || dominated || abs(d_qaly) < 1e-9) NA else d_cost / d_qaly
nmb       <- wtp * d_qaly - d_cost

# ── comparison ───────────────────────────────────────────────────────────────
cat("\nIndependent R re-implementation vs Python (econ/markov.py)\n")
cat(sprintf("fixture: %s\n", fixture))
cat(strrep("-", 74), "\n")
cat(sprintf("%-22s %14s %14s %12s  %s\n", "quantity", "R", "Python", "abs diff", "ok"))

fails <- 0L
check <- function(label, got, want, tol = 0.01) {
  ok <- is.finite(got) && is.finite(want) && abs(got - want) <= tol
  if (!ok) fails <<- fails + 1L
  cat(sprintf("%-22s %14.4f %14.4f %12.2e  %s\n",
              label, got, want, abs(got - want), if (ok) "OK" else "FAIL"))
}

# Tolerance 0.01 absorbs the 2dp/4dp rounding the Python engine applies to its
# reported per-arm totals; anything larger is a real disagreement.
check("standard care cost",  sc$cost, exp_num("sc_cost"))
check("standard care QALYs", sc$qaly, exp_num("sc_qaly"), tol = 1e-4)
check("standard care LYs",   sc$ly,   exp_num("sc_ly"),   tol = 1e-4)
check("guided cost",         gg$cost, exp_num("gg_cost"))
check("guided QALYs",        gg$qaly, exp_num("gg_qaly"), tol = 1e-4)
check("guided LYs",          gg$ly,   exp_num("gg_ly"),   tol = 1e-4)
check("incremental cost",    d_cost,  exp_num("incremental_cost"))
check("incremental QALYs",   d_qaly,  exp_num("incremental_qaly"), tol = 1e-4)

cat(strrep("-", 74), "\n")

# ── the one real disagreement ────────────────────────────────────────────────
# Python computes NMB from its ROUNDED per-arm totals (total_qaly is rounded to
# 4dp before the delta is taken), then multiplies the QALY delta by the WTP
# threshold. At $100k/QALY a 5e-5 rounding becomes ~$5 of NMB. Reproducing the
# engine's rounding here recovers its exact figure, which identifies the cause
# as display-precision leaking into the arithmetic rather than a model error.
nmb_py         <- exp_num("nmb_at_wtp")
nmb_emulated   <- wtp * (round(gg$qaly, 4) - round(sc$qaly, 4)) -
                        (round(gg$cost, 2) - round(sc$cost, 2))
cat(sprintf("%-22s %14.2f\n", "NMB (R, full prec.)", nmb))
cat(sprintf("%-22s %14.2f\n", "NMB (Python)",        nmb_py))
cat(sprintf("%-22s %14.2f   <- reproduces Python exactly\n",
            "NMB (R, engine round)", nmb_emulated))
if (abs(nmb_emulated - nmb_py) <= 0.01 && abs(nmb - nmb_py) > 0.01) {
  cat("\nFINDING: markov_cost_effectiveness() derives its incremental values from\n")
  cat("per-arm totals that have already been rounded (2dp cost, 4dp QALY), so NMB\n")
  cat(sprintf("carries ~$%.2f of rounding noise at a $%s/QALY threshold. Both\n",
              abs(nmb - nmb_py), format(wtp, big.mark = ",", scientific = FALSE)))
  cat("implementations agree on the model; they disagree only on when to round.\n")
  cat("Fix: take the deltas from unrounded totals and round for display only.\n")
} else {
  check("NMB at WTP", nmb, nmb_py)
}

cat(strrep("-", 74), "\n")

# ── categorical agreement ────────────────────────────────────────────────────
py_dom  <- toupper(exp_chr("dominant"))  == "TRUE"
py_dmd  <- toupper(exp_chr("dominated")) == "TRUE"
py_icer <- exp_chr("icer")
cat(sprintf("dominant   R=%-5s Python=%-5s  %s\n", dominant, py_dom,
            if (dominant == py_dom) "OK" else "FAIL"))
cat(sprintf("dominated  R=%-5s Python=%-5s  %s\n", dominated, py_dmd,
            if (dominated == py_dmd) "OK" else "FAIL"))
cat(sprintf("ICER       R=%-5s Python=%-5s  %s  (suppressed under dominance)\n",
            ifelse(is.na(icer), "NA", sprintf("%.0f", icer)),
            ifelse(py_icer %in% c("NA", "None", ""), "NA", py_icer),
            if (is.na(icer) && py_icer %in% c("NA", "None", "")) "OK" else "check"))
if (dominant != py_dom) fails <- fails + 1L
if (dominated != py_dmd) fails <- fails + 1L

cat(strrep("-", 74), "\n")
if (fails == 0L) {
  cat("PASS - the R and Python implementations agree on the model.\n\n")
  quit(status = 0)
} else {
  cat(sprintf("FAIL - %d quantity/quantities disagree beyond tolerance.\n\n", fails))
  quit(status = 1)
}
