# 06_simulation.R
# Purpose: Simulation study: how much prediction error
#          does each broken assumption cause? Small version: two scenarios,
#          compared on the reliability part of the Brier score.
#            1. Sanity run: at the real test states, Y ~ Bernoulli(h(b, a)).
#               This is exactly the fitted chain's outcome distribution, so
#               any reliability left is sampling noise (repeated N_SANITY
#               times). Run at L = 100 (baseline) and L = 1000.
#            2. A4 (one unit per event) broken, at L = 1000: simulate the two
#               queues in SHARES from each test state. Event types are drawn
#               with probabilities (alpha, beta, alpha, beta) as in the model;
#               sizes are drawn from the real sizes (train days, one-tick
#               spread; additions and removals separately). Removal sizes
#               include the depletion chunks: in each PX interval where a side
#               emptied, the shares that side held at the start of the
#               interval (the same depletions 03_estimate.R counts). A queue
#               empties when a removal is at least its size. Queues are capped
#               at K x L shares, as in the chain. Score the BASELINE h at
#               L = 1000 on the simulated outcomes. Paths still running after
#               MAX_STEPS events are dropped and their share reported.
#          L = 1000 (K = 16, the 04_validate.R sensitivity setting) keeps queues at
#          a few units, so h is not pinned near 0.5 and paths end quickly.
#          Scenarios A1-A3 are not run in this version.
#          Reliability: Murphy decomposition over 10 equal-count bins of the
#          predictions, Brier ~ REL - RES + UNC (same bins as the calibration
#          plot in 04_validate.R).
# Inputs:  output/params_MSFT.rds, data/processed/samples_MSFT.csv,
#          data/processed/events_MSFT.csv, data/processed/lob1_MSFT.rds,
#          R/model.R
# Outputs: output/simulation_MSFT.csv  Brier, REL, RES, UNC per scenario, and
#          for the A4 rows the mean change of one queue per event (shares)

suppressMessages(library(data.table))
source("R/model.R")

TICKER    <- "MSFT"
SEED      <- 348
N_SANITY  <- 20           # repeats of the sanity run
L_A4      <- 1000         # shares per unit in the A4 scenario
MAX_STEPS <- 20000        # events per path before it is dropped

params  <- readRDS(sprintf("output/params_%s.rds", TICKER))
samples <- fread(sprintf("data/processed/samples_%s.csv", TICKER))
events  <- fread(sprintf("data/processed/events_%s.csv", TICKER))
events[, day := as.IDate(day)]
test   <- samples[split == "test"]
alpha  <- params$alpha_corr
beta   <- params$beta_corr
cap_sh <- params$K * params$L                    # 15,500 shares, as in 04_validate.R
K_A4   <- ceiling(cap_sh / L_A4)                 # 16

# ---- Murphy decomposition ------------------------------------------------------
# murphy(p, y): Brier score and its reliability, resolution and uncertainty
# parts over equal-count bins of p. REL - RES + UNC equals the Brier score up to
# the spread of p within bins.
murphy <- function(p, y, bins = 10) {
  br  <- unique(quantile(p, seq(0, 1, length.out = bins + 1), names = FALSE))
  bin <- if (length(br) > 1) cut(p, br, include.lowest = TRUE, labels = FALSE) else rep(1L, length(p))
  g   <- data.table(bin, p, y)[, .(n = .N, pbar = mean(p), ybar = mean(y)), by = bin]
  M   <- length(y)
  c(brier = mean((p - y)^2),
    rel   = sum(g$n * (g$pbar - g$ybar)^2) / M,
    res   = sum(g$n * (g$ybar - mean(y))^2) / M,
    unc   = mean(y) * (1 - mean(y)))
}

# Known answers: a constant forecast has REL = (p - ybar)^2 and RES = 0.
local({
  y <- c(0, 1, 1, 0, 1)
  m <- murphy(rep(0.5, 5), y)
  stopifnot(abs(m[["rel"]] - 0.01) < 1e-12, m[["res"]] == 0,
            abs(m[["brier"]] - 0.25) < 1e-12, abs(m[["unc"]] - 0.24) < 1e-12)
})

# h at each test state for unit size L (cap K x L shares, as in 04_validate.R).
h_test <- function(L) {
  K <- ceiling(cap_sh / L)
  H <- h_table(K, alpha, beta)
  H[cbind(pmin(ceiling(test$b_sh / L), K), pmin(ceiling(test$a_sh / L), K))]
}
p100  <- h_test(params$L)
p1000 <- h_test(L_A4)

# ---- 1. Sanity run ---------------------------------------------------------------
set.seed(SEED)
sanity <- rbindlist(lapply(c(params$L, L_A4), function(L) {
  p <- if (L == params$L) p100 else p1000
  rbindlist(lapply(seq_len(N_SANITY), function(r)
    as.data.table(as.list(murphy(p, rbinom(length(p), 1, p))))))[, L := L]
}))
cat(sprintf("\n==== Sanity run: Y ~ Bernoulli(h), %d repeats, %d test states ====\n",
            N_SANITY, nrow(test)))
print(sanity[, .(rel_mean = mean(rel), rel_max = max(rel), brier_mean = mean(brier)), by = L],
      digits = 4)

# ---- 2. A4: real event sizes, L = 1000 ---------------------------------------------
# Sizes in shares from the train days, one-tick spread (same events as
# 05_assumptions.R uses for the A4 test).
sz <- events[day %in% params$train_days & ev != "PX" & spread_ticks == 1]
add_sz <- abs(sz[ev %in% c("BA", "AA"), size])
rem_nonpx <- abs(sz[ev %in% c("BR", "AR"), size])

# Depletion chunks. A removal that empties a queue is recorded as PX, not BR/AR,
# so without these the removal sizes would be cut off at "less than the whole
# queue". Same definition as 03_estimate.R: in a PX interval that starts at one
# tick, the best bid fell (bid emptied) or the best ask rose (ask emptied); the
# chunk is what that side held at the start of the interval.
book <- readRDS(sprintf("data/processed/lob1_%s.rds", TICKER))[day %in% params$train_days]
book[, c("p_bid_px", "p_ask_px", "p_bid_sz", "p_ask_sz") :=
         shift(.(bid_px, ask_px, bid_sz, ask_sz)), by = day]
px <- book[!is.na(p_bid_px) & (p_ask_px - p_bid_px) == 100 &
             (bid_px != p_bid_px | ask_px != p_ask_px)]
dep_sz <- c(px[bid_px < p_bid_px, p_bid_sz], px[ask_px > p_ask_px, p_ask_sz])
rem_sz <- c(rem_nonpx, dep_sz)
stopifnot(length(dep_sz) == sum(params$n_depletions))   # matches 03_estimate.R

mean_change <- alpha * mean(add_sz) - beta * mean(rem_sz)
cat(sprintf("\nA4 sizes (shares): additions n %d, mean %.0f, median %.0f\n",
            length(add_sz), mean(add_sz), median(add_sz)))
cat(sprintf("  removals without depletions: n %d, mean %.0f, median %.0f\n",
            length(rem_nonpx), mean(rem_nonpx), median(rem_nonpx)))
cat(sprintf("  depletion chunks:            n %d, mean %.0f, median %.0f\n",
            length(dep_sz), mean(dep_sz), median(dep_sz)))
cat(sprintf("  all removals (used):         n %d, mean %.0f, median %.0f\n",
            length(rem_sz), mean(rem_sz), median(rem_sz)))
cat(sprintf("Mean change of one queue per event: %+.1f shares (without depletions: %+.1f)\n",
            mean_change, alpha * mean(add_sz) - beta * mean(rem_nonpx)))

# All paths move together, one event per step; finished paths drop out.
cap  <- K_A4 * L_A4
b    <- pmin(test$b_sh, cap)
a    <- pmin(test$a_sh, cap)
Ysim <- rep(NA_integer_, nrow(test))
steps <- rep(NA_integer_, nrow(test))
live <- seq_len(nrow(test))
for (s in seq_len(MAX_STEPS)) {
  n   <- length(live)
  if (n == 0) break
  typ <- sample.int(4, n, replace = TRUE, prob = c(alpha, beta, alpha, beta))  # BA BR AA AR
  add <- typ %in% c(1, 3)
  amt <- numeric(n)
  amt[add]  <- add_sz[sample.int(length(add_sz), sum(add),  replace = TRUE)]
  amt[!add] <- rem_sz[sample.int(length(rem_sz), sum(!add), replace = TRUE)]
  bid <- typ <= 2
  i <- live[bid];  b[i] <- pmin(b[i] + ifelse(add[bid],  amt[bid],  -amt[bid]),  cap)
  i <- live[!bid]; a[i] <- pmin(a[i] + ifelse(add[!bid], amt[!bid], -amt[!bid]), cap)
  down <- b[live] <= 0
  up   <- a[live] <= 0 & !down
  Ysim[live[down]] <- 0L
  Ysim[live[up]]   <- 1L
  steps[live[down | up]] <- s
  live <- live[!(down | up)]
}
done <- !is.na(Ysim)
cat(sprintf("A4 paths: %d of %d absorbed within %d events (%.2f%% unabsorbed); median %d events\n",
            sum(done), length(done), MAX_STEPS, 100 * mean(!done), median(steps, na.rm = TRUE)))

# ---- Results -------------------------------------------------------------------
row <- function(scenario, L, m, extra = list())
  as.data.table(c(list(scenario = scenario, L = L), as.list(m), extra))
res <- rbind(
  row("real data", params$L, murphy(p100, test$Y)),
  row("real data", L_A4,     murphy(p1000, test$Y)),
  sanity[, .(scenario = "sanity", L, brier = mean(brier), rel = mean(rel),
             res = mean(res), unc = mean(unc)), by = L][, .(scenario, L, brier, rel, res, unc)],
  row("A4 real sizes", L_A4, murphy(p1000[done], Ysim[done])),
  # Simulated paths are in shares, so the same outcomes also score the L = 100 h.
  row("A4 real sizes", params$L, murphy(p100[done], Ysim[done])),
  fill = TRUE)
res[, unabsorbed := c(NA, NA, NA, NA, mean(!done), mean(!done))]
res[, mean_change_sh := c(NA, NA, NA, NA, mean_change, mean_change)]
cat("\n==== Baseline h scored on each scenario (sanity: mean of repeats) ====\n")
print(res, digits = 4)
cat(sprintf("Share up: real %.4f, A4 simulated %.4f\n", mean(test$Y), mean(Ysim[done])))

fwrite(res, sprintf("output/simulation_%s.csv", TICKER))
cat(sprintf("\nWrote output/simulation_%s.csv\n", TICKER))
