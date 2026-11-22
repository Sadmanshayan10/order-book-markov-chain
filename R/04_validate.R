# 04_validate.R
# Purpose: First predictive result on the TEST days.
#          Predictors of Y = 1{next mid-price move is up}:
#            null      p = 0.5
#            naive     p = I = b / (b + a)
#            logistic  logit p = c0 + c1 I, fitted on train (Gould & Bonart)
#            chain     p = h(b, a) at alpha-hat, beta-hat (baseline: corrected
#                      for depletions; uncorrected reported alongside)
#          Scored with Brier, log loss, AUC. Plus sensitivity of the chain to
#          L (100, 500, 1000) and to K (100, 155, 250).
#          t(b, a) is NOT used: with alpha > beta it is driven by the cap K
#          (notes/decisions.md, 2026-11-13).
# Inputs:  data/processed/samples_<TICKER>.csv, output/params_<TICKER>.rds,
#          R/model.R
# Outputs: output/results_<TICKER>.csv      predictor comparison (test days)
#          output/sensitivity_<TICKER>.csv  chain under other L and K
#          output/brier_diff_<TICKER>.csv   Brier differences with day-bootstrap CIs
#          output/heatmap_h_<TICKER>.png    h(b, a) at the baseline estimates
#          output/calibration_<TICKER>.png  calibration plot (test days)

suppressMessages(library(data.table))
source("R/model.R")

TICKER <- "MSFT"
samples <- fread(sprintf("data/processed/samples_%s.csv", TICKER))
params  <- readRDS(sprintf("output/params_%s.rds", TICKER))

# ---- Scoring functions --------------------------------------------------------
# score(p, y): named vector of Brier score, log loss and AUC.
#   Brier    mean (p - y)^2; always predicting 0.5 scores exactly 0.25.
#   Log loss -mean(y log p + (1 - y) log(1 - p)); p clipped away from 0 and 1
#            so a single confident miss cannot make it infinite.
#   AUC      probability that a random up-move gets a higher p than a random
#            down-move (ties count 1/2), via the Mann-Whitney rank formula.
score <- function(p, y) {
  pc <- pmin(pmax(p, 1e-15), 1 - 1e-15)
  n1 <- as.numeric(sum(y == 1))  # numeric: n1 * n0 overflows integer range
  n0 <- as.numeric(sum(y == 0))
  auc <- (sum(rank(p)[y == 1]) - n1 * (n1 + 1) / 2) / (n1 * n0)
  c(brier    = mean((p - y)^2),
    log_loss = -mean(y * log(pc) + (1 - y) * log(1 - pc)),
    auc      = auc)
}

# calib_table(p, y, bins): bins by predicted probability (equal-count bins),
# mean prediction vs observed share of up-moves in each.
calib_table <- function(p, y, bins = 10) {
  br  <- unique(quantile(p, seq(0, 1, length.out = bins + 1), names = FALSE))
  bin <- if (length(br) > 1) cut(p, br, include.lowest = TRUE, labels = FALSE) else rep(1L, length(p))
  data.table(bin = bin, p = p, y = y)[
    , .(mean_pred = mean(p), obs_freq = mean(y), count = .N), keyby = bin]
}

# Known answers: stop if any fails.
local({
  y <- c(0, 1, 1, 0, 1)
  s <- score(rep(0.5, 5), y)
  stopifnot(abs(s[["brier"]] - 0.25) < 1e-12,
            abs(s[["log_loss"]] - log(2)) < 1e-12,
            abs(s[["auc"]] - 0.5) < 1e-12)
  s <- score(y, y)
  stopifnot(s[["brier"]] == 0, s[["auc"]] == 1)
  stopifnot(abs(score(c(.2, .4, .6, .8), c(0, 1, 0, 1))[["auc"]] - 0.75) < 1e-12)
})
cat("score() checks passed: p = 0.5 gives Brier 0.25, log loss log 2, AUC 0.5\n")

# ---- Predictions --------------------------------------------------------------
train <- samples[split == "train"]
test  <- samples[split == "test"]
K <- params$K

# h at each sample's state: index the K x K table with (b, a) pairs.
h_at <- function(H, b, a) H[cbind(b, a)]
H_corr <- h_table(K, params$alpha_corr, params$beta_corr)
H_raw  <- h_table(K, params$alpha,      params$beta)

logit_fit <- glm(Y ~ I, family = binomial, data = train)

preds <- list(
  "null (0.5)"              = rep(0.5, nrow(test)),
  "naive I"                 = test$I,
  "logistic on I (train)"   = as.numeric(predict(logit_fit, test, type = "response")),
  "chain, corrected alpha"  = h_at(H_corr, test$b, test$a),
  "chain, uncorrected alpha"= h_at(H_raw,  test$b, test$a))

results <- rbindlist(lapply(names(preds), function(nm) {
  p <- preds[[nm]]
  as.data.table(as.list(c(score(p, test$Y),
                          p_sd = sd(p), p_min = min(p), p_max = max(p))))[, predictor := nm]
}))
results[, brier_skill := 1 - brier / 0.25]
setcolorder(results, c("predictor", "brier", "brier_skill", "log_loss", "auc", "p_sd", "p_min", "p_max"))

cat(sprintf("\n==== Test days: %d samples, share up %.4f ====\n", nrow(test), mean(test$Y)))
cat(sprintf("Baseline chain: K = %d, L = %d, alpha %.4f, beta %.4f\n",
            K, params$L, params$alpha_corr, params$beta_corr))
cat(sprintf("Logistic fit on train: c0 = %.3f, c1 = %.3f\n", coef(logit_fit)[1], coef(logit_fit)[2]))
print(results, digits = 4)

# ---- Sensitivity: L and K ------------------------------------------------------
# L: rebuild b, a from the share columns; keep the baseline cap in SHARES
#    (K x L = 15,500, the 02_events.R 99th pct) so K = ceiling(15,500 / L). alpha, beta do not depend on L (they count
#    events, not shares), so they stay at the baseline estimates.
# K: baseline L = 100, cap at 100, 155, 250.
cap_sh <- params$K * params$L
chain_score <- function(L, K) {
  b <- pmin(ceiling(test$b_sh / L), K)
  a <- pmin(ceiling(test$a_sh / L), K)
  H <- h_table(K, params$alpha_corr, params$beta_corr)
  p <- h_at(H, b, a)
  list(p = p, row = as.data.table(as.list(c(L = L, K = K, score(p, test$Y),
                                             p_sd = sd(p), p_min = min(p), p_max = max(p)))))
}
sens_L <- lapply(c(100, 500, 1000), function(L) chain_score(L, ceiling(cap_sh / L)))
sens_K <- lapply(c(100, 155, 250), function(K) chain_score(100, K))
sens <- rbindlist(c(lapply(sens_L, `[[`, "row"), lapply(sens_K, `[[`, "row")))
sens[, check := c(rep("L", 3), rep("K", 3))]
sens <- unique(sens, by = c("L", "K"))
cat(sprintf("\n==== Chain sensitivity (cap = %.0f shares) ====\n", cap_sh))
print(sens, digits = 4)
pK <- lapply(sens_K, `[[`, "p")
cat(sprintf("Largest change in any test prediction: K 100 vs 155 = %.2e, K 250 vs 155 = %.2e\n",
            max(abs(pK[[1]] - pK[[2]])), max(abs(pK[[3]] - pK[[2]]))))

# ---- Uncertainty: per-day scores and a day-level bootstrap --------------------
# Samples within a day are strongly dependent (consecutive seconds share the
# same queues and often the same next move), so resampling individual rows
# would give intervals far too narrow. Instead resample whole TEST DAYS with
# replacement (B = 2000) and recompute the Brier difference each time.
# With only 7 test days this is coarse; the per-day table shows the raw spread.
pred_dt <- data.table(day = test$day, Y = test$Y,
                      null = preds[["null (0.5)"]], naive = preds[["naive I"]],
                      logistic = preds[["logistic on I (train)"]],
                      chain = preds[["chain, corrected alpha"]])
# Per-day sum of squared errors and count, so any resample of days can be
# scored by adding rows up.
day_sse <- pred_dt[, .(n = .N,
                       null = sum((null - Y)^2), naive = sum((naive - Y)^2),
                       logistic = sum((logistic - Y)^2), chain = sum((chain - Y)^2)), by = day]
cat("\n==== Brier score by test day ====\n")
print(day_sse[, .(day, n, null = null / n, naive = naive / n,
                  logistic = logistic / n, chain = chain / n)], digits = 4)

set.seed(348)
B <- 2000
boot <- t(replicate(B, {
  d <- day_sse[sample.int(.N, .N, replace = TRUE)]
  n <- sum(d$n)
  c(chain_minus_null     = (sum(d$chain) - sum(d$null)) / n,
    chain_minus_naive    = (sum(d$chain) - sum(d$naive)) / n,
    chain_minus_logistic = (sum(d$chain) - sum(d$logistic)) / n)
}))
full <- with(day_sse, c(chain_minus_null     = (sum(chain) - sum(null)) / sum(n),
                        chain_minus_naive    = (sum(chain) - sum(naive)) / sum(n),
                        chain_minus_logistic = (sum(chain) - sum(logistic)) / sum(n)))
brier_diff <- data.table(comparison = names(full), estimate = full,
                         ci_low  = apply(boot, 2, quantile, 0.025),
                         ci_high = apply(boot, 2, quantile, 0.975),
                         days_chain_better = with(day_sse, c(sum(chain < null), sum(chain < naive),
                                                             sum(chain < logistic))))
cat(sprintf("\n==== Brier difference, chain minus other (negative = chain better), day bootstrap B = %d ====\n", B))
print(brier_diff, digits = 3)

# ---- Save tables --------------------------------------------------------------
fwrite(results,    sprintf("output/results_%s.csv", TICKER))
fwrite(sens,       sprintf("output/sensitivity_%s.csv", TICKER))
fwrite(brier_diff, sprintf("output/brier_diff_%s.csv", TICKER))

# ---- Heatmap of h(b, a) ---------------------------------------------------------
# Plotted up to 125 units (95th pct of the larger queue in the samples); the
# remaining states up to K behave the same way. Diverging palette centred on
# 0.5 so "up" and "down" read symmetrically; contours at 0.3/0.4/0.5/0.6/0.7.
png(sprintf("output/heatmap_h_%s.png", TICKER), width = 1600, height = 1400, res = 220)
show <- 1:125
pal <- hcl.colors(41, "Blue-Red 3")
image(show, show, H_corr[show, show], col = pal, zlim = c(0, 1),
      xlab = "bid queue b (units of 100 shares)", ylab = "ask queue a (units of 100 shares)",
      main = sprintf("h(b, a) = P(next move up), %s, alpha = %.4f, beta = %.4f, K = %d",
                     TICKER, params$alpha_corr, params$beta_corr, K),
      cex.main = 0.8, useRaster = TRUE)
contour(show, show, H_corr[show, show], levels = c(.3, .4, .5, .6, .7), add = TRUE,
        labcex = 0.7, col = "grey20")
legend("topright", c("h near 0: next move down", "h = 0.5", "h near 1: next move up"),
       fill = pal[c(1, 21, 41)], bg = "white", cex = 0.7)
invisible(dev.off())

# ---- Calibration plot (test days) --------------------------------------------
# Each point is one of 10 equal-count bins of the predictions; a well-calibrated
# predictor sits on the 45-degree line. Null is a single point at (0.5, share up).
# Drawn at its printed size (6 x 3.5 in) with a 10 pt base font, so axis and
# legend text are 9-10 pt in the report. No title: the report caption has it.
png(sprintf("output/calibration_%s.png", TICKER), width = 6, height = 3.5, units = "in",
    res = 300, pointsize = 10)
par(mar = c(3.6, 3.6, 0.6, 11), mgp = c(2.3, 0.7, 0), cex.axis = 0.9)
plot(c(0, 1), c(0.2, 0.9), type = "n", xlab = "mean predicted P(up)",
     ylab = "observed share up", las = 1)
abline(0, 1, lty = 2, col = "grey50")
show_pred <- c("naive I", "logistic on I (train)", "chain, corrected alpha")
show_lab  <- c("naive imbalance I", "logistic on I", "chain h(b, a)")
cols <- c("#E69F00", "#0072B2", "#D55E00")
for (i in seq_along(show_pred)) {
  ct <- calib_table(preds[[show_pred[i]]], test$Y, 10)
  lines(ct$mean_pred, ct$obs_freq, type = "b", pch = 19, col = cols[i], lwd = 2)
}
points(0.5, mean(test$Y), pch = 4, cex = 1.5, lwd = 2)
legend(par("usr")[2] + 0.03, par("usr")[4], c(show_lab, "null (0.5)", "perfect calibration"),
       col = c(cols, "black", "grey50"), pch = c(19, 19, 19, 4, NA),
       lty = c(1, 1, 1, NA, 2), lwd = c(2, 2, 2, NA, 1), bty = "n", cex = 0.9, xpd = TRUE)
invisible(dev.off())

cat(sprintf("\nSaved output/results_%s.csv, sensitivity_%s.csv, heatmap_h_%s.png, calibration_%s.png\n",
            TICKER, TICKER, TICKER, TICKER))
