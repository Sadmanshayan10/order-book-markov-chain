# 05_assumptions.R
# Purpose: Test the chain's assumptions on the TRAIN days:
#            A1  event types independent (Markov): 4 x 4 table of consecutive
#                events, chi-square test of independence (9 df)
#            A2  same probabilities all day: alpha-hat by 30-minute bucket,
#                chi-square test of homogeneity
#            A3  same probabilities at every queue size: P(removal | event on
#                that side) by queue-size bucket, chi-square test
#            A4  one unit per event: distribution of |net change| / L
#            Timing (not needed for h): gaps between events vs exponential
#          With ~600,000 events every test rejects; the useful output is the
#          SIZE of each departure (Cramer's V and the estimated probabilities
#          themselves). Snapshot caveat: events are net changes per 100 ms, so
#          A1 cannot see the order of events inside one interval and timing
#          gaps are multiples of 100 ms.
# Inputs:  data/processed/events_<TICKER>.csv, data/processed/lob1_<TICKER>.rds,
#          output/params_<TICKER>.rds
# Outputs: output/assumptions_<TICKER>.csv      one row per test
#          output/assumption_tables_<TICKER>.rds  A1 transition matrix, A2 and
#                                                 A3 tables (inputs for 06)
#          output/assumptions_<TICKER>.png      four-panel figure

suppressMessages(library(data.table))

TICKER <- "MSFT"
TICK   <- 100
params <- readRDS(sprintf("output/params_%s.rds", TICKER))
L <- params$L
K <- params$K
train_days <- params$train_days

events <- fread(sprintf("data/processed/events_%s.csv", TICKER))
events[, day := as.IDate(day)]
events <- events[day %in% train_days]
book <- readRDS(sprintf("data/processed/lob1_%s.rds", TICKER))[day %in% train_days]

# Cramer's V: chi-square rescaled to [0, 1] so departures can be compared
# across tables of different size and sample size. 0 = no association.
cramers_v <- function(tab) {
  ct <- suppressWarnings(chisq.test(tab))
  v  <- sqrt(ct$statistic / (sum(tab) * (min(dim(tab)) - 1)))
  list(chisq = unname(ct$statistic), df = unname(ct$parameter),
       p = ct$p.value, v = unname(v))
}

# ---- Side-event table from the snapshots ---------------------------------------
# Same definitions as 03_estimate.R (intervals starting at one tick), with the
# queue size at the START of the interval, which A3 needs. Depletions (PX
# intervals where a side emptied) are included as removals, as in the
# corrected estimate.
book[, c("p_bid_px", "p_ask_px", "p_bid_sz", "p_ask_sz") :=
         shift(.(bid_px, ask_px, bid_sz, ask_sz)), by = day]
iv <- book[!is.na(p_bid_px) & (p_ask_px - p_bid_px) == TICK]
iv[, px := bid_px != p_bid_px | ask_px != p_ask_px]
side_ev <- rbindlist(list(
  iv[!px & bid_sz != p_bid_sz, .(day, time, side = "bid", start = p_bid_sz, removal = bid_sz < p_bid_sz)],
  iv[!px & ask_sz != p_ask_sz, .(day, time, side = "ask", start = p_ask_sz, removal = ask_sz < p_ask_sz)],
  iv[px & bid_px < p_bid_px,   .(day, time, side = "bid", start = p_bid_sz, removal = TRUE)],
  iv[px & ask_px > p_ask_px,   .(day, time, side = "ask", start = p_ask_sz, removal = TRUE)]))
side_ev[, units := pmin(ceiling(start / L), K)]
alpha_all <- side_ev[, mean(!removal) / 2]
cat(sprintf("Train side events (incl. depletions): %d; alpha-hat %.4f (params: %.4f)\n",
            nrow(side_ev), alpha_all, params$alpha_corr))

# ---- A1: consecutive event types ------------------------------------------------
# Work at the level of 100 ms intervals. A PX resets the queues, so pairs never
# cross a PX or a day boundary (seg increments at each PX). When BOTH sides
# changed in one interval we cannot tell which came first, and any fixed
# ordering (e.g. bid first) would bias the table, so we keep only pairs of
# consecutive event intervals that each contain exactly ONE event.
setorder(events, day, time, -direction)
events[, seg := cumsum(ev == "PX"), by = day]
iv_ev <- events[ev != "PX", .(n = .N, ev = ev[1]), by = .(day, seg, time)]
iv_ev[, `:=`(next_ev = shift(ev, -1), next_n = shift(n, -1)), by = .(day, seg)]
pairs <- iv_ev[!is.na(next_ev) & n == 1 & next_n == 1]
lv <- c("BA", "BR", "AA", "AR")
tab1 <- table(factor(pairs$ev, lv), factor(pairs$next_ev, lv))
a1 <- cramers_v(tab1)
cond1 <- prop.table(tab1, 1)
cat(sprintf("\n==== A1: P(next event | current event), %d pairs ====\n", sum(tab1)))
print(round(cond1, 3))
cat("Marginal (what independence predicts for every row):\n")
print(round(prop.table(colSums(tab1)), 3))

# ---- A2: half-hour buckets -----------------------------------------------------
side_ev[, bucket := pmin((time - 36000000) %/% 1800000, 10)]   # 0 = 10:00-10:30, ..., 10 = 15:00-15:30
tab2 <- side_ev[, table(bucket, removal)]
a2 <- cramers_v(tab2)
by_bucket <- side_ev[, .(n = .N, alpha = mean(!removal) / 2), keyby = bucket]
by_bucket[, `:=`(start = sprintf("%02d:%02d", 10 + (bucket * 30) %/% 60, (bucket * 30) %% 60),
                 ci = 1.96 * sqrt(alpha * (0.5 - alpha) / n))]
cat("\n==== A2: alpha-hat by half-hour (train) ====\n")
print(by_bucket[, .(start, n, alpha = round(alpha, 4), ci = round(ci, 4))])

# ---- A3: queue size ------------------------------------------------------------
# Under the model P(removal | an event on this side) = beta / (alpha + beta),
# the same at every queue size. Buckets of 10 units (1-10, 11-20, ...).
side_ev[, size_bucket := 10 * ((units - 1) %/% 10) + 1]
tab3 <- side_ev[, table(size_bucket, removal)]
a3 <- cramers_v(tab3)
by_size <- side_ev[, .(n = .N, p_removal = mean(removal)), keyby = .(side, size_bucket)]
p_rem_model <- 1 - 2 * alpha_all   # beta / (alpha + beta) = 2 beta
cat(sprintf("\n==== A3: P(removal | event on that side) by queue size (model: %.3f everywhere) ====\n",
            p_rem_model))
print(dcast(by_size, size_bucket ~ side, value.var = "p_removal")[, lapply(.SD, round, 3)])

# ---- A4: units per event ---------------------------------------------------------
# |net change| / L for BA/BR/AA/AR (train, one-tick spread). The model moves
# exactly one unit per event. Depletions: how many units vanished in the
# interval where the queue emptied (the start size).
u <- events[ev != "PX" & spread_ticks == 1, size / L]
dep_units <- iv[px & (bid_px < p_bid_px | ask_px > p_ask_px),
                fifelse(bid_px < p_bid_px, p_bid_sz, p_ask_sz) / L]
a4 <- c(share_1unit = mean(u <= 1.5), share_ge5 = mean(u >= 5), share_ge10 = mean(u >= 10),
        q50 = median(u), q90 = quantile(u, .9, names = FALSE), q99 = quantile(u, .99, names = FALSE),
        dep_median_units = median(dep_units))
cat("\n==== A4: units moved per event (model: always 1) ====\n")
print(round(a4, 3))

# ---- Timing ---------------------------------------------------------------------
# Gaps between consecutive intervals that contain at least one event, within a
# day. Exponential gaps have coefficient of variation sd/mean = 1; bursty
# trading gives CV > 1. Gaps are multiples of 100 ms (snapshot resolution).
ev_times <- unique(events[ev != "PX", .(day, time)])
gaps <- ev_times[, diff(time) / 1000, by = day]$V1
cv <- sd(gaps) / mean(gaps)
cat(sprintf("\n==== Timing: %d gaps, mean %.3f s, CV %.2f (exponential: 1), share of gaps = 100 ms: %.3f ====\n",
            length(gaps), mean(gaps), cv, mean(gaps == 0.1)))

# ---- Summary table and saved inputs for the simulation study ----------------------
summary_tab <- data.table(
  test = c("A1 Markov (consecutive types)", "A2 time of day", "A3 queue size",
           "A4 one unit per event", "Timing (exponential gaps)"),
  statistic = c(a1$chisq, a2$chisq, a3$chisq, a4[["q50"]], cv),
  stat_name = c("chi-square", "chi-square", "chi-square", "median units", "CV of gaps"),
  df = c(a1$df, a2$df, a3$df, NA, NA),
  p_value = c(a1$p, a2$p, a3$p, NA, NA),
  effect = c(a1$v, a2$v, a3$v, 1 - a4[["share_1unit"]], cv - 1),
  effect_name = c("Cramer's V", "Cramer's V", "Cramer's V", "share of events > 1 unit", "CV - 1"))
cat("\n==== Summary ====\n")
print(summary_tab, digits = 3)
fwrite(summary_tab, sprintf("output/assumptions_%s.csv", TICKER))
saveRDS(list(A1_transition = cond1, A1_counts = tab1, A2_by_bucket = by_bucket,
             A3_by_size = by_size, A4_units = a4, timing_cv = cv),
        sprintf("output/assumption_tables_%s.rds", TICKER))

# ---- Figure ----------------------------------------------------------------------
png(sprintf("output/assumptions_%s.png", TICKER), width = 2400, height = 2000, res = 220)
par(mfrow = c(2, 2), mar = c(4.5, 4.5, 3, 1))

# A1: P(next | current) as a table-like heatmap; each cell shows the value,
# with the marginal (what independence predicts) in brackets in the header.
marg <- prop.table(colSums(tab1))
image(1:4, 1:4, t(cond1[4:1, ]), col = hcl.colors(30, "Oranges", rev = TRUE), zlim = c(0, 0.6),
      axes = FALSE, xlab = "next event", ylab = "current event",
      main = "A1: P(next | current); independence = column marginal")
axis(1, 1:4, sprintf("%s (%.2f)", lv, marg), cex.axis = 0.8)
axis(2, 1:4, rev(lv), las = 1)
# t(cond1) runs row by row (current BA first), so x cycles fastest.
text(rep(1:4, times = 4), rep(4:1, each = 4), sprintf("%.2f", as.numeric(t(cond1))), cex = 0.9)

# A2: alpha by half-hour with 95% CI and the pooled estimate
with(by_bucket, {
  plot(bucket, alpha, pch = 19, xaxt = "n", ylim = range(alpha - ci, alpha + ci, 0.25),
       xlab = "half-hour starting", ylab = "alpha-hat", main = "A2: alpha-hat through the day")
  axis(1, at = bucket, labels = start, cex.axis = 0.7)
  segments(bucket, alpha - ci, bucket, alpha + ci)
  abline(h = alpha_all, lty = 2, col = "#0072B2")
  abline(h = 0.25, lty = 3, col = "grey50")
})
legend("bottomright", c("pooled", "alpha = beta"), lty = c(2, 3), col = c("#0072B2", "grey50"),
       bty = "n", cex = 0.8)

# A3: removal probability by queue size, both sides
plot(range(by_size$size_bucket), range(by_size$p_removal, p_rem_model), type = "n",
     xlab = "queue size at start (units, bucket start)", ylab = "P(removal | event)",
     main = "A3: removal probability by queue size")
for (s in c("bid", "ask")) with(by_size[side == s],
  lines(size_bucket, p_removal, type = "b", pch = 19,
        col = if (s == "bid") "#0072B2" else "#D55E00"))
abline(h = p_rem_model, lty = 2)
legend("topright", c("bid", "ask", "model (constant)"), col = c("#0072B2", "#D55E00", "black"),
       lty = c(1, 1, 2), pch = c(19, 19, NA), bty = "n", cex = 0.8)

# Timing: QQ plot of gaps against the exponential with the same mean
qq_p <- ppoints(2000)
plot(qexp(qq_p, 1 / mean(gaps)), quantile(gaps, qq_p), pch = 20, cex = 0.5,
     xlab = "exponential quantile (s)", ylab = "observed gap quantile (s)",
     main = sprintf("Timing: gaps vs exponential (CV = %.2f)", cv), log = "xy")
abline(0, 1, lty = 2)
invisible(dev.off())

cat(sprintf("\nSaved output/assumptions_%s.csv, assumption_tables_%s.rds, assumptions_%s.png\n",
            TICKER, TICKER, TICKER))
