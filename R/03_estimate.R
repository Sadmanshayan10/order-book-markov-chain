# 03_estimate.R
# Purpose: (1) Split trading days into train (first ~60%) and test (rest).
#          (2) Estimate alpha, beta on train:
#                alpha-hat = (n_BA + n_AA) / (2n),  beta-hat = (n_BR + n_AR) / (2n),
#              n = number of BA/BR/AA/AR events, both
#                - uncorrected: events as classified in 02_events.R, and
#                - corrected: plus one BR (AR) for each PX interval in which
#                  the bid (ask) emptied, i.e. the best bid fell (best ask
#                  rose). Depletions are otherwise invisible to beta-hat.
#              Also lambda-hat, mu-hat per second and removal counts by
#              queue size (for the theta extension).
#          (3) Build evaluation samples on a 1-second grid where the spread
#              is one tick: b, a (units of L, capped at K), I = b / (b + a),
#              Y = 1 if the NEXT mid-price change is up, else 0.
#          Only intervals that START with a one-tick spread are used for
#          estimation, matching the samples and the model's assumption.
# Inputs:  data/processed/events_<TICKER>.csv  (02_events.R)
#          data/processed/lob1_<TICKER>.rds    (01_parse.R; for PX depletions,
#                                               queue sizes and the 1 s grid)
# Outputs: data/processed/samples_<TICKER>.csv
#          output/params_<TICKER>.rds

suppressMessages(library(data.table))

TICKER    <- "MSFT"
TICK      <- 100
L         <- 100          # shares per unit (notes/decisions.md, 2026-11-13)
K         <- 155          # cap: 99th pct queue size / L
TRAIN_FR  <- 0.6
WIN_FROM  <- 10.0 * 3600 * 1000
WIN_TO    <- 15.5 * 3600 * 1000

events <- fread(sprintf("data/processed/events_%s.csv", TICKER))
book   <- readRDS(sprintf("data/processed/lob1_%s.rds", TICKER))
events[, day := as.IDate(day)]

# ---- Train / test split by day ------------------------------------------------
# Whole days, in time order, never at random. 18 days -> 11 train, 7 test.
all_days   <- sort(unique(book$day))
n_train    <- round(TRAIN_FR * length(all_days))
train_days <- all_days[seq_len(n_train)]
test_days  <- setdiff(all_days, train_days)
cat(sprintf("Days: %d train (%s to %s), %d test (%s to %s)\n",
            length(train_days), min(train_days), max(train_days),
            length(test_days), min(as.IDate(test_days)), max(as.IDate(test_days))))

# ---- Interval table from the snapshots ---------------------------------------
# Same pairing as 02_events.R: each snapshot with the one before it, same day.
# Needed here for the start-of-interval spread and for PX depletions.
book[, c("p_bid_px", "p_ask_px", "p_bid_sz", "p_ask_sz") :=
         shift(.(bid_px, ask_px, bid_sz, ask_sz)), by = day]
iv <- book[!is.na(p_bid_px) & day %in% train_days]
iv[, start_1tick := (p_ask_px - p_bid_px) == TICK]

# ---- Event counts on train ---------------------------------------------------
# For BA/BR/AA/AR rows prices did not change, so the end-of-interval
# spread_ticks in the events file is also the start-of-interval spread.
tr <- events[day %in% train_days & ev != "PX" & spread_ticks == 1]
cnt <- tr[, .N, by = ev]
n_ev <- setNames(cnt$N, cnt$ev)[c("BA", "BR", "AA", "AR")]

# PX intervals starting at one tick: which side emptied?
#   best bid fell -> bid queue depleted (a BR the event table misses)
#   best ask rose -> ask queue depleted (an AR it misses)
#   best bid rose / best ask fell -> from a one-tick spread this can only happen
#     after the OTHER side emptied in the same interval: the whole quote shifts
#     a tick. It is the follow-on of a depletion, not a separate event.
px <- iv[start_1tick & (bid_px != p_bid_px | ask_px != p_ask_px)]
dep_bid <- px[, sum(bid_px < p_bid_px)]
dep_ask <- px[, sum(ask_px > p_ask_px)]
improve <- px[, sum(bid_px > p_bid_px | ask_px < p_ask_px)]

# ---- alpha-hat, beta-hat with 95% CI ------------------------------------------
# Under the model each event is i.i.d. with probabilities (alpha, beta, alpha,
# beta), so 2 alpha-hat = (n_BA + n_AA) / n is a sample proportion p-hat and
#   Var(alpha-hat) = p (1 - p) / (4 n) = alpha (1/2 - alpha) / n.
# Consecutive events are dependent, so this interval is too narrow anyway.
est <- function(n_add, n_rem) {
  n <- n_add + n_rem
  a <- n_add / (2 * n)
  se <- sqrt(a * (0.5 - a) / n)
  list(alpha = a, beta = 0.5 - a, alpha_ci = a + c(-1, 1) * 1.96 * se, n = n)
}
n_add <- n_ev[["BA"]] + n_ev[["AA"]]
n_rem <- n_ev[["BR"]] + n_ev[["AR"]]
e_raw  <- est(n_add, n_rem)
e_corr <- est(n_add, n_rem + dep_bid + dep_ask)

# Rates per second (rate MLE = count / time), per queue:
# lambda = additions per queue per second, mu = removals per queue per second.
T_sec <- iv[, sum(start_1tick)] * 0.1
rates <- function(n_add, n_rem) c(lambda = n_add / (2 * T_sec), mu = n_rem / (2 * T_sec))
r_raw  <- rates(n_add, n_rem)
r_corr <- rates(n_add, n_rem + dep_bid + dep_ask)

cat("\n==== Train event counts (intervals starting at one tick) ====\n")
print(n_ev)
cat(sprintf("PX intervals: %d  (bid emptied %d, ask emptied %d; other side followed, quote shifted a tick: %d)\n",
            nrow(px), dep_bid, dep_ask, improve))
cat(sprintf("Train time at one tick: %.0f s\n", T_sec))

cat("\n==== Estimates ====\n")
show <- function(lab, e, r) cat(sprintf(
  "%-12s alpha %.4f [%.4f, %.4f]  beta %.4f  rho = alpha/beta %.3f  lambda %.2f/s  mu %.2f/s  n = %d\n",
  lab, e$alpha, e$alpha_ci[1], e$alpha_ci[2], e$beta, e$alpha / e$beta, r[["lambda"]], r[["mu"]], e$n))
show("uncorrected", e_raw, r_raw)
show("corrected", e_corr, r_corr)

# Per-day alpha-hat (both versions): a first look at assumption A2.
by_day_ev <- events[day %in% train_days & ev != "PX" & spread_ticks == 1,
                    .(add = sum(ev %in% c("BA", "AA")), rem = sum(ev %in% c("BR", "AR"))), by = day]
by_day_dp <- px[, .(dep = sum(bid_px < p_bid_px) + sum(ask_px > p_ask_px)), by = day]
by_day <- by_day_dp[by_day_ev, on = "day"]
by_day[is.na(dep), dep := 0]
by_day[, `:=`(alpha_raw = add / (2 * (add + rem)), alpha_corr = add / (2 * (add + rem + dep)))]
cat("\nPer-day alpha-hat on train:\n")
print(by_day[, .(day, alpha_raw = round(alpha_raw, 4), alpha_corr = round(alpha_corr, 4))])

# ---- Removals by queue size (input for the theta extension) -------------------
# For each side and each no-PX interval starting at one tick: queue size at the
# START (units, capped at K, bucketed by 10) and whether that side lost shares.
# removals / intervals per bucket estimates how removal frequency grows with
# queue size.
np <- iv[start_1tick & bid_px == p_bid_px & ask_px == p_ask_px]
to_units <- function(sh) pmin(ceiling(sh / L), K)
rbs <- rbindlist(list(
  np[, .(size = to_units(p_bid_sz), removal = bid_sz < p_bid_sz)],
  np[, .(size = to_units(p_ask_sz), removal = ask_sz < p_ask_sz)]))
removals_by_size <- rbs[, .(intervals = .N, removals = sum(removal)),
                        by = .(size_bucket = 10 * ((size - 1) %/% 10) + 1)][order(size_bucket)]
removals_by_size[, rate := removals / intervals]

# ---- Samples on the 1-second grid ---------------------------------------------
# Y needs the next mid-price change, which can be up to minutes later, so we
# look ahead over ALL snapshots of the day (not just the 1 s grid).
#   chg_idx: positions where the mid differs from the previous snapshot.
#   For a sample at position k, the next change is the first chg_idx > k;
#   Y = 1 if the mid there is above the mid at k.
# Samples with no later change that day get Y = NA and are dropped.
book[, mid2 := bid_px + ask_px]  # 2 x mid, in integer price units
samples <- book[, {
  k_all   <- seq_len(.N)
  chg_idx <- which(mid2 != shift(mid2))
  grid    <- which(time %% 1000 == 0 & spread_ticks == 1)
  nxt     <- chg_idx[findInterval(grid, chg_idx) + 1]   # first change after each grid point
  .(t = time[grid], b_sh = bid_sz[grid], a_sh = ask_sz[grid],
    Y = fifelse(is.na(nxt), NA_integer_, as.integer(mid2[nxt] > mid2[grid])),
    wait_ms = time[nxt] - time[grid])
}, by = day]
n_grid <- nrow(samples)
samples <- samples[!is.na(Y)]
samples[, `:=`(b = to_units(b_sh), a = to_units(a_sh))]
samples[, I := b / (b + a)]
samples[, split := fifelse(day %in% train_days, "train", "test")]
setcolorder(samples, c("day", "t", "b", "a", "I", "Y", "split", "b_sh", "a_sh", "wait_ms"))

cat(sprintf("\n==== Samples ====\n1 s grid points at one tick: %d; dropped %d with no later mid change; kept %d\n",
            n_grid, n_grid - nrow(samples), nrow(samples)))
print(samples[, .(n = .N, share_up = round(mean(Y), 4),
                  capped_b = sum(b_sh > K * L), capped_a = sum(a_sh > K * L),
                  med_wait_s = median(wait_ms) / 1000), by = split])
cat("Y vs naive imbalance on train (share up by I bucket):\n")
print(samples[split == "train", .(n = .N, share_up = round(mean(Y), 3)),
              by = .(I_bucket = cut(I, seq(0, 1, 0.1), include.lowest = TRUE))][order(I_bucket)])

# ---- Save ---------------------------------------------------------------------
fwrite(samples, sprintf("data/processed/samples_%s.csv", TICKER))
params <- list(
  ticker = TICKER, L = L, K = K,
  alpha = e_raw$alpha, beta = e_raw$beta, alpha_ci = e_raw$alpha_ci,
  alpha_corr = e_corr$alpha, beta_corr = e_corr$beta, alpha_ci_corr = e_corr$alpha_ci,
  n_events = e_raw$n, n_depletions = c(bid = dep_bid, ask = dep_ask),
  lambda = r_raw[["lambda"]], mu = r_raw[["mu"]],
  lambda_corr = r_corr[["lambda"]], mu_corr = r_corr[["mu"]],
  window = c("10:00", "15:30"),
  train_days = train_days, test_days = as.IDate(test_days),
  removals_by_size = removals_by_size)
dir.create("output", showWarnings = FALSE)
saveRDS(params, sprintf("output/params_%s.rds", TICKER))
cat(sprintf("\nSaved data/processed/samples_%s.csv (%d rows) and output/params_%s.rds\n",
            TICKER, nrow(samples), TICKER))
