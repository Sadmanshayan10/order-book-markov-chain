# 02_events.R
# Purpose: Turn 100 ms level-1 snapshots into the model's four event types.
#          For each interval between consecutive snapshots of the same day:
#            - best bid or best ask price changed           -> PX
#            - otherwise net change in bid size > 0 / < 0   -> BA / BR
#            - otherwise net change in ask size > 0 / < 0   -> AA / AR
#          One interval can give a bid event AND an ask event (two rows).
#          LIMITATION: several real events can fall in one 100 ms interval;
#          we only see the NET change, so e.g. +500 then -300 shows up as one
#          BA of 200 shares. Event counts are therefore lower bounds.
#          Removals are cross-checked against the market-order list: a sell
#          MO should shrink the bid, a buy MO the ask.
#          Prints the summary used to choose L and K.
# Inputs:  data/processed/lob1_<TICKER>.rds, data/processed/mo_<TICKER>.rds
# Outputs: data/processed/events_<TICKER>.csv (only once L_SHARES is set)
#          Columns follow the standard event-table layout with these adaptations for
#          snapshot data: `day` added; `type` = number of market orders on
#          that side during the interval (no LOBSTER message types exist);
#          `direction` = +1 bid side, -1 ask side, or the sign of the
#          mid-price move for PX rows.

suppressMessages(library(data.table))

TICKER   <- "MSFT"
TICK     <- 100
L_SHARES <- 100        # shares per unit (confirmed 2026-11-13; see notes/decisions.md)

book <- readRDS(sprintf("data/processed/lob1_%s.rds", TICKER))
mo   <- readRDS(sprintf("data/processed/mo_%s.rds", TICKER))

# ---- Intervals: pair each snapshot with the one before it (same day) ----------
# shift() lags each column by one row within a day, so row m holds the book at
# the end of the interval and the p* columns hold the book at its start.
# The first snapshot of each day (10:00:00.000) has no predecessor and is
# dropped: it only serves as a starting book.
book[, c("p_bid_px", "p_ask_px", "p_bid_sz", "p_ask_sz") :=
         shift(.(bid_px, ask_px, bid_sz, ask_sz)), by = day]
iv <- book[!is.na(p_bid_px)]
iv[, px_change := bid_px != p_bid_px | ask_px != p_ask_px]
iv[, d_bid := bid_sz - p_bid_sz]
iv[, d_ask := ask_sz - p_ask_sz]
iv[, d_mid := ((bid_px + ask_px) - (p_bid_px + p_ask_px)) / 2]

# ---- Market orders per interval -----------------------------------------------
# Interval ending at snapshot time t covers (t - 100, t], so an MO at time s
# would naively belong to the interval ending at ceiling(s / 100) * 100.
# Checked empirically: the book only reacts to an MO in the snapshot
# AFTER that one. With no offset, 20% of sell MOs coincide with the bid
# shrinking or a price move; with a one-bin offset, 92.5% do (buy MOs and the
# ask: 21% vs 93%). So the snapshot stamped t shows the book as of t - 100,
# and each MO is assigned to the interval one bin later.
MO_LAG_MS <- 100
mo[, time := ceiling(time / 100) * 100 + MO_LAG_MS]
mo_iv <- mo[, .(n_mo_buy  = sum(mo_type == -1), vol_mo_buy  = sum(volume[mo_type == -1]),
                n_mo_sell = sum(mo_type ==  1), vol_mo_sell = sum(volume[mo_type ==  1])),
            by = .(day, time)]
iv <- mo_iv[iv, on = .(day, time)]
for (col in c("n_mo_buy", "vol_mo_buy", "n_mo_sell", "vol_mo_sell")) {
  set(iv, which(is.na(iv[[col]])), col, 0)
}

# ---- Long event table ---------------------------------------------------------
# One row per PX interval, plus one row per side whose size changed in a
# no-PX interval. Book columns are the state at the END of the interval.
keep <- c("day", "time", "bid_px", "bid_sz", "ask_px", "ask_sz", "spread_ticks")
ev_px  <- iv[px_change == TRUE,
             c(mget(keep), .(ev = "PX", size = NA_real_, direction = sign(d_mid),
                             type = n_mo_buy + n_mo_sell))]
ev_bid <- iv[px_change == FALSE & d_bid != 0,
             c(mget(keep), .(ev = fifelse(d_bid > 0, "BA", "BR"), size = abs(d_bid),
                             direction = 1, type = n_mo_sell))]
ev_ask <- iv[px_change == FALSE & d_ask != 0,
             c(mget(keep), .(ev = fifelse(d_ask > 0, "AA", "AR"), size = abs(d_ask),
                             direction = -1, type = n_mo_buy))]
events <- rbindlist(list(ev_px, ev_bid, ev_ask))
setorder(events, day, time, -direction)

# ---- Summary: per day and overall --------------------------------------------
# Snapshots are evenly spaced, so "% of snapshots" = "% of time".
q <- function(x, p) as.numeric(quantile(x, p, names = FALSE))
summ <- function(b, i, e) {
  ch <- e[ev != "PX", size]
  data.table(intervals    = nrow(i),
             pct_nochange = 100 * i[, mean(!px_change & d_bid == 0 & d_ask == 0)],
             pct_PX       = 100 * i[, mean(px_change)],
             pct_both     = 100 * i[, mean(!px_change & d_bid != 0 & d_ask != 0)],
             n_events     = sum(e$ev != "PX"),
             pct_1tick    = 100 * b[, mean(spread_ticks == 1)],
             bid_p50 = q(b$bid_sz, .5), bid_p90 = q(b$bid_sz, .9), bid_p99 = q(b$bid_sz, .99),
             ask_p50 = q(b$ask_sz, .5), ask_p90 = q(b$ask_sz, .9), ask_p99 = q(b$ask_sz, .99),
             med_change = median(ch))
}
by_day <- rbindlist(lapply(split(seq_len(nrow(book)), book$day), function(r) {
  d <- book$day[r[1]]
  cbind(day = d, summ(book[r], iv[day == d], events[day == d]))
}))
overall <- cbind(day = "ALL", summ(book, iv, events))

options(width = 200)
cat("\n==== Per day and overall (sizes in shares) ====\n")
print(rbind(by_day[, day := as.character(day)], overall), digits = 3)

cat("\n==== Event counts (overall) ====\n")
print(events[, .N, by = ev][order(ev)])

cat("\n==== Size of one nonzero net change (shares), BA/BR/AA/AR pooled ====\n")
ch <- events[ev != "PX", size]
print(round(quantile(ch, c(.10, .25, .5, .75, .9, .99))))
cat(sprintf("share that are multiples of 100: %.1f%%\n", 100 * mean(ch %% 100 == 0)))

# ---- Cross-check removals against the market-order list ----------------------
# Only no-PX intervals: a sell MO that empties the bid shows up as PX instead.
cat("\n==== Cross-check against market orders (no-PX intervals) ====\n")
x <- iv[px_change == FALSE]
cat(sprintf("Intervals with a sell MO: %d. Bid net change <0: %.1f%%, =0: %.1f%%, >0: %.1f%%\n",
            x[n_mo_sell > 0, .N], 100 * x[n_mo_sell > 0, mean(d_bid < 0)],
            100 * x[n_mo_sell > 0, mean(d_bid == 0)], 100 * x[n_mo_sell > 0, mean(d_bid > 0)]))
cat(sprintf("Intervals with a buy MO:  %d. Ask net change <0: %.1f%%, =0: %.1f%%, >0: %.1f%%\n",
            x[n_mo_buy > 0, .N], 100 * x[n_mo_buy > 0, mean(d_ask < 0)],
            100 * x[n_mo_buy > 0, mean(d_ask == 0)], 100 * x[n_mo_buy > 0, mean(d_ask > 0)]))
cat(sprintf("BR events with a sell MO in the interval: %.1f%%  (rest: cancellations only)\n",
            100 * events[ev == "BR", mean(type > 0)]))
cat(sprintf("AR events with a buy MO in the interval:  %.1f%%\n",
            100 * events[ev == "AR", mean(type > 0)]))
cat(sprintf("Market orders that fall in a PX interval: %.1f%%\n",
            100 * iv[, sum((n_mo_buy + n_mo_sell)[px_change]) / sum(n_mo_buy + n_mo_sell)]))
cat(sprintf("Sell-MO intervals where bid fell by at least the MO volume: %.1f%%\n",
            100 * x[n_mo_sell > 0, mean(-d_bid >= vol_mo_sell)]))

# ---- Save (needs L) -----------------------------------------------------------
if (is.na(L_SHARES)) {
  cat("\nL_SHARES not set yet: review the summary, set L, then rerun to save events.\n")
} else {
  events[, units := fifelse(ev == "PX", NA_real_, pmax(1, round(size / L_SHARES)))]
  setcolorder(events, c("day", "time", "type", "direction", "size", "bid_px", "bid_sz",
                        "ask_px", "ask_sz", "spread_ticks", "ev", "units"))
  fwrite(events, sprintf("data/processed/events_%s.csv", TICKER))
  cat(sprintf("\nSaved %d events to data/processed/events_%s.csv\n", nrow(events), TICKER))
}
