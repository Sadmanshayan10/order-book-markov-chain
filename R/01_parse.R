# 01_parse.R
# Purpose: Read the Cartea-Jaimungal-Penalva Nov 2014 NASDAQ files for one
#          ticker and keep what the model needs:
#            - the level-1 book (best bid/ask price and size) every 100 ms,
#              10:00-15:30 only;
#            - the market-order list, same window.
#          Each raw file is one day, a MATLAB v5 struct LOB with fields
#            EventTime  ms after midnight, one snapshot every 100 ms
#            BuyPrice, SellPrice, BuyVolume, SellVolume   20 levels, $ x 10000
#            NumberMO, VolumeMO   market orders per 100 ms (col 1 buy, col 2 sell)
#            MO         one row per market order: time, bid, ask, bid vol,
#                       ask vol, avg price, volume, type (-1 buy, +1 sell), ...
# Inputs:  data/raw/<TICKER>_YYYYMMDD.mat
# Outputs: data/processed/lob1_<TICKER>.rds  level-1 snapshots, all days
#          data/processed/mo_<TICKER>.rds    market orders, all days

suppressMessages({
  library(R.matlab)
  library(data.table)
})

TICKER   <- "MSFT"
TICK     <- 100                 # $0.01 in the file's units ($ x 10000)
WIN_FROM <- 10.0 * 3600 * 1000  # 10:00:00.000 in ms after midnight
WIN_TO   <- 15.5 * 3600 * 1000  # 15:30:00.000

# ---- Read one day -------------------------------------------------------------
# readMat() returns the struct as a list-like array; we name its 8 fields and
# keep level 1 only (column 1 of the price/volume matrices).
read_day <- function(path) {
  m   <- readMat(path)$LOB
  lob <- setNames(lapply(seq_len(dim(m)[1]), function(i) m[[i]]), dimnames(m)[[1]])
  day <- as.IDate(sub(".*_(\\d{8})\\.mat$", "\\1", path), format = "%Y%m%d")

  book <- data.table(day    = day,
                     time   = as.numeric(lob$EventTime),
                     bid_px = lob$BuyPrice[, 1],  bid_sz = lob$BuyVolume[, 1],
                     ask_px = lob$SellPrice[, 1], ask_sz = lob$SellVolume[, 1])

  mo <- data.table(day = day,
                   time   = lob$MO[, 1],
                   bid_px = lob$MO[, 2], ask_px = lob$MO[, 3],
                   avg_px = lob$MO[, 6],
                   volume = lob$MO[, 7],
                   mo_type = as.integer(lob$MO[, 8]))  # -1 buy (lifts ask), +1 sell (hits bid)
  list(book = book, mo = mo)
}

# ---- Read all days ------------------------------------------------------------
files <- sort(list.files("data/raw", sprintf("^%s_\\d{8}\\.mat$", TICKER), full.names = TRUE))
stopifnot(length(files) > 0)
cat(sprintf("Reading %d %s files\n", length(files), TICKER))
days <- lapply(files, read_day)
book <- rbindlist(lapply(days, `[[`, "book"))
mo   <- rbindlist(lapply(days, `[[`, "mo"))

# ---- Trading window -----------------------------------------------------------
# Keep snapshots from 10:00:00.000 to 15:30:00.000 inclusive. The 10:00 snapshot
# is the starting book; each later snapshot closes one 100 ms interval.
# Market orders are kept if they fall inside an interval, i.e. (10:00, 15:30].
n_raw <- nrow(book)
book  <- book[time >= WIN_FROM & time <= WIN_TO]
mo    <- mo[time > WIN_FROM & time <= WIN_TO]

# ---- Cleaning ------------------------------------------------------------------
# No auctions/halts/hidden-order rows exist in this format. What can still go
# wrong is an empty side (size 0 or non-positive price) or a crossed/locked
# book. We count each problem before dropping it so the report can state it.
bad_empty   <- book[, bid_sz <= 0 | ask_sz <= 0 | bid_px <= 0 | ask_px <= 0]
bad_crossed <- book[, ask_px <= bid_px]
bad_offgrid <- book[, bid_px %% TICK != 0 | ask_px %% TICK != 0]
cat(sprintf("Snapshots: %d raw, %d in window; dropping %d empty-side, %d crossed/locked, %d off-tick-grid\n",
            n_raw, nrow(book), sum(bad_empty), sum(bad_crossed), sum(bad_offgrid)))
book <- book[!(bad_empty | bad_crossed | bad_offgrid)]

book[, spread_ticks := as.integer(round((ask_px - bid_px) / TICK))]
setkey(book, day, time)
setkey(mo, day, time)

# ---- Save ---------------------------------------------------------------------
dir.create("data/processed", showWarnings = FALSE)
saveRDS(book, sprintf("data/processed/lob1_%s.rds", TICKER))
saveRDS(mo,   sprintf("data/processed/mo_%s.rds",   TICKER))
cat(sprintf("Saved %d snapshots and %d market orders over %d days\n",
            nrow(book), nrow(mo), uniqueN(book$day)))
