# Decisions

| Date | Decision | Reason |
|------|----------|--------|
| 2026-11-13 | Data source: Cartea-Jaimungal-Penalva Nov 2014 NASDAQ files (.mat, 100 ms level-1 snapshots), MSFT first; INTC not analysed | LOBSTER samples are no longer free; CJP data is free and comes with a textbook |
| 2026-11-13 | Events = net size change per 100 ms interval: BA/BR/AA/AR if best prices unchanged, PX if a best price changed; one interval can give a bid and an ask event | Data are snapshots, not messages; several real events can fall in one interval, so counts are lower bounds (limitation) |
| 2026-11-13 | `MO_LAG_MS = 100`: assign each market order to the interval one 100 ms bin after its timestamp | With no offset only 20% of sell MOs coincide with the bid shrinking or a price move; with one bin 92.5% do (buy MOs/ask: 21% vs 93%). The snapshot stamped t shows the book as of t - 100 |
| 2026-11-13 | events_MSFT.csv columns: `day` added; `type` = number of MOs on that side in the interval; `direction` = +1 bid side, -1 ask side, or sign of the mid move for PX rows | No LOBSTER message types or directions exist in snapshot data |
| 2026-11-13 | L = 100 shares, K = 155 (baseline) | L = one round lot = median nonzero net change (109); K = 99th pct queue (~15,400 shares) / L. Rule fixed before seeing results; L = 500, 1000 as sensitivity in 04_validate.R |
| 2026-11-13 | Estimate alpha, beta with and without counting one removal on the emptied side in each PX interval | Depletions are recorded as PX, not BR/AR, which biases beta-hat down |
| 2026-11-13 | Treat t(b, a) as qualitative only; K-sensitivity of h at K = 100, 155, 250 in 04_validate.R | alpha-hat > beta-hat, so queues drift to the cap and t at large queues is driven by K (orders of magnitude between K = 50 and 155); h is stable |
| 2026-11-13 | Do not commit data/processed/events_MSFT.csv (52 MB) or lob1_MSFT.rds | Over the 20 MB limit we set for committed files; rebuilt from raw data by 01_parse.R + 02_events.R in ~30 s |
| 2026-11-13 | Train = first 11 days (Nov 3-17), test = last 7 (Nov 18-26); n = number of BA/BR/AA/AR events; CI for alpha uses Var = alpha(1/2 - alpha)/n | Split by day per plan; counting events makes alpha + beta = 1/2 |
| 2026-11-13 | samples_MSFT.csv adds b_sh, a_sh (shares) and wait_ms; not committed (22 MB) | Shares allow the L = 500, 1000 sensitivity without rebuilding; size over the 20 MB rule |
| 2026-11-13 | Chain baseline uses the depletion-corrected alpha-hat; uncorrected reported alongside. Logistic regression on I (Gould & Bonart 2016) added as the fitted benchmark | Corrected estimate is less biased; logistic is the project's main comparison |
| 2026-11-13 | L-sensitivity keeps the baseline cap in shares (K x L = 15,500), so L = 500 -> K = 31, L = 1000 -> K = 16; alpha, beta unchanged across L | Same 99th-pct rule as the baseline; alpha, beta count events, not shares |
| 2026-11-13 | Brier differences get 95% intervals from a bootstrap over whole test days (B = 2000), plus per-day scores | Samples within a day are dependent; row-level resampling would be far too narrow. Only 7 test days, so intervals are coarse |
| 2026-11-13 | Assumption tests run on train days; depletions counted as removals; A1 uses only pairs of consecutive intervals that each contain exactly one event | Matches the estimation sample; when both sides change in one 100 ms interval their order is unknown, and any fixed order biased the table (bid share of "next" events was 64%) |
