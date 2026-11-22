# tests.R
# Purpose: Tests for the functions in R/model.R. Run from the project root:
#            Rscript R/tests.R
#          Fixed points:
#            1. P is a valid transition matrix; K = 2 Q and R match the hand-derived matrices
#            2. K = 2, alpha = 1/6, beta = 1/3 hand values of h and tau
#            3. h(b, b) = 1/2 (and h(b, a) + h(a, b) = 1) for any K, alpha
#            4. the solver reproduces 1D Gambler's Ruin (1 - r^i) / (1 - r^N)
#            5. a K = 30 reference table (3 decimals)
#            6. brute-force Monte Carlo agrees within Monte Carlo error
#            7. h_table() is fast enough at a realistic K
# Inputs:  R/model.R
# Outputs: none (prints PASS/FAIL per test; errors if any test fails)

source("R/model.R")

# ---- Tiny test harness --------------------------------------------------------
n_pass <- 0
n_fail <- 0
check <- function(name, ok) {
  ok <- isTRUE(all(ok))
  cat(if (ok) "PASS  " else "FAIL  ", name, "\n", sep = "")
  if (ok) n_pass <<- n_pass + 1 else n_fail <<- n_fail + 1
}
# Largest absolute difference, ignoring names/dimnames.
max_diff <- function(x, y) max(abs(as.numeric(x) - as.numeric(y)))

# Read a K x K table back out in bid-major state order
# (1,1), (1,2), (2,1), (2,2), ... (bid-major).
as_state_vector <- function(M) as.numeric(t(M))

# ---- 1. Structure of P --------------------------------------------------------
for (K in c(1, 2, 5, 30)) {
  P <- build_P(K, 1/6, 1/3)
  check(sprintf("P rows sum to 1 (K = %d)", K), max(abs(rowSums(P) - 1)) < 1e-12)
}
P2 <- build_P(2, 1/6, 1/3)
Q_hand <- (1/6) * rbind(c(0, 1, 1, 0),
                        c(2, 1, 0, 1),
                        c(2, 0, 1, 1),
                        c(0, 2, 2, 2))
R_hand <- (1/3) * rbind(c(1, 1),
                        c(0, 1),
                        c(1, 0),
                        c(0, 0))
check("K = 2: Q matches the hand-derived matrix", max_diff(P2[1:4, 1:4], Q_hand) < 1e-12)
check("K = 2: R matches the hand-derived matrix", max_diff(P2[1:4, 5:6], R_hand) < 1e-12)

# ---- 2. K = 2 hand values -----------------------------------------------------
h2 <- h_table(2, 1/6, 1/3)
t2 <- t_table(2, 1/6, 1/3)
check("K = 2: h = (1/2, 3/10, 7/10, 1/2)",
      max_diff(as_state_vector(h2), c(1/2, 3/10, 7/10, 1/2)) < 1e-12)
check("K = 2: tau = (1.95, 2.85, 2.85, 4.35)",
      max_diff(as_state_vector(t2), c(1.95, 2.85, 2.85, 4.35)) < 1e-12)
check("K = 2: h[2, 1] = 0.7 is entry [b = 2, a = 1]", abs(h2[2, 1] - 0.7) < 1e-12)

# ---- 3. Symmetry: h(b, b) = 1/2 and h(b, a) + h(a, b) = 1 ----------------------
# Same alpha, beta on both sides, so swapping bid and ask swaps U and D.
# Tolerance 1e-6, not 1e-10: when alpha > beta the queues drift UP to the cap,
# tau reaches ~6e9 events at K = 30 and I - Q has condition number ~3e12, so
# floating-point error in h near the cap is ~3e-7. That is far below anything
# that matters for prediction, but it would fail a 1e-10 check.
for (K in c(5, 30)) {
  for (alpha in c(1/6, 1/4, 1/3)) {
    H <- h_table(K, alpha, 1/2 - alpha)
    check(sprintf("h(b, b) = 1/2              (K = %d, alpha = %.3f)", K, alpha),
          max(abs(diag(H) - 1/2)) < 1e-6)
    check(sprintf("h(b, a) + h(a, b) = 1      (K = %d, alpha = %.3f)", K, alpha),
          max(abs(H + t(H) - 1)) < 1e-6)
  }
}

# ---- 4. 1D Gambler's Ruin -----------------------------------------------------
# States 0, 1, ..., N; 0 and N absorbing; win $1 w.p. p, lose $1 w.p. 1 - p.
# Row/column for state s is s + 1. We hand this matrix to the SAME generic
# absorb_solve() that h_table() uses.
gr_P <- function(N, p) {
  nm <- as.character(0:N)
  rows <- 2:N  # transient states 1, ..., N - 1
  sparseMatrix(i = c(rows, rows, 1, N + 1),
               j = c(rows + 1, rows - 1, 1, N + 1),
               x = c(rep(p, N - 1), rep(1 - p, N - 1), 1, 1),
               dims = c(N + 1, N + 1), dimnames = list(nm, nm))
}
gr_formula <- function(i, N, p) {
  r <- (1 - p) / p
  if (abs(r - 1) < 1e-12) i / N else (1 - r^i) / (1 - r^N)
}
N <- 10
for (p in c(0.4, 0.5, 0.6)) {
  g_solver <- absorb_solve(gr_P(N, p), absorbing = c(as.character(N), "0"))$H[, 1]
  check(sprintf("Gambler's Ruin, N = %d, p = %.1f, all i", N, p),
        max_diff(g_solver, gr_formula(1:(N - 1), N, p)) < 1e-10)
}
g3 <- absorb_solve(gr_P(5, 0.4), absorbing = c("5", "0"))$H["3", 1]
check("Gambler's Ruin, textbook example i = 3, N = 5, p = 0.4 gives 0.360",
      abs(g3 - 0.360) < 5e-4)

# ---- 5. Reference table (K = 30, rounded to 3 decimals) -----------------------
tab76 <- data.frame(b = c(2, 3, 4), a = c(1, 1, 2),
                    a14 = c(0.698, 0.791, 0.703),   # alpha = beta = 1/4
                    a16 = c(0.722, 0.831, 0.755),   # alpha = 1/6, beta = 1/3
                    a13 = c(0.653, 0.708, 0.598))   # alpha = 1/3, beta = 1/6
for (col in c("a14", "a16", "a13")) {
  alpha <- c(a14 = 1/4, a16 = 1/6, a13 = 1/3)[[col]]
  H <- h_table(30, alpha, 1/2 - alpha)
  got <- H[cbind(tab76$b, tab76$a)]
  check(sprintf("Reference table, K = 30, alpha = %.3f", alpha),
        max(abs(got - tab76[[col]])) <= 5e-4 + 1e-12)
}

# ---- 6. Brute-force Monte Carlo -----------------------------------------------
# Simulates the event rules directly (NOT by sampling rows of build_P), so a
# mistake in build_P would show up as a disagreement here.
# All paths are advanced together; finished paths are dropped from `alive`.
simulate_hit <- function(b0, a0, K, alpha, beta, n_paths) {
  b <- rep(b0, n_paths)
  a <- rep(a0, n_paths)
  steps <- integer(n_paths)
  up    <- logical(n_paths)
  alive <- rep(TRUE, n_paths)
  while (any(alive)) {
    idx  <- which(alive)
    move <- sample.int(4, length(idx), replace = TRUE,
                       prob = c(alpha, beta, alpha, beta))  # BA, BR, AA, AR
    bi <- pmin(b[idx] + (move == 1) - (move == 2), K)       # grow at K: stay
    ai <- pmin(a[idx] + (move == 3) - (move == 4), K)
    b[idx] <- bi
    a[idx] <- ai
    steps[idx] <- steps[idx] + 1L
    up[idx[ai == 0]] <- TRUE                                # ask emptied -> U
    alive[idx[ai == 0 | bi == 0]] <- FALSE
  }
  list(h = mean(up),    se_h = sqrt(mean(up) * (1 - mean(up)) / n_paths),
       tau = mean(steps), se_tau = sd(steps) / sqrt(n_paths))
}

set.seed(348)
K_mc <- 5
for (sc in list(c(b = 2, a = 1, alpha = 1/6),
                c(b = 4, a = 2, alpha = 1/6),
                c(b = 3, a = 1, alpha = 1/4))) {
  alpha <- sc[["alpha"]]
  beta  <- 1/2 - alpha
  sol   <- solve_chain(K_mc, alpha, beta)
  mc    <- simulate_hit(sc[["b"]], sc[["a"]], K_mc, alpha, beta, n_paths = 20000)
  h_ex  <- sol$h[sc[["b"]], sc[["a"]]]
  t_ex  <- sol$tau[sc[["b"]], sc[["a"]]]
  lab   <- sprintf("(%d,%d), K = %d, alpha = %.3f", sc[["b"]], sc[["a"]], K_mc, alpha)
  cat(sprintf("      MC %s: h %.4f vs exact %.4f (se %.4f); tau %.3f vs exact %.3f (se %.3f)\n",
              lab, mc$h, h_ex, mc$se_h, mc$tau, t_ex, mc$se_tau))
  check(sprintf("Monte Carlo h   within 4 se, %s", lab), abs(mc$h - h_ex) < 4 * mc$se_h)
  check(sprintf("Monte Carlo tau within 4 se, %s", lab), abs(mc$tau - t_ex) < 4 * mc$se_tau)
}

# ---- 7. Speed (under 5 s at a realistic K) ------------------------------------
# K = 100 means 10,000 transient states (the baseline is K = 155).
el <- system.time(H100 <- h_table(100, 1/6, 1/3))[["elapsed"]]
cat(sprintf("      h_table(K = 100) took %.2f s\n", el))
check("h_table(K = 100) runs in under 5 seconds", el < 5)

# ---- Summary ------------------------------------------------------------------
cat(sprintf("\n%d passed, %d failed\n", n_pass, n_fail))
if (n_fail > 0) stop("some tests failed")
