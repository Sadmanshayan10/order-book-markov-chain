# model.R
# Purpose: Shared model definition functions used by the other scripts.
#          The chain X_n = (B_n, A_n) lives on
#            S = {(b, a) : 1 <= b, a <= K}  u  {U, D}.
#          At each event: bid +1 w.p. alpha, bid -1 w.p. beta (b = 1 -> D),
#          ask +1 w.p. alpha, ask -1 w.p. beta (a = 1 -> U), 2 alpha + 2 beta = 1.
#          A "grow" at the cap K leaves the state where it is.
#          We solve (I - Q) [H | tau] = [R | 1] for
#            h(b, a)   = P(X_tau = U | X_0 = (b, a))
#            tau(b, a) = E[tau | X_0 = (b, a)]
#          without ever forming N = (I - Q)^{-1}.
# Inputs:  none (sourced by other scripts)
# Outputs: none (defines build_P, absorb_solve, h_table, t_table)

library(Matrix)

# ---- State order --------------------------------------------------------------
# Transient states are ordered (1,1), (1,2), ..., (1,K), (2,1), ..., (K,K):
# bid-major, the same order as the K = 2 hand example.
# U and D come last, at positions K^2 + 1 and K^2 + 2.
state_index <- function(b, a, K) (b - 1) * K + a

state_names <- function(K) {
  g <- expand.grid(a = 1:K, b = 1:K)  # first column varies fastest -> bid-major
  c(sprintf("(%d,%d)", g$b, g$a), "U", "D")
}

# ---- Parameter checks ---------------------------------------------------------
# beta > 0 is needed: if queues never shrink, U and D are never reached, the
# (b, a) states are not transient and I - Q is singular.
check_params <- function(K, alpha, beta) {
  stopifnot(length(K) == 1, K >= 1, K == round(K),
            length(alpha) == 1, length(beta) == 1,
            alpha >= 0, beta > 0,
            abs(2 * alpha + 2 * beta - 1) < 1e-8)
}

# ---- One-step transition matrix P ---------------------------------------------
# Returns a sparse (K^2 + 2) x (K^2 + 2) matrix with named rows/columns.
# Each transient row has (at most) four non-zero entries, one per event type.
# Two moves that land on the same state (e.g. both "grow" moves at (K, K))
# are added together: sparseMatrix() sums duplicate (i, j) pairs.
build_P <- function(K, alpha, beta) {
  check_params(K, alpha, beta)
  n_tr <- K^2
  iU <- n_tr + 1
  iD <- n_tr + 2

  g <- expand.grid(a = 1:K, b = 1:K)
  b <- g$b
  a <- g$a
  from <- state_index(b, a, K)  # = 1:n_tr, written out for readability

  to_bid_up   <- state_index(pmin(b + 1, K), a, K)            # capped: stays at K
  to_bid_down <- ifelse(b == 1, iD, state_index(b - 1, a, K))  # bid empties -> D
  to_ask_up   <- state_index(b, pmin(a + 1, K), K)            # capped: stays at K
  to_ask_down <- ifelse(a == 1, iU, state_index(b, a - 1, K))  # ask empties -> U

  i <- c(rep(from, 4), iU, iD)
  j <- c(to_bid_up, to_bid_down, to_ask_up, to_ask_down, iU, iD)
  x <- c(rep(c(alpha, beta, alpha, beta), each = n_tr), 1, 1)  # U, D absorbing

  nm <- state_names(K)
  sparseMatrix(i = i, j = j, x = x, dims = c(n_tr + 2, n_tr + 2),
               dimnames = list(nm, nm))
}

# ---- Generic absorption solver ------------------------------------------------
# Works for ANY absorbing chain given as a (sparse or dense) matrix P with
# row/column names, plus the names of its absorbing states. This is what lets
# tests.R feed it a 1D Gambler's Ruin chain.
#   Q = transient -> transient block, R = transient -> absorbing block.
#   First-step analysis:  H = Q H + R  and  tau = 1 + Q tau,  i.e.
#     (I - Q) H = R,   (I - Q) tau = 1.
# Both systems share the matrix I - Q, so we solve them in one call with the
# right-hand side [R | 1]. solve() on a sparse matrix does a sparse LU
# factorisation; nothing is inverted.
absorb_solve <- function(P, absorbing) {
  ab <- match(absorbing, rownames(P))
  stopifnot(!anyNA(ab))
  tr <- setdiff(seq_len(nrow(P)), ab)

  Q <- P[tr, tr, drop = FALSE]
  R <- P[tr, ab, drop = FALSE]
  I_minus_Q <- Diagonal(length(tr)) - Q

  rhs <- cbind(as.matrix(R), 1)
  sol <- as.matrix(solve(I_minus_Q, rhs))
  dimnames(sol) <- list(rownames(P)[tr], c(absorbing, "tau"))

  list(H   = sol[, absorbing, drop = FALSE],  # absorption probabilities
       tau = sol[, "tau"])                     # expected events to absorption
}

# ---- K x K tables --------------------------------------------------------------
# solve_chain() does the work once; h_table() and t_table() pick out one piece.
# Results come back in bid-major order, so matrix(..., byrow = TRUE) puts
# state (b, a) at entry [b, a].
solve_chain <- function(K, alpha, beta) {
  sol <- absorb_solve(build_P(K, alpha, beta), absorbing = c("U", "D"))
  to_grid <- function(v) {
    matrix(v, nrow = K, ncol = K, byrow = TRUE,
           dimnames = list(b = 1:K, a = 1:K))
  }
  list(h = to_grid(sol$H[, "U"]), tau = to_grid(sol$tau))
}

# h_table(K, alpha, beta)[b, a] = h(b, a), the probability the next move is up.
h_table <- function(K, alpha, beta) solve_chain(K, alpha, beta)$h

# t_table(K, alpha, beta)[b, a] = tau(b, a), expected number of EVENTS (not
# seconds) until the price moves.
t_table <- function(K, alpha, beta) solve_chain(K, alpha, beta)$tau
