## test_propagation.R -- ground-truth validation of R/propagation.R.
##
## A null result on real data is only interpretable if the machinery provably
## recovers signal that IS there. Three checks:
##   1. sign algebra on a hand-built chain (inhibition must flip)
##   2. planted signal on a random graph (must recover, with a significant p)
##   3. negative control (pure noise must score ~0.5 and be non-significant)
##
##   Rscript test_propagation.R
source("R/propagation.R")
msg <- function(...) message(sprintf(...))
set.seed(42)

## ---------------------------------------------------------------------------
## 1. Sign algebra: A --| B --| C  and  A --> D --| E
## ---------------------------------------------------------------------------
W <- build_signed_matrix(source = c("A", "B", "A", "D"),
                         target = c("B", "C", "D", "E"),
                         sign   = c(-1, -1,  1, -1))
p <- propagate(W, c(A = 1), lambda = 0.5, k_hops = 4)
msg("chain: B=%.4f C=%.4f D=%.4f E=%.4f", p["B"], p["C"], p["D"], p["E"])
stopifnot(p["B"] < 0,   # activation into an inhibiting edge
          p["C"] > 0,   # double inhibition restores sign
          p["D"] > 0,   # activation preserves
          p["E"] < 0)   # then inhibited
msg("PASS 1: sign algebra")

## ---------------------------------------------------------------------------
## 2. Planted signal on a random signed graph
## ---------------------------------------------------------------------------
N <- 600; M <- 3000
nodes <- sprintf("G%03d", seq_len(N))
src <- sample(nodes, M, replace = TRUE)
tgt <- sample(nodes, M, replace = TRUE)
sgn <- sample(c(1, -1), M, replace = TRUE, prob = c(0.8, 0.2))  # OmniPath-like
Wr <- build_signed_matrix(src, tgt, sgn, weight = runif(M, 0.3, 1), nodes = nodes)

seeds <- stats::setNames(sample(c(-1, 1), 25, replace = TRUE), sample(nodes, 25))
truth <- propagate(Wr, seeds, lambda = 0.5, k_hops = 4)
bin <- make_degree_binner(Wr)

run_null <- function(obs, n = 300) {
  vapply(seq_len(n), function(i)
    score_prediction(propagate(Wr, shuffle_seeds(seeds, bin), lambda = 0.5, k_hops = 4),
                     obs, exclude = names(seeds))$rho_local, 0)
}

for (nsd in c(0.02, 0.1, 0.5)) {
  obs <- stats::setNames(as.numeric(truth + stats::rnorm(length(truth), 0, nsd)),
                         names(truth))
  sc <- score_prediction(propagate(Wr, seeds, lambda = 0.5, k_hops = 4),
                         obs, exclude = names(seeds))
  nl <- run_null(obs)
  pe <- (sum(nl >= sc$rho_local) + 1) / (length(nl) + 1)
  msg("planted (noise_sd=%.2f): rho_local=%.3f rho_conf=%.3f sign_all=%.3f sign_conf=%.3f p=%.4f",
      nsd, sc$rho_local, sc$rho_conf, sc$sign_acc, sc$sign_acc_conf, pe)
  ## rho_conf / sign_acc_conf are the sensitive statistics; rho_local dilutes
  ## quickly once noise swamps the many near-zero predictions.
  if (nsd <= 0.02) stopifnot(sc$rho_conf > 0.90, sc$sign_acc_conf > 0.90, pe < 0.01)
  if (nsd <= 0.10) stopifnot(sc$rho_conf > 0.50, sc$sign_acc_conf > 0.70, pe < 0.01)
}
msg("PASS 2: planted signal recovered and significant vs degree-matched null")

## ---------------------------------------------------------------------------
## 3. Negative control: observed unrelated to the graph
## ---------------------------------------------------------------------------
accs <- ps <- numeric(30)
for (i in seq_len(30)) {
  obs <- stats::setNames(stats::rnorm(length(truth)), names(truth))
  sc <- score_prediction(propagate(Wr, seeds, lambda = 0.5, k_hops = 4),
                         obs, exclude = names(seeds))
  nl <- run_null(obs, n = 100)
  accs[i] <- sc$sign_acc_conf
  ps[i] <- (sum(nl >= sc$rho_local) + 1) / (length(nl) + 1)
}
msg("negative control: mean sign_acc_conf=%.3f (expect ~0.50), frac p<0.05=%.2f (expect ~0.05)",
    mean(accs), mean(ps < 0.05))
stopifnot(abs(mean(accs) - 0.5) < 0.05, mean(ps < 0.05) < 0.15)
msg("PASS 3: negative control uninformative, null calibrated")

msg("ALL PROPAGATION CHECKS PASSED")
