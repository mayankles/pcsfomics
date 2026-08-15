#' Edge selection frequency under prize perturbation
#'
#' PCSF is an argmax, so a single optimal solution is brittle: small changes in
#' prizes can reroute whole branches. Selection FREQUENCY across perturbed runs
#' is the quantity worth reporting. An edge recovered in 95 of 100 runs is a
#' claim; an edge appearing in one optimal solution is an anecdote.
#'
#' Note this perturbs PRIZES, whereas \code{PCSF::PCSF_rand()} perturbs EDGE
#' COSTS. The choice should follow which input you trust less. Prizes derived
#' from a ~30-sample cohort are considerably less stable than curation-based
#' confidence scores from OmniPath, so prize noise is usually the more honest
#' perturbation here. Doing both, and reporting edges stable under each, is
#' stronger still.
#'
#' @param g igraph with edge attribute \code{cost}
#' @param prizes named numeric vector
#' @param omega,mu solver parameters
#' @param n_runs number of perturbed replicates
#' @param noise_sd multiplicative Gaussian noise SD
#' @param seed RNG seed
#' @return data.frame of source, target, frequency (descending), with a
#'   \code{node_frequency} attribute
#' @export
robust_pcsf <- function(g, prizes, omega = 3, mu = 0.005,
                        n_runs = 100, noise_sd = 0.10, seed = 1) {
  set.seed(seed)
  cache <- new_sp_cache(g)
  etab <- new.env(parent = emptyenv())
  ntab <- new.env(parent = emptyenv())

  for (i in seq_len(n_runs)) {
    noisy <- pmax(0, prizes * (1 + stats::rnorm(length(prizes), 0, noise_sd)))
    names(noisy) <- names(prizes)
    sub <- solve_pcsf(g, noisy, omega = omega, mu = mu, cache = cache)$subnetwork
    if (igraph::vcount(sub) == 0L) next
    for (nmv in igraph::V(sub)$name)
      assign(nmv, (if (exists(nmv, ntab, inherits = FALSE)) get(nmv, ntab) else 0) + 1, ntab)
    if (igraph::ecount(sub) > 0L) {
      el <- igraph::as_edgelist(sub)
      keys <- apply(el, 1, function(r) paste(sort(r), collapse = "|"))
      for (k in keys)
        assign(k, (if (exists(k, etab, inherits = FALSE)) get(k, etab) else 0) + 1, etab)
    }
  }

  keys <- ls(etab)
  if (length(keys) == 0L)
    return(data.frame(source = character(0), target = character(0),
                      frequency = numeric(0), stringsAsFactors = FALSE))

  parts <- do.call(rbind, strsplit(keys, "|", fixed = TRUE))
  out <- data.frame(source = parts[, 1], target = parts[, 2],
                    frequency = vapply(keys, function(k) get(k, etab), 0) / n_runs,
                    row.names = NULL, stringsAsFactors = FALSE)
  out <- out[order(-out$frequency), , drop = FALSE]
  nf <- vapply(ls(ntab), function(k) get(k, ntab), 0) / n_runs
  attr(out, "node_frequency") <- nf
  out
}


#' Reassign prizes among genes of comparable degree
#'
#' Naive prize shuffling is anticonservative on interaction networks.
#' Differentially expressed genes skew toward well-studied, highly-connected
#' proteins, so reassigning prizes uniformly moves them onto low-degree nodes
#' and makes any degree-correlated structure look significant. Binning by degree
#' removes that confound.
#'
#' @param g igraph
#' @param prizes named numeric vector
#' @param n_bins number of degree strata
#' @return shuffled named numeric vector over all vertices
#' @export
degree_matched_shuffle <- function(g, prizes, n_bins = 10) {
  nm <- igraph::V(g)$name
  d <- igraph::degree(g)[nm]
  p <- ifelse(is.na(prizes[nm]), 0, prizes[nm]); names(p) <- nm

  br <- stats::quantile(rank(d, ties.method = "first"),
                        probs = seq(0, 1, length.out = n_bins + 1))
  bin <- cut(rank(d, ties.method = "first"), breaks = unique(br),
             include.lowest = TRUE, labels = FALSE)

  out <- p
  for (b in unique(bin[!is.na(bin)])) {
    idx <- which(bin == b)
    out[idx] <- sample(p[idx])
  }
  out
}


#' Degree-matched permutation test
#'
#' Reports empirical p-values for subnetwork compactness and Steiner fraction.
#'
#' A caveat learned from benchmarking: node count is a weak test statistic.
#' Decoy terminals inflate the observed and null subnetworks alike, so this can
#' come back non-significant even when the recovered subnetwork is enriched for
#' real signal. \code{objective} and \code{concordance enrichment within the
#' subnetwork} discriminate better; both are returned so you can compare.
#'
#' @param g igraph
#' @param prizes named numeric vector
#' @param omega,mu solver parameters
#' @param n_perm permutations
#' @param seed RNG seed
#' @export
permutation_test <- function(g, prizes, omega = 3, mu = 0.005,
                             n_perm = 200, seed = 1) {
  set.seed(seed)
  cache <- new_sp_cache(g)
  obs <- solve_pcsf(g, prizes, omega = omega, mu = mu, cache = cache)
  obs_n <- igraph::vcount(obs$subnetwork)
  obs_s <- obs$info$steiner
  obs_o <- obs$info$objective

  nn <- ns <- no <- numeric(n_perm)
  for (i in seq_len(n_perm)) {
    sh <- degree_matched_shuffle(g, prizes)
    r <- solve_pcsf(g, sh, omega = omega, mu = mu, cache = cache)
    nn[i] <- igraph::vcount(r$subnetwork)
    ns[i] <- r$info$steiner
    no[i] <- r$info$objective
  }

  list(observed_nodes = obs_n,
       observed_steiner = obs_s,
       observed_objective = obs_o,
       null_nodes_mean = mean(nn),
       null_steiner_mean = mean(ns),
       null_objective_mean = mean(no),
       p_compactness = (sum(nn <= obs_n) + 1) / (n_perm + 1),
       p_steiner     = (sum(ns >= obs_s) + 1) / (n_perm + 1),
       p_objective   = (sum(no >= obs_o) + 1) / (n_perm + 1),
       result = obs)
}
