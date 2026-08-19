## propagation.R -- signed propagation over a directed regulatory network, and
## scoring of the predicted state against observed data.
##
## The premise: a subnetwork is a causal hypothesis, not just a connectivity
## claim. Given seed nodes with known states, the signed directed edges predict
## the state of everything downstream. How well that prediction matches observed
## expression is a "goodness" score for the subnetwork -- a model-selection
## criterion grounded in mechanism rather than in size or objective value.
##
## Validated against planted ground truth in test_propagation.R (sign algebra,
## planted-signal recovery, and a calibrated negative control). Requires the
## Matrix package.

#' Signed directed influence matrix from an edge table
#'
#' \code{W[i, j]} is the signed, confidence-weighted influence of source i on
#' target j. Note this is DIRECTED and must be built from the raw interaction
#' table -- \code{load_omnipath()} returns an undirected graph with reciprocal
#' edges collapsed, which is right for PCSF and wrong here, since A->B and B->A
#' are different edges that may carry different signs.
#'
#' Duplicate (i,j) pairs are summed, so contradictory activating/inhibiting
#' records for the same ordered pair partially cancel.
#'
#' @param source,target character vectors of node names
#' @param sign +1 activation, -1 inhibition, 0 to drop the edge
#' @param weight edge confidence, recycled if length 1
#' @param nodes optional node universe; defaults to observed nodes
#' @return sparse Matrix with dimnames
#' @export
build_signed_matrix <- function(source, target, sign, weight = 1, nodes = NULL) {
  if (!requireNamespace("Matrix", quietly = TRUE))
    stop("build_signed_matrix() needs the Matrix package.")
  keep <- !is.na(source) & !is.na(target) & sign != 0 & source != target
  source <- source[keep]; target <- target[keep]; sign <- sign[keep]
  weight <- if (length(weight) == 1L) rep(weight, length(source)) else weight[keep]
  if (is.null(nodes)) nodes <- sort(unique(c(source, target)))
  idx <- stats::setNames(seq_along(nodes), nodes)
  Matrix::sparseMatrix(i = idx[source], j = idx[target], x = sign * weight,
                       dims = c(length(nodes), length(nodes)),
                       dimnames = list(nodes, nodes))
}


#' Damped signed propagation
#'
#' \deqn{\hat{s} = \sum_{k=1}^{K} \lambda^k (W_{norm}^T)^k s_0}
#'
#' Column normalisation makes a target's incoming influence a weighted mean over
#' its regulators rather than a sum, so a gene with fifty regulators does not
#' automatically dominate. Finite K with \code{lambda < 1} keeps feedback loops
#' from ringing, which matters because regulatory networks are dense in cycles.
#'
#' The seed vector itself is not accumulated -- only its downstream consequences
#' -- since seeds are excluded from scoring anyway.
#'
#' @param W signed influence matrix from \code{build_signed_matrix}
#' @param seeds named numeric vector of seed states
#' @param nodeset optional node subset to restrict propagation to
#' @param lambda damping per hop
#' @param k_hops number of hops
#' @param normalise column-normalise incoming influence
#' @return named numeric vector of predicted states
#' @export
propagate <- function(W, seeds, nodeset = NULL, lambda = 0.5, k_hops = 4,
                      normalise = TRUE) {
  nodes <- rownames(W)
  if (!is.null(nodeset)) {
    nodes <- intersect(nodes, nodeset)
    W <- W[nodes, nodes, drop = FALSE]
  }
  if (normalise) {
    cin <- Matrix::colSums(abs(W)); cin[cin == 0] <- 1
    W <- W %*% Matrix::Diagonal(x = 1 / cin)
  }
  s <- stats::setNames(rep(0, length(nodes)), nodes)
  common <- intersect(names(seeds), nodes)
  s[common] <- seeds[common]
  acc <- stats::setNames(rep(0, length(nodes)), nodes)
  for (k in seq_len(k_hops)) {
    s <- as.numeric(lambda * Matrix::crossprod(W, s))  # target_j <- sum_i W[i,j] s_i
    acc <- acc + s
  }
  acc
}


#' Score predicted states against observed statistics
#'
#' Sign accuracy is reported twice, and the distinction matters. Propagation
#' reaches many genes with a predicted magnitude near zero, where the predicted
#' SIGN is arbitrary; grading those dilutes the metric toward 0.5. On planted
#' ground truth with almost no noise, accuracy across all reached genes read
#' 0.72 while the confident tail read 0.98. \code{sign_acc_conf} restricts to
#' the top \code{1 - conf_quantile} fraction by |prediction| and is the number
#' to trust.
#'
#' @param pred named numeric vector of predicted states
#' @param observed named numeric vector of observed statistics (e.g. DE t)
#' @param exclude nodes to drop from evaluation, normally the seeds
#' @param sig optional named logical, TRUE where the observed change is
#'   significant; sign is only graded there, since a null gene's sign is noise
#' @param min_abs magnitude above which a node counts as reached
#' @param conf_quantile quantile of |prediction| defining the confident tail
#' @return list of correlation and accuracy statistics
#' @export
score_prediction <- function(pred, observed, exclude = character(0),
                             sig = NULL, min_abs = 1e-9, conf_quantile = 0.75) {
  eval_genes <- setdiff(names(observed), exclude)
  p <- stats::setNames(rep(0, length(eval_genes)), eval_genes)
  common <- intersect(names(pred), eval_genes)
  p[common] <- pred[common]
  obs <- observed[eval_genes]

  reached <- names(p)[abs(p) > min_abs]
  rho_g <- suppressWarnings(stats::cor(p, obs, method = "spearman"))
  rho_l <- if (length(reached) >= 10)
    suppressWarnings(stats::cor(p[reached], obs[reached], method = "spearman")) else NA_real_

  hits <- if (is.null(sig)) reached else reached[sig[reached] %in% TRUE]
  acc_all <- if (length(hits) >= 10) mean(sign(p[hits]) == sign(obs[hits])) else NA_real_

  conf <- character(0)
  if (length(hits) >= 10) {
    thr <- stats::quantile(abs(p[hits]), conf_quantile)
    conf <- hits[abs(p[hits]) >= thr]
  }
  acc_conf <- if (length(conf) >= 10) mean(sign(p[conf]) == sign(obs[conf])) else NA_real_
  rho_conf <- if (length(conf) >= 10)
    suppressWarnings(stats::cor(p[conf], obs[conf], method = "spearman")) else NA_real_

  list(rho_global = rho_g, rho_local = rho_l, rho_conf = rho_conf,
       sign_acc = acc_all, sign_acc_conf = acc_conf,
       n_reached = length(reached), n_sign_scored = length(hits),
       n_conf = length(conf))
}


#' Degree bins for seed relocation
#'
#' Naive seed relocation is anticonservative: seeds sit on well-studied,
#' high-degree genes, and moving them onto low-degree nodes guarantees weaker
#' propagation in the null. Binning by total degree removes that confound. Same
#' reasoning as \code{degree_matched_shuffle()} in robustness.R.
#'
#' @export
make_degree_binner <- function(W, n_bins = 10) {
  A <- abs(W) > 0
  deg <- Matrix::colSums(A) + Matrix::rowSums(A)
  bin <- cut(rank(deg, ties.method = "first"), breaks = n_bins, labels = FALSE)
  names(bin) <- rownames(W)
  bin
}

#' Relocate seeds within degree bins, states riding along
#' @export
shuffle_seeds <- function(seeds, bin) {
  nodes <- names(bin)
  bins <- bin[names(seeds)]
  nm_new <- character(length(seeds))
  for (b in unique(bins[!is.na(bins)])) {
    at <- which(bins == b)
    nm_new[at] <- sample(nodes[bin == b], length(at))
  }
  ok <- nm_new != ""
  stats::setNames(as.numeric(seeds)[ok], nm_new[ok])
}
