#' Solve Prize-Collecting Steiner Forest
#'
#' Objective (maximised):
#' \deqn{\sum_{v \in V'} p_v - \sum_{e \in E'} c_e - \omega k}
#'
#' Two backends:
#' \describe{
#'   \item{\code{"builtin"}}{KMB Steiner heuristic followed by an exact
#'     bottom-up tree DP for prize pruning. Pure R + igraph, no compilation.}
#'   \item{\code{"pcsf"}}{Delegates to \code{PCSF::PCSF()}, the Goemans-Williamson
#'     heuristic of Akhmedov et al. Requires the GitHub-only PCSF package
#'     (\code{remotes::install_github("murodzhon/PCSF")}), which needs Boost.}
#' }
#'
#' Running both on identical inputs is a genuine robustness check: if two
#' different approximation schemes converge on the same core subnetwork, that is
#' worth a sentence in your methods. If they disagree sharply, your parameters
#' are in a regime where the solution is not well determined.
#'
#' @param g igraph object with edge attribute \code{cost}
#' @param prizes named numeric vector
#' @param omega per-tree cost; higher gives fewer, larger components
#' @param mu hub penalty, \code{prize_v <- prize_v - mu * degree_v}. The PCSF
#'   package authors recommend 1e-4 to 5e-2 for biological networks, chosen by
#'   inspecting the resulting Steiner/terminal ratio.
#' @param backend "builtin" or "pcsf"
#' @param cache optional environment from \code{new_sp_cache} to reuse
#'   shortest-path distances across resampling runs
#' @return list(subnetwork = igraph, info = list)
#' @export
solve_pcsf <- function(g, prizes, omega = 3, mu = 0.005,
                       backend = c("builtin", "pcsf"), cache = NULL) {
  backend <- match.arg(backend)

  deg <- igraph::degree(g)
  vn <- igraph::V(g)$name
  p <- ifelse(is.na(prizes[vn]), 0, prizes[vn])
  names(p) <- vn
  p <- p - mu * deg[vn]
  p <- p[p > 0]

  if (length(p) == 0L) {
    warning("All prizes are zero after the hub penalty.")
    return(list(subnetwork = igraph::make_empty_graph(directed = FALSE),
                info = list(n_trees = 0, objective = 0, terminals = 0, steiner = 0)))
  }

  if (backend == "pcsf") {
    if (!requireNamespace("PCSF", quietly = TRUE))
      stop("Install the PCSF package or use backend = 'builtin'.")
    sub <- PCSF::PCSF(g, p, w = omega, b = 1, mu = mu)
    return(list(subnetwork = sub, info = .summarise(sub, p, omega)))
  }

  if (is.null(cache)) cache <- new_sp_cache(g)

  comp <- igraph::components(g)
  pieces <- list()
  for (ci in seq_len(comp$no)) {
    members <- vn[comp$membership == ci]
    terms <- intersect(names(p), members)
    if (length(terms) == 0L) next
    if (length(terms) == 1L) {
      pieces[[length(pieces) + 1L]] <- igraph::make_graph(edges = character(0),
                                                          isolates = terms,
                                                          directed = FALSE)
      next
    }
    tr <- .kmb_steiner(g, terms, cache)
    tr <- .prune_tree_dp(tr, p)
    if (igraph::vcount(tr) > 0L) pieces[[length(pieces) + 1L]] <- tr
  }

  merged <- if (length(pieces) == 0L) igraph::make_empty_graph(directed = FALSE)
            else .merge_pieces(pieces)

  # forest split: an extra tree costs omega, so cut edges pricier than that
  if (igraph::ecount(merged) > 0L) {
    ec <- igraph::E(merged)$cost
    ec[is.na(ec)] <- 1
    drop <- which(ec > omega)
    if (length(drop)) merged <- igraph::delete_edges(merged, igraph::E(merged)[drop])
  }

  # Reclaim terminals that pruning disconnected. In the rooted PCSF formulation
  # any node with prize > omega pays for itself as a singleton tree, so dropping
  # it would be wrong -- it was cut from a connection, not judged worthless.
  reclaim <- setdiff(names(p)[p > omega], igraph::V(merged)$name)
  if (length(reclaim)) merged <- igraph::add_vertices(merged, length(reclaim),
                                                      name = reclaim)

  iso <- igraph::V(merged)$name[igraph::degree(merged) == 0]
  kill <- iso[is.na(p[iso]) | ifelse(is.na(p[iso]), 0, p[iso]) <= omega]
  if (length(kill)) merged <- igraph::delete_vertices(merged, kill)

  if (igraph::vcount(merged) > 0L) {
    nm <- igraph::V(merged)$name
    igraph::V(merged)$prize <- ifelse(is.na(prizes[nm]), 0, prizes[nm])
    igraph::V(merged)$terminal <- nm %in% names(p)
  }
  list(subnetwork = merged, info = .summarise(merged, p, omega))
}


# Combine the per-component pieces into one forest graph.
#
# igraph::union() CANNOT be used here: when two operands carry the same edge
# attribute (they all have `cost`), it renames them to cost_1, cost_2, ... and
# leaves `cost` NA or absent. Downstream that silently turned every edge cost
# into the NA->1 fallback, which wrecked both the forest-split (cost > omega
# test) and the objective (the primary test statistic). The pieces are drawn
# from distinct connected components of g, so they are vertex-disjoint and their
# edges never collide -- assembling one edge list preserves costs exactly.
.merge_pieces <- function(pieces) {
  verts <- unique(unlist(lapply(pieces, function(p) igraph::V(p)$name)))
  attrs <- c("cost", "is_directed", "is_stimulation", "is_inhibition")
  edfs <- lapply(pieces, function(p) {
    if (igraph::ecount(p) == 0L) return(NULL)
    el <- igraph::as_edgelist(p)
    d <- data.frame(from = el[, 1], to = el[, 2], stringsAsFactors = FALSE)
    for (a in attrs)
      if (!is.null(igraph::edge_attr(p, a))) d[[a]] <- igraph::edge_attr(p, a)
    d
  })
  edf <- do.call(rbind, edfs)
  vdf <- data.frame(name = verts, stringsAsFactors = FALSE)
  if (is.null(edf))
    edf <- data.frame(from = character(0), to = character(0), stringsAsFactors = FALSE)
  igraph::graph_from_data_frame(edf, directed = FALSE, vertices = vdf)
}


.summarise <- function(g, p, omega) {
  if (igraph::vcount(g) == 0L)
    return(list(n_trees = 0, objective = 0, terminals = 0, steiner = 0))
  nm <- igraph::V(g)$name
  k <- igraph::components(g)$no
  ec <- if (igraph::ecount(g) > 0L) igraph::E(g)$cost else numeric(0)
  ec[is.na(ec)] <- 1
  term <- nm %in% names(p)
  list(n_trees = k,
       objective = sum(ifelse(is.na(p[nm]), 0, p[nm])) - sum(ec) - omega * k,
       terminals = sum(term),
       steiner = sum(!term))
}


#' Shortest-path distance cache
#'
#' The metric closure depends only on the graph, not the prizes, so distances
#' computed during one run are reusable by every subsequent run. This is what
#' makes hundreds of noise-perturbation or permutation replicates tractable.
#'
#' @export
new_sp_cache <- function(g) {
  e <- new.env(parent = emptyenv())
  e$g <- g
  e$d <- new.env(parent = emptyenv())
  e
}

.dists <- function(cache, src) {
  if (!exists(src, envir = cache$d, inherits = FALSE)) {
    dd <- igraph::distances(cache$g, v = src, weights = igraph::E(cache$g)$cost)[1, ]
    assign(src, dd, envir = cache$d)
  }
  get(src, envir = cache$d, inherits = FALSE)
}


# KMB heuristic: metric closure -> MST -> path expansion -> MST -> leaf strip
.kmb_steiner <- function(g, terminals, cache) {
  terminals <- intersect(terminals, igraph::V(g)$name)
  if (length(terminals) < 2L)
    return(igraph::make_graph(edges = character(0), isolates = terminals,
                              directed = FALSE))

  dm <- do.call(rbind, lapply(terminals, function(t) .dists(cache, t)[terminals]))
  rownames(dm) <- terminals; colnames(dm) <- terminals
  dm[!is.finite(dm)] <- NA

  # mode = "max" rather than "undirected": as of igraph 1.6.0 the latter warns
  # on matrices it cannot verify as exactly symmetric, which floating-point
  # distance matrices routinely are not.
  cl <- igraph::graph_from_adjacency_matrix(dm, mode = "max",
                                            weighted = TRUE, diag = FALSE)
  if (igraph::ecount(cl) == 0L)
    return(igraph::make_graph(edges = character(0), isolates = terminals,
                              directed = FALSE))
  mst_cl <- igraph::mst(cl, weights = igraph::E(cl)$weight)

  el <- igraph::as_edgelist(mst_cl)
  keep <- character(0)
  for (i in seq_len(nrow(el))) {
    pth <- igraph::shortest_paths(g, from = el[i, 1], to = el[i, 2],
                                  weights = igraph::E(g)$cost,
                                  output = "epath")$epath[[1]]
    keep <- union(keep, as.character(as.numeric(pth)))
  }
  if (length(keep) == 0L)
    return(igraph::make_graph(edges = character(0), isolates = terminals,
                              directed = FALSE))

  expanded <- igraph::subgraph.edges(g, igraph::E(g)[as.numeric(keep)],
                                     delete.vertices = TRUE)
  tree <- igraph::mst(expanded, weights = igraph::E(expanded)$cost)

  repeat {
    d <- igraph::degree(tree)
    leaves <- igraph::V(tree)$name[d == 1 & !(igraph::V(tree)$name %in% terminals)]
    if (length(leaves) == 0L) break
    tree <- igraph::delete_vertices(tree, leaves)
  }
  tree
}


# Exact bottom-up pruning of a fixed tree.
#   net(v) = p(v) + sum_children max(0, net(child) - cost(v, child))
# A child subtree is kept iff it contributes positive net value. Given the tree
# topology from KMB, this pruning step is optimal, not heuristic.
.prune_tree_dp <- function(tree, p) {
  n <- igraph::vcount(tree)
  if (n == 0L) return(tree)
  nm <- igraph::V(tree)$name
  pv <- ifelse(is.na(p[nm]), 0, p[nm]); names(pv) <- nm

  root <- nm[which.max(pv)]
  # Explicit BFS rather than igraph::bfs(). The father vertex sequence returned
  # by igraph::bfs() does not survive as_ids() -- it yields the vertex names in
  # index order rather than each vertex's parent -- and indexing it directly
  # fails on the root's NA parent. A hand-rolled traversal on a tree is cheap
  # and unambiguous.
  father <- stats::setNames(rep(NA_character_, n), nm)
  seen <- stats::setNames(rep(FALSE, n), nm)
  order_nm <- character(0)
  queue <- root; seen[[root]] <- TRUE
  while (length(queue)) {
    v <- queue[1]; queue <- queue[-1]
    order_nm <- c(order_nm, v)
    for (w in igraph::as_ids(igraph::neighbors(tree, v))) {
      if (!seen[[w]]) { seen[[w]] <- TRUE; father[[w]] <- v; queue <- c(queue, w) }
    }
  }

  eid_of <- function(a, b) igraph::get.edge.ids(tree, c(a, b))
  cost_of <- function(a, b) {
    id <- eid_of(a, b)
    if (id == 0L) return(1)
    cc <- igraph::E(tree)$cost[id]
    if (is.na(cc)) 1 else cc
  }

  net <- setNames(rep(NA_real_, n), nm)
  keep <- list()
  for (v in rev(order_nm)) {                    # children before parents
    tot <- pv[[v]]
    kids <- names(father)[!is.na(father) & father == v]
    for (w in kids) {
      gain <- net[[w]] - cost_of(v, w)
      if (isTRUE(gain > 0)) { tot <- tot + gain; keep[[paste(v, w)]] <- TRUE }
      else keep[[paste(v, w)]] <- FALSE
    }
    net[[v]] <- tot
  }

  drop_e <- integer(0)
  for (v in order_nm) {
    f <- father[[v]]
    if (is.na(f)) next
    if (!isTRUE(keep[[paste(f, v)]])) {
      id <- eid_of(f, v)
      if (id != 0L) drop_e <- c(drop_e, id)
    }
  }
  if (length(drop_e)) tree <- igraph::delete_edges(tree, igraph::E(tree)[drop_e])

  comp <- igraph::components(tree)
  rid <- comp$membership[[root]]
  igraph::induced_subgraph(tree, igraph::V(tree)[comp$membership == rid])
}
