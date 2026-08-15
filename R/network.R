#' Load an OmniPath interaction table into a weighted igraph object
#'
#' Edge cost is \code{1 - rank_normalise(confidence)}, floored at
#' \code{min_cost}. Well-curated interactions are cheap to traverse, so the
#' solver crosses them freely; weakly-supported ones are expensive and only get
#' used when the prize payoff justifies it.
#'
#' Direction and sign are carried through as edge attributes for post-hoc
#' interpretation. The PCSF objective itself is undirected -- worth remembering
#' when you write the methods section, since it means a recovered edge asserts
#' connection, not causal flow.
#'
#' Get the table with:
#' \preformatted{
#'   library(OmnipathR)
#'   ia <- import_omnipath_interactions()
#'   write.table(ia, "omnipath.tsv", sep = "\t", row.names = FALSE, quote = FALSE)
#' }
#'
#' @param path tab-separated OmniPath export
#' @param source_col,target_col,confidence_col column names
#' @param min_cost,max_cost cost range
#' @return igraph object with edge attribute \code{cost}
#' @export
load_omnipath <- function(path,
                          source_col = "source_genesymbol",
                          target_col = "target_genesymbol",
                          confidence_col = "curation_effort",
                          min_cost = 0.05,
                          max_cost = 1.0) {
  df <- utils::read.delim(path, stringsAsFactors = FALSE)
  df <- df[!is.na(df[[source_col]]) & !is.na(df[[target_col]]), , drop = FALSE]
  df <- df[df[[source_col]] != "" & df[[target_col]] != "", , drop = FALSE]
  df <- df[df[[source_col]] != df[[target_col]], , drop = FALSE]
  
  if (!is.null(df[[confidence_col]])) {
    conf <- rank_normalise(ifelse(is.na(df[[confidence_col]]), 0,
                                  df[[confidence_col]]))
  } else {
    warning(sprintf("No '%s' column; using uniform edge costs.", confidence_col))
    conf <- rep(0.5, nrow(df))
  }
  df$cost <- pmin(pmax(max_cost * (1 - conf), min_cost), max_cost)
  
  # collapse reciprocal duplicates, keeping the best-supported (cheapest) edge
  key <- apply(cbind(df[[source_col]], df[[target_col]]), 1,
               function(r) paste(sort(r), collapse = "|"))
  df <- df[order(df$cost), , drop = FALSE]
  key <- key[order(df$cost)]
  df <- df[!duplicated(key), , drop = FALSE]
  
  edges <- data.frame(from = df[[source_col]], to = df[[target_col]],
                      cost = df$cost, stringsAsFactors = FALSE)
  for (a in c("is_directed", "is_stimulation", "is_inhibition")) {
    if (!is.null(df[[a]])) edges[[a]] <- df[[a]]
  }
  igraph::graph_from_data_frame(edges, directed = FALSE)
}


#' Suggested starting value for omega
#'
#' omega is a per-tree admission charge, so it only does useful work when it
#' sits inside the prize distribution. Far below it, every terminal is retained
#' and the forest is trivially large; far above it, nothing is admitted.
#'
#' @export
suggest_omega <- function(prizes, quantile = 0.90) {
  v <- prizes[prizes > 0]
  if (length(v) == 0L) return(1)
  as.numeric(stats::quantile(v, quantile))
}


#' Component size distribution of a solution
#'
#' Singleton components restate what the prize vector already said. Multi-node
#' components are where the network contributed information, so report those.
#'
#' @export
component_summary <- function(g) {
  if (igraph::vcount(g) == 0L) {
    return(list(n_components = 0, n_singletons = 0, n_multinode = 0,
                largest = 0, sizes = integer(0)))
  }
  cs <- igraph::components(g)$csize
  list(n_components = length(cs),
       n_singletons = sum(cs == 1),
       n_multinode  = sum(cs > 1),
       largest      = max(cs),
       sizes        = sort(cs[cs > 1], decreasing = TRUE))
}


#' Valid `fields` values for the OmniPath interactions query
#'
#' The server rejects unknown fields outright rather than ignoring them, and it
#' answers with an HTML error page reading "Something is not entirely good."
#' rather than an HTTP error status, so a bad field looks like a malformed table
#' rather than a failed request. Validating client-side turns that into a clear
#' message.
#'
#' Live list: https://omnipathdb.org/queries/interactions
#'
#' Note that `n_references` and `n_resources` appear in OmnipathR's output but
#' are NOT server fields -- OmnipathR computes them from the `references` and
#' `sources` strings after download. \code{fetch_omnipath()} does the same.
#'
#' @export
OMNIPATH_FIELDS <- c(
  "curation_effort", "databases", "datasets", "dorothea_level",
  "dorothea_methods", "entity_type", "evidences", "extra_attrs",
  "ncbi_tax_id", "ncbi_tax_id_source", "ncbi_tax_id_target",
  "organism", "references", "resources", "sources", "type")

#' Valid `datasets` values for the OmniPath interactions query
#' @export
OMNIPATH_DATASETS <- c(
  "collectri", "dorothea", "kinaseextra", "ligrecextra", "lncrna_mrna",
  "mirnatarget", "omnipath", "pathwayextra", "small_molecule", "tf_mirna",
  "tf_target")


#' Download OmniPath interactions directly from the web service
#'
#' Bypasses \code{OmnipathR::import_omnipath_interactions()}, which is prone to
#' a failure unrelated to the interaction data itself. Before running a query,
#' OmnipathR validates the \code{organisms} parameter by translating names to
#' NCBI taxon IDs, and that lookup scrapes the Ensembl species HTML page with
#' rvest. When the page shape changes, \code{html_element()} returns an
#' \code{xml_missing} object and \code{html_table()} raises
#' "no applicable method for \'html_table\' applied to an object of class
#' xml_missing". The OmniPath web service is fine; only the organism name
#' lookup is broken. Requesting the taxon ID numerically over plain HTTP avoids
#' the code path entirely.
#'
#' @param datasets character vector, see \code{OMNIPATH_DATASETS}. "omnipath" is
#'   the curated causal signalling network. Add "collectri" or "dorothea" for
#'   TF-target edges, or "pathwayextra" / "kinaseextra" / "ligrecextra" for
#'   larger, less curated sets.
#' @param organism NCBI taxon ID as an integer. The service serves 9606 (human),
#'   10090 (mouse) and 10116 (rat). Pass the number, not a name.
#' @param fields extra columns, see \code{OMNIPATH_FIELDS}
#' @param license "academic" or "commercial"
#' @param destfile optional path to write the raw TSV
#' @return data.frame of interactions, with n_references and n_resources
#'   computed client-side
#' @export
fetch_omnipath <- function(datasets = "omnipath",
                           organism = 9606L,
                           fields = c("curation_effort", "sources", "references"),
                           license = "academic",
                           destfile = NULL) {
  stopifnot(is.numeric(organism))
  
  bad_f <- setdiff(fields, OMNIPATH_FIELDS)
  if (length(bad_f))
    stop("Not valid OmniPath fields: ", paste(bad_f, collapse = ", "),
         "\n  Valid: ", paste(OMNIPATH_FIELDS, collapse = ", "),
         "\n  (n_references and n_resources are computed here, not requested.)")
  
  bad_d <- setdiff(datasets, OMNIPATH_DATASETS)
  if (length(bad_d))
    stop("Not valid OmniPath datasets: ", paste(bad_d, collapse = ", "),
         "\n  Valid: ", paste(OMNIPATH_DATASETS, collapse = ", "))
  
  if (!as.integer(organism) %in% c(9606L, 10090L, 10116L))
    warning("The service serves taxa 9606, 10090 and 10116; ",
            organism, " will be orthology-translated or rejected.")
  
  url <- sprintf(
    "https://omnipathdb.org/interactions?genesymbols=yes&datasets=%s&organisms=%d&fields=%s&license=%s",
    paste(datasets, collapse = ","), as.integer(organism),
    paste(fields, collapse = ","), license)
  
  message("GET ", url)
  df <- tryCatch(
    utils::read.delim(url, sep = "\t", stringsAsFactors = FALSE,
                      quote = "", comment.char = ""),
    error = function(e) stop("Could not reach the OmniPath service: ",
                             conditionMessage(e)))
  
  # The service answers bad queries with 200 OK and a short HTML/text error
  # page, so a malformed request surfaces as a one-column table rather than an
  # HTTP failure. Catch that explicitly.
  if (ncol(df) <= 2 || !"source_genesymbol" %in% names(df)) {
    stop("OmniPath returned an error page rather than a table.\n",
         "  Response began: ", paste(head(names(df), 3), collapse = " "), "\n",
         "  This normally means an invalid parameter value. Check the valid\n",
         "  values at https://omnipathdb.org/queries/interactions and compare\n",
         "  against the URL above.")
  }
  
  # OmnipathR reports these two, but they are derived, not served. Recompute so
  # downstream code that expects them still works.
  if (!is.null(df$references))
    df$n_references <- vapply(strsplit(as.character(df$references), ";"),
                              function(x) length(unique(x[x != ""])), 0L)
  if (!is.null(df$sources))
    df$n_resources <- vapply(strsplit(as.character(df$sources), ";"),
                             function(x) length(unique(x[x != ""])), 0L)
  
  message(sprintf("%d interactions, %d unique gene symbols",
                  nrow(df),
                  length(unique(c(df$source_genesymbol, df$target_genesymbol)))))
  
  if (!is.null(destfile)) {
    utils::write.table(df, destfile, sep = "\t", row.names = FALSE, quote = FALSE)
    message("wrote ", destfile)
  }
  df
}


#' Build the network graph directly, skipping the intermediate file
#'
#' Convenience wrapper: \code{fetch_omnipath()} then the same cost model as
#' \code{load_omnipath()}.
#'
#' @inheritParams fetch_omnipath
#' @param min_cost,max_cost edge cost range
#' @export
omnipath_graph <- function(datasets = "omnipath", organism = 9606L,
                           min_cost = 0.05, max_cost = 1.0, destfile = NULL) {
  df <- fetch_omnipath(datasets = datasets, organism = organism,
                       destfile = destfile)
  omnipath_df_to_graph(df, min_cost = min_cost, max_cost = max_cost)
}


#' Convert an OmniPath interaction data.frame to a weighted igraph
#'
#' Split out from \code{load_omnipath()} so it works on a data.frame from any
#' source -- the direct download, OmnipathR if it happens to be working, or a
#' file you saved earlier.
#'
#' @export
omnipath_df_to_graph <- function(df,
                                 source_col = "source_genesymbol",
                                 target_col = "target_genesymbol",
                                 confidence_col = "curation_effort",
                                 min_cost = 0.05, max_cost = 1.0) {
  df <- df[!is.na(df[[source_col]]) & !is.na(df[[target_col]]), , drop = FALSE]
  df <- df[df[[source_col]] != "" & df[[target_col]] != "", , drop = FALSE]
  df <- df[df[[source_col]] != df[[target_col]], , drop = FALSE]
  
  # A minority of OmniPath records carry complexes as COMP:... or multi-subunit
  # symbols joined by underscores. These are not gene symbols and will never
  # match a prize, so drop them rather than letting them sit as isolated nodes.
  bad <- grepl("^COMP:", df[[source_col]]) | grepl("^COMP:", df[[target_col]])
  if (any(bad)) {
    message(sprintf("dropping %d complex records", sum(bad)))
    df <- df[!bad, , drop = FALSE]
  }
  
  if (!is.null(df[[confidence_col]])) {
    conf <- rank_normalise(ifelse(is.na(df[[confidence_col]]), 0,
                                  df[[confidence_col]]))
  } else {
    warning(sprintf("No '%s' column; using uniform edge costs.", confidence_col))
    conf <- rep(0.5, nrow(df))
  }
  df$cost <- pmin(pmax(max_cost * (1 - conf), min_cost), max_cost)
  
  ord <- order(df$cost)
  df <- df[ord, , drop = FALSE]
  key <- paste(pmin(df[[source_col]], df[[target_col]]),
               pmax(df[[source_col]], df[[target_col]]), sep = "|")
  df <- df[!duplicated(key), , drop = FALSE]
  
  edges <- data.frame(from = df[[source_col]], to = df[[target_col]],
                      cost = df$cost, stringsAsFactors = FALSE)
  for (a in c("is_directed", "is_stimulation", "is_inhibition")) {
    if (!is.null(df[[a]])) edges[[a]] <- df[[a]]
  }
  g <- igraph::graph_from_data_frame(edges, directed = FALSE)
  message(sprintf("graph: %d nodes, %d edges", igraph::vcount(g), igraph::ecount(g)))
  g
}