#' Rank-normalise a numeric vector to (0, 1]
#'
#' Rank-normalisation decouples prize magnitude from cohort size and assay
#' platform. This is the property that makes subnetworks from a 30-sample
#' cohort and a 500-sample TCGA cohort comparable; \code{-log10(p)} does not
#' have it, because p-values scale with n. Preserve this if you modify the
#' prize function.
#'
#' @param x numeric vector
#' @return numeric vector on (0, 1], largest input mapped to 1
#' @export
rank_normalise <- function(x) {
  if (length(x) == 0L) return(numeric(0))
  r <- rank(-x, ties.method = "average")
  (length(x) - r + 1) / length(x)
}


#' Rank-based expression evidence
#'
#' @param de data.frame with gene / effect size / adjusted p-value columns
#' @param gene_col,stat_col,padj_col column names
#' @param padj_cutoff,lfc_cutoff significance and effect-size thresholds
#' @return data.frame: gene, expr_score, expr_direction
#' @export
expression_scores <- function(de,
                              gene_col = "gene",
                              stat_col = "log2FC",
                              padj_col = "padj",
                              padj_cutoff = 0.05,
                              lfc_cutoff = 0) {
  d <- de[stats::complete.cases(de[, c(gene_col, stat_col, padj_col)]), , drop = FALSE]
  d <- d[d[[padj_col]] <= padj_cutoff & abs(d[[stat_col]]) >= lfc_cutoff, , drop = FALSE]
  if (nrow(d) == 0L) {
    warning("No genes passed expression thresholds.")
    return(data.frame(gene = character(0), expr_score = numeric(0),
                      expr_direction = numeric(0), stringsAsFactors = FALSE))
  }
  # collapse duplicate symbols by strongest effect
  d <- d[order(-abs(d[[stat_col]])), , drop = FALSE]
  d <- d[!duplicated(d[[gene_col]]), , drop = FALSE]

  data.frame(
    gene           = as.character(d[[gene_col]]),
    expr_score     = rank_normalise(abs(d[[stat_col]])),
    expr_direction = sign(d[[stat_col]]),
    stringsAsFactors = FALSE
  )
}


#' Rank-based methylation evidence, split by genomic context
#'
#' Promoter and gene-body methylation are scored separately because their
#' expected relationship to expression has opposite sign. Pooling them cancels
#' signal. The default label sets match the \code{UCSC_RefGene_Group} annotation
#' in \code{IlluminaHumanMethylation450kanno.ilmn12.hg19}.
#'
#' Within a gene x context, aggregation uses the maximum |delta beta| rather
#' than a sum. Summing rewards genes covered by many probes, which tracks CpG
#' density rather than regulatory relevance.
#'
#' @param dmr data.frame with gene / region / delta beta / adjusted p columns
#' @param gene_col,region_col,delta_col,padj_col column names
#' @param padj_cutoff,min_delta thresholds
#' @param promoter_labels,body_labels region label vocabularies
#' @return data.frame: gene, prom_score, prom_direction, body_score, body_direction
#' @export
methylation_scores <- function(dmr,
                               gene_col = "gene",
                               region_col = "region",
                               delta_col = "delta_beta",
                               padj_col = "padj",
                               padj_cutoff = 0.05,
                               min_delta = 0.1,
                               promoter_labels = c("TSS1500", "TSS200", "5'UTR",
                                                   "1stExon", "promoter"),
                               body_labels = c("Body", "3'UTR", "ExonBnd", "body")) {
  d <- dmr[stats::complete.cases(dmr[, c(gene_col, delta_col, padj_col)]), , drop = FALSE]
  d <- d[d[[padj_col]] <= padj_cutoff & abs(d[[delta_col]]) >= min_delta, , drop = FALSE]

  one <- function(labels, prefix) {
    s <- d[d[[region_col]] %in% labels, , drop = FALSE]
    if (nrow(s) == 0L) {
      out <- data.frame(gene = character(0), a = numeric(0), b = numeric(0),
                        stringsAsFactors = FALSE)
    } else {
      s <- s[order(-abs(s[[delta_col]])), , drop = FALSE]
      s <- s[!duplicated(s[[gene_col]]), , drop = FALSE]
      out <- data.frame(gene = as.character(s[[gene_col]]),
                        a = rank_normalise(abs(s[[delta_col]])),
                        b = sign(s[[delta_col]]),
                        stringsAsFactors = FALSE)
    }
    names(out) <- c("gene", paste0(prefix, "_score"), paste0(prefix, "_direction"))
    out
  }

  merge(one(promoter_labels, "prom"), one(body_labels, "body"),
        by = "gene", all = TRUE)
}


#' Concordance bonus for mechanistically coherent multi-modal signal
#'
#' This is where the integration actually happens. Adding expression and
#' methylation evidence together treats a gene with strong-but-incoherent signal
#' the same as a gene whose two modalities tell a consistent story. This term
#' does not.
#'
#' Promoter hypermethylation with downregulation (or the reverse) is the
#' canonical silencing axis and receives full credit. Gene-body hypermethylation
#' with upregulation is the expected direction but rests on weaker evidence, so
#' it receives partial credit.
#'
#' The bonus is scaled by the geometric mean of the two contributing scores, so
#' a gene earns a large bonus only when both modalities are strong AND coherent.
#'
#' @param m merged score data.frame
#' @param promoter_weight,body_weight relative credit
#' @return numeric vector
#' @export
concordance_bonus <- function(m, promoter_weight = 1, body_weight = 0.5) {
  z <- function(v) { v <- m[[v]]; if (is.null(v)) 0 else ifelse(is.na(v), 0, v) }
  es <- z("expr_score"); ed <- z("expr_direction")

  prom <- pmax(pmin(-sign(z("prom_direction") * ed), 1), 0) *
    sqrt(z("prom_score") * es)
  body <- pmax(pmin(sign(z("body_direction") * ed), 1), 0) *
    sqrt(z("body_score") * es)

  promoter_weight * prom + body_weight * body
}


#' Prize configuration
#'
#' \code{beta} scales prizes against edge costs and is the single most
#' consequential numeric parameter -- but note that in benchmarking, changing
#' the prize *definition* (the weights) moved precision far more than sweeping
#' beta did.
#'
#' @export
prize_config <- function(w_expr = 1, w_prom = 1, w_body = 0.4,
                         w_conc = 1.5, beta = 1) {
  list(w_expr = w_expr, w_prom = w_prom, w_body = w_body,
       w_conc = w_conc, beta = beta)
}


#' Combine modality scores into a prize vector
#'
#' @param expr output of \code{expression_scores}
#' @param meth output of \code{methylation_scores}
#' @param config output of \code{prize_config}
#' @return data.frame sorted by prize, with all component scores retained
#' @export
build_prizes <- function(expr, meth, config = prize_config()) {
  m <- merge(expr, meth, by = "gene", all = TRUE)
  m$concordance <- concordance_bonus(m)

  z <- function(v) { v <- m[[v]]; if (is.null(v)) 0 else ifelse(is.na(v), 0, v) }
  m$prize <- config$beta * (
    config$w_expr * z("expr_score") +
    config$w_prom * z("prom_score") +
    config$w_body * z("body_score") +
    config$w_conc * m$concordance
  )
  m[order(-m$prize), , drop = FALSE]
}


#' Named prize vector, dropping zeros
#'
#' The format \code{PCSF::PCSF()} expects for its \code{terminals} argument.
#' @export
prize_vector <- function(prizes_df) {
  v <- prizes_df$prize
  names(v) <- prizes_df$gene
  v[v > 0]
}
