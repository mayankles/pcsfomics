## tcga_input.R -- fetch TCGA data and produce pcsfomics prize inputs.
##
## Two layers, deliberately separated:
##   fetch_*   talk to the GDC via TCGAbiolinks. Slow, network-bound, cached.
##   diff_*    pure computation on matrices. Fast, testable, no network.
##
## Keeping them apart means you can iterate on the analysis (which you will do
## many times) without re-downloading (which you want to do once).


# ---------------------------------------------------------------------------
# Barcode utilities
# ---------------------------------------------------------------------------

#' Patient identifier from a TCGA barcode
#'
#' TCGA-A7-A0CE-01A-11R-A00Z-07 -> TCGA-A7-A0CE
tcga_patient <- function(barcode) substr(barcode, 1, 12)

#' Sample type code from a TCGA barcode
#'
#' Characters 14-15. 01 = primary solid tumour, 11 = solid tissue normal,
#' 06 = metastatic. Anything else is usually not what you want in a
#' tumour-vs-normal contrast.
tcga_sample_type <- function(barcode) substr(barcode, 14, 15)

#' Label samples as tumour / normal / other
#'
#' @param barcodes character vector of TCGA barcodes
#' @param tumour_codes,normal_codes sample type codes to accept
#' @return factor with levels normal, tumour, other
tcga_group <- function(barcodes,
                       tumour_codes = c("01"),
                       normal_codes = c("11")) {
  st <- tcga_sample_type(barcodes)
  out <- rep("other", length(st))
  out[st %in% tumour_codes] <- "tumour"
  out[st %in% normal_codes] <- "normal"
  factor(out, levels = c("normal", "tumour", "other"))
}

#' Keep one aliquot per patient per group
#'
#' TCGA frequently has several aliquots for the same patient and sample type.
#' Leaving them in inflates n and violates independence -- the replicates are
#' technical, not biological. Ordering by barcode makes the choice deterministic
#' and reproducible, which matters more than which particular aliquot wins.
dedup_aliquots <- function(barcodes) {
  key <- paste(tcga_patient(barcodes), tcga_sample_type(barcodes))
  ord <- order(barcodes)                       # deterministic tie-break
  keep <- ord[!duplicated(key[ord])]
  barcodes[sort(keep)]                         # restore original order
}


# ---------------------------------------------------------------------------
# Fetch layer (requires TCGAbiolinks + GDC network access)
# ---------------------------------------------------------------------------

#' Download and prepare TCGA RNA-seq counts
#'
#' @param project e.g. "TCGA-LUAD"
#' @param cache_dir directory for the prepared .rds
#' @return SummarizedExperiment of STAR raw counts
fetch_expression <- function(project, cache_dir = "cache") {
  dir.create(cache_dir, showWarnings = FALSE, recursive = TRUE)
  f <- file.path(cache_dir, paste0(project, "_expr.rds"))
  if (file.exists(f)) return(readRDS(f))

  if (!requireNamespace("TCGAbiolinks", quietly = TRUE))
    stop("Install TCGAbiolinks: BiocManager::install('TCGAbiolinks')")

  q <- TCGAbiolinks::GDCquery(
    project = project,
    data.category = "Transcriptome Profiling",
    data.type = "Gene Expression Quantification",
    workflow.type = "STAR - Counts")
  TCGAbiolinks::GDCdownload(q)
  se <- TCGAbiolinks::GDCprepare(q)
  saveRDS(se, f)
  se
}

#' Download and prepare TCGA 450K methylation beta values
#'
#' Note this pulls the harmonised beta-value files, not IDATs. If you want
#' control over normalisation (noob, funnorm) you need the raw IDATs and
#' sesame/minfi instead -- worth considering, since GDC betas come from a fixed
#' pipeline you cannot adjust.
#'
#' @param project e.g. "TCGA-LUAD"
#' @param platform "Illumina Human Methylation 450" or "... 27"
fetch_methylation <- function(project, platform = "Illumina Human Methylation 450",
                              cache_dir = "cache") {
  dir.create(cache_dir, showWarnings = FALSE, recursive = TRUE)
  f <- file.path(cache_dir, paste0(project, "_meth.rds"))
  if (file.exists(f)) return(readRDS(f))

  if (!requireNamespace("TCGAbiolinks", quietly = TRUE))
    stop("Install TCGAbiolinks: BiocManager::install('TCGAbiolinks')")

  q <- TCGAbiolinks::GDCquery(
    project = project,
    data.category = "DNA Methylation",
    data.type = "Methylation Beta Value",
    platform = platform)
  TCGAbiolinks::GDCdownload(q)
  se <- TCGAbiolinks::GDCprepare(q)
  saveRDS(se, f)
  se
}


# ---------------------------------------------------------------------------
# Gene identifier mapping
# ---------------------------------------------------------------------------

#' Map Ensembl gene IDs to HGNC symbols
#'
#' STAR count matrices are indexed by versioned Ensembl IDs
#' (ENSG00000141510.16). OmniPath is indexed by gene symbol. This mismatch is
#' the single most common reason a network integration silently loses half its
#' terminals, so it is worth checking the mapping rate rather than assuming it.
#'
#' The SummarizedExperiment from GDCprepare carries gene_name in rowData; that
#' is preferred over a biomaRt round-trip because it matches the exact
#' annotation build the counts were generated against.
#'
#' @param se SummarizedExperiment, or NULL if supplying ids directly
#' @param ids character vector of Ensembl IDs (used when se is NULL)
#' @return character vector of symbols, NA where unmapped
ensembl_to_symbol <- function(se = NULL, ids = NULL) {
  if (!is.null(se)) {
    rd <- SummarizedExperiment::rowData(se)
    for (col in c("gene_name", "external_gene_name", "hgnc_symbol")) {
      if (!is.null(rd[[col]])) return(as.character(rd[[col]]))
    }
    ids <- rownames(se)
  }
  warning("No symbol column found in rowData; returning stripped Ensembl IDs. ",
          "Map these yourself before matching against OmniPath.")
  sub("\\..*$", "", ids)
}


# ---------------------------------------------------------------------------
# Differential expression
# ---------------------------------------------------------------------------

#' Differential expression, tumour vs normal
#'
#' Uses limma-voom by default. The paired option blocks on patient, which is a
#' large power gain when adjacent normals exist -- but TCGA normals are scarce
#' (often 10-60 per tumour type versus hundreds of tumours), and blocking
#' restricts the analysis to complete pairs. Whether that trade is worth it
#' depends on how many pairs you actually have, which \code{summarise_design}
#' reports before you commit.
#'
#' @param counts integer matrix, genes x samples
#' @param group factor with levels normal, tumour
#' @param patient patient IDs, required when paired = TRUE
#' @param symbols gene symbols aligned to rownames(counts)
#' @param paired block on patient
#' @param min_count filterByExpr threshold
#' @return data.frame: gene, log2FC, padj, plus AveExpr and t
diff_expression <- function(counts, group, patient = NULL, symbols = NULL,
                            paired = FALSE, min_count = 10) {
  stopifnot(ncol(counts) == length(group))
  keep_s <- group %in% c("normal", "tumour")
  counts <- counts[, keep_s, drop = FALSE]
  group <- droplevels(factor(group[keep_s], levels = c("normal", "tumour")))
  if (!is.null(patient)) patient <- patient[keep_s]

  if (paired) {
    if (is.null(patient)) stop("paired = TRUE requires patient IDs.")
    complete <- names(which(tapply(group, patient,
                                   function(x) length(unique(x)) == 2)))
    if (length(complete) < 3)
      stop(sprintf("Only %d complete tumour/normal pairs; use paired = FALSE.",
                   length(complete)))
    sel <- patient %in% complete
    counts <- counts[, sel, drop = FALSE]
    group <- droplevels(group[sel])
    patient <- droplevels(factor(patient[sel]))
    design <- stats::model.matrix(~ patient + group)
    coef_name <- "grouptumour"
  } else {
    design <- stats::model.matrix(~ group)
    coef_name <- "grouptumour"
  }

  d <- edgeR::DGEList(counts = counts, group = group)
  keep <- edgeR::filterByExpr(d, design = design, min.count = min_count)
  d <- d[keep, , keep.lib.sizes = FALSE]
  d <- edgeR::calcNormFactors(d)

  v <- limma::voom(d, design)
  fit <- limma::eBayes(limma::lmFit(v, design))
  tt <- limma::topTable(fit, coef = coef_name, number = Inf, sort.by = "none")

  gene <- if (!is.null(symbols)) symbols[keep] else rownames(d)
  out <- data.frame(gene = as.character(gene),
                    log2FC = tt$logFC,
                    padj = tt$adj.P.Val,
                    AveExpr = tt$AveExpr,
                    t = tt$t,
                    stringsAsFactors = FALSE)
  out <- out[!is.na(out$gene) & out$gene != "", , drop = FALSE]
  # collapse duplicate symbols by strongest effect, matching what the prize
  # function would do anyway -- doing it here keeps the table honest
  out <- out[order(-abs(out$log2FC)), , drop = FALSE]
  out[!duplicated(out$gene), , drop = FALSE]
}


# ---------------------------------------------------------------------------
# Differential methylation
# ---------------------------------------------------------------------------

#' Beta to M-value
#'
#' Testing is done on M-values because beta values are bounded on [0,1] and
#' strongly heteroscedastic at the extremes, which violates the constant-variance
#' assumption of a linear model (Du et al. 2010). Effect sizes are reported back
#' on the beta scale, where "a 30% methylation difference" is interpretable.
beta_to_m <- function(beta, eps = 1e-3) {
  b <- pmin(pmax(beta, eps), 1 - eps)
  log2(b / (1 - b))
}

#' Differential methylation, tumour vs normal, probe level
#'
#' @param beta numeric matrix, probes x samples, values on [0,1]
#' @param group factor with levels normal, tumour
#' @param patient patient IDs, required when paired = TRUE
#' @param paired block on patient
#' @return data.frame: probe, delta_beta, padj, t
diff_methylation <- function(beta, group, patient = NULL, paired = FALSE) {
  stopifnot(ncol(beta) == length(group))
  keep_s <- group %in% c("normal", "tumour")
  beta <- beta[, keep_s, drop = FALSE]
  group <- droplevels(factor(group[keep_s], levels = c("normal", "tumour")))
  if (!is.null(patient)) patient <- patient[keep_s]

  if (paired) {
    if (is.null(patient)) stop("paired = TRUE requires patient IDs.")
    complete <- names(which(tapply(group, patient,
                                   function(x) length(unique(x)) == 2)))
    if (length(complete) < 3)
      stop(sprintf("Only %d complete pairs; use paired = FALSE.", length(complete)))
    sel <- patient %in% complete
    beta <- beta[, sel, drop = FALSE]
    group <- droplevels(group[sel])
    patient <- droplevels(factor(patient[sel]))
    design <- stats::model.matrix(~ patient + group)
  } else {
    design <- stats::model.matrix(~ group)
  }

  ok <- rowSums(is.na(beta)) == 0
  beta <- beta[ok, , drop = FALSE]

  m <- beta_to_m(beta)
  fit <- limma::eBayes(limma::lmFit(m, design))
  tt <- limma::topTable(fit, coef = "grouptumour", number = Inf, sort.by = "none")

  # effect size on the beta scale, always the unpaired group means so the
  # number means the same thing regardless of the design used for testing
  db <- rowMeans(beta[, group == "tumour", drop = FALSE]) -
        rowMeans(beta[, group == "normal", drop = FALSE])

  data.frame(probe = rownames(beta),
             delta_beta = as.numeric(db),
             padj = tt$adj.P.Val,
             t = tt$t,
             stringsAsFactors = FALSE)
}


#' Expand probe-level results to gene x genomic-context rows
#'
#' The 450K annotation packs multiple gene/region assignments into a single
#' semicolon-delimited field, and the two fields are POSITIONALLY ALIGNED --
#' the third entry of UCSC_RefGene_Name goes with the third entry of
#' UCSC_RefGene_Group. Splitting them independently silently scrambles the
#' gene-to-region mapping, which is exactly the annotation the promoter/body
#' split depends on. This function unpacks them pairwise.
#'
#' Duplicate gene/region pairs from one probe (common, since a probe can hit
#' several transcripts of the same gene) are collapsed.
#'
#' @param dm output of \code{diff_methylation}
#' @param anno data.frame with probe, gene_field, region_field columns
#' @return data.frame: gene, region, delta_beta, padj, probe
annotate_probes <- function(dm, anno,
                            probe_col = "probe",
                            gene_field = "UCSC_RefGene_Name",
                            region_field = "UCSC_RefGene_Group") {
  a <- anno[anno[[probe_col]] %in% dm[[probe_col]], , drop = FALSE]
  genes <- strsplit(as.character(a[[gene_field]]), ";", fixed = TRUE)
  regs  <- strsplit(as.character(a[[region_field]]), ";", fixed = TRUE)

  n <- vapply(genes, length, 0L)
  keep <- n > 0 & n == vapply(regs, length, 0L)
  if (any(!keep & n > 0))
    warning(sprintf("%d probes dropped: gene and region fields differ in length.",
                    sum(!keep & n > 0)))

  long <- data.frame(
    probe  = rep(a[[probe_col]][keep], n[keep]),
    gene   = unlist(genes[keep], use.names = FALSE),
    region = unlist(regs[keep], use.names = FALSE),
    stringsAsFactors = FALSE)
  long <- long[long$gene != "" & !is.na(long$gene), , drop = FALSE]
  long <- long[!duplicated(long[, c("probe", "gene", "region")]), , drop = FALSE]

  out <- merge(long, dm, by = "probe")
  out[, c("gene", "region", "delta_beta", "padj", "probe")]
}


#' Drop probes that should not inform a differential analysis
#'
#' Sex chromosomes confound any cohort with unbalanced sex. Cross-reactive
#' probes map to multiple genomic locations and produce signal that does not
#' belong to the annotated gene. SNP-overlapping probes report genotype rather
#' than methylation. Chen et al. (2013) and Zhou et al. (2017) publish the
#' standard exclusion lists; pass them in via \code{bad_probes}.
filter_probes <- function(beta, anno, probe_col = "probe",
                          chr_col = "chr", bad_probes = character(0),
                          drop_sex = TRUE) {
  keep <- rownames(beta)
  if (drop_sex && !is.null(anno[[chr_col]])) {
    sex <- anno[[probe_col]][anno[[chr_col]] %in% c("chrX", "chrY", "X", "Y")]
    keep <- setdiff(keep, sex)
  }
  keep <- setdiff(keep, bad_probes)
  beta[rownames(beta) %in% keep, , drop = FALSE]
}


# ---------------------------------------------------------------------------
# Design reporting
# ---------------------------------------------------------------------------

#' Report the design before committing to it
#'
#' Prints tumour/normal counts per modality and, critically, the number of
#' patients with BOTH assays. For a project built on cross-modal integration,
#' that overlap is the real sample size -- not the per-assay totals.
#'
#' @export
summarise_design <- function(expr_barcodes, meth_barcodes) {
  eg <- tcga_group(expr_barcodes); mg <- tcga_group(meth_barcodes)
  ep <- tcga_patient(expr_barcodes); mp <- tcga_patient(meth_barcodes)

  cat("expression : ", sum(eg == "tumour"), " tumour, ",
      sum(eg == "normal"), " normal\n", sep = "")
  cat("methylation: ", sum(mg == "tumour"), " tumour, ",
      sum(mg == "normal"), " normal\n", sep = "")

  et <- unique(ep[eg == "tumour"]); mt <- unique(mp[mg == "tumour"])
  en <- unique(ep[eg == "normal"]); mn <- unique(mp[mg == "normal"])
  cat("patients with both assays (tumour): ", length(intersect(et, mt)), "\n", sep = "")
  cat("patients with both assays (normal): ", length(intersect(en, mn)), "\n", sep = "")

  pe <- intersect(unique(ep[eg == "tumour"]), unique(ep[eg == "normal"]))
  pm <- intersect(unique(mp[mg == "tumour"]), unique(mp[mg == "normal"]))
  cat("complete tumour/normal pairs -- expression: ", length(pe),
      ", methylation: ", length(pm), "\n", sep = "")

  invisible(list(expr_pairs = pe, meth_pairs = pm,
                 both_tumour = intersect(et, mt)))
}
