## run_tcga_pipeline.R -- TCGA download through to a PCSF subnetwork.
##
##   Rscript inst/examples/run_tcga_pipeline.R TCGA-LUAD
##
## Requires TCGAbiolinks and network access to the GDC. The first run downloads
## several GB and takes a while; results are cached as .rds in cache/, so
## re-runs are fast.
##
## Strategy note: run this on a large TCGA cohort FIRST to choose beta, omega
## and mu, then freeze those values and apply them to your small cohort.
## Selecting parameters on ~30 samples will overfit, and you will not be able to
## tell that it did.

# Installs required
# install.packages("BiocManager")
# install.packages("igraph")
# install.packages("stringdist")
# BiocManager::install(c("limma","edgeR","SummarizedExperiment"))
# BiocManager::install(c("sesame","sesameData"))

rm(list=ls())
setwd("~/Documents/GitHub/pcsfomics")

suppressMessages({
  library(limma); library(edgeR); library(igraph)
  library(SummarizedExperiment)
})
for (f in list.files("R", pattern="*.R",full.names = TRUE)) source(f)

args <- commandArgs(trailingOnly = TRUE)
PROJECT  <- if (length(args) >= 1) args[1] else "TCGA-THCA"
OUT      <- if (length(args) >= 2) args[2] else "results"
PAIRED   <- TRUE      # block on patient where adjacent normals exist
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

# ---------------------------------------------------------------------------
# Thresholds (the knobs to sweep later).
#
# Expression previously used no fold-change floor (lfc_cutoff = 0), so every
# FDR-significant gene became a prize terminal -- ~19k on a large TCGA cohort,
# where tumour-vs-normal genuinely shifts most of the transcriptome. That count
# is what breaks the downstream solver: the KMB metric closure is O(terminals^2)
# and, worse, permutation_test relocates prizes onto fresh uncached nodes each
# permutation, so cost scales with terminals x permutations. An effect-size
# floor on expression cuts terminals to the low thousands and makes both the
# solve and the null test tractable.
#
# These are surfaced here, not buried in library defaults, because the plan is
# to vary them and choose the most parsimonious / stable network. Keeping them
# in one block makes that sweep a loop over this script's inputs.
DE_PADJ    <- 0.05
DE_LFC     <- 1.0     # |log2FC| floor; 1.0 == 2-fold. Was effectively 0.
METH_PADJ  <- 0.05
METH_DELTA <- 0.1     # |delta beta| floor (methylation_scores default)

# OmniPath is a protein-coding signalling network, so non-coding DEGs (lncRNA,
# pseudogenes, small RNAs -- ~40% of the significant list on THCA) can never map
# into it. They only dilute the rank-normalised prizes and inflate the terminal
# count. Restricting expression evidence to protein_coding is therefore correct
# FOR THIS NETWORK. It is a pipeline knob, not a library default, because other
# networks (planned) may warrant a different biotype scope. Coverage is reported
# by biotype every run (biotype_coverage.tsv) so the effect stays monitorable.
PROTEIN_CODING_ONLY <- TRUE

data_dir=".data"
if (!dir.exists(data_dir)){dir.create(data_dir, recursive=T)}
message("project: ", PROJECT)

# ---------------------------------------------------------------------------
# Timing. Every expensive stage is wrapped so it is obvious where wall-clock
# goes -- the permutation and robustness steps scale with terminal count, so
# their timing is also the signal for whether thresholds need tightening.
# ---------------------------------------------------------------------------
.timings <- list()
timeit <- function(label, expr) {
  t0 <- Sys.time()
  val <- force(expr)
  dt  <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  .timings[[label]] <<- dt
  message(sprintf("[timing] %-26s %9.1f s", label, dt))
  val
}

# ---------------------------------------------------------------------------
# 1. Fetch
# ---------------------------------------------------------------------------
options(timeout = 3000)
se_expr <- timeit("fetch expression", fetch_expression(PROJECT, cache_dir = data_dir))
se_meth <- timeit("fetch methylation", fetch_methylation(PROJECT, cache_dir = data_dir))

message("\n--- design ---")
design_info <- summarise_design(colnames(se_expr), colnames(se_meth))

# ---------------------------------------------------------------------------
# 2. Differential expression
# ---------------------------------------------------------------------------
message("\n--- differential expression ---")
expr_bc <- dedup_aliquots(colnames(se_expr))
se_expr <- se_expr[, expr_bc]

counts  <- SummarizedExperiment::assay(se_expr, "unstranded")
symbols <- ensembl_to_symbol(se_expr)
grp     <- tcga_group(expr_bc)
pat     <- tcga_patient(expr_bc)

# fall back to unpaired if there are too few complete pairs to block on
n_pairs <- length(design_info$expr_pairs)
use_paired <- PAIRED && n_pairs >= 10
if (PAIRED && !use_paired)
  message("only ", n_pairs, " complete pairs; falling back to unpaired")

de <- timeit("diff_expression",
             diff_expression(counts, grp, pat, symbols = symbols, paired = use_paired))
# Console message: "calcNormFactors has been renamed to normLibSizes"

# Genuinely many genes clear FDR at TCGA sample sizes -- tumour-vs-normal shifts
# most of the transcriptome, so FDR alone is not "too many", it is expected. The
# fold-change floor below (DE_LFC) is what keeps the terminal set sane. Report
# both counts so the effect of the floor is visible.
n_sig  <- sum(de$padj < DE_PADJ, na.rm = TRUE)
n_sigl <- sum(de$padj < DE_PADJ & abs(de$log2FC) >= DE_LFC, na.rm = TRUE)
message(nrow(de), " genes tested, ", n_sig, " at FDR < ", DE_PADJ,
        ", ", n_sigl, " also |log2FC| >= ", DE_LFC)

# Biotype map from the STAR SummarizedExperiment (gene_name -> gene_type). Used
# both to restrict prizes to protein_coding and to report coverage by biotype.
biotype <- stats::setNames(
  as.character(SummarizedExperiment::rowData(se_expr)$gene_type),
  as.character(SummarizedExperiment::rowData(se_expr)$gene_name))
de$biotype <- biotype[de$gene]

sig_mask  <- de$padj < DE_PADJ & abs(de$log2FC) >= DE_LFC & !is.na(de$padj)
sig_bt    <- ifelse(is.na(de$biotype[sig_mask]), "unknown", de$biotype[sig_mask])
n_pc_sig  <- sum(sig_bt == "protein_coding")
message(sprintf("significant DEGs: %d total, %d protein_coding (%.0f%%), %d non-coding/other",
                length(sig_bt), n_pc_sig, 100 * n_pc_sig / length(sig_bt),
                length(sig_bt) - n_pc_sig))

utils::write.table(de, file.path(OUT, "de_table.tsv"), sep = "\t",
                   row.names = FALSE, quote = FALSE)

# ---------------------------------------------------------------------------
# 3. Differential methylation
# ---------------------------------------------------------------------------
message("\n--- differential methylation ---")
meth_bc <- dedup_aliquots(colnames(se_meth))
se_meth <- se_meth[, meth_bc]

beta  <- SummarizedExperiment::assay(se_meth)
mgrp  <- tcga_group(meth_bc)
mpat  <- tcga_patient(meth_bc)

# Probe annotation. GDCprepare attaches gene and position columns to rowData;
# the IlluminaHumanMethylation450kanno package is richer (CpG island context,
# which is worth stratifying on) if you have it installed.
rd <- as.data.frame(SummarizedExperiment::rowData(se_meth))
anno <- data.frame(
  probe = rownames(se_meth),
  UCSC_RefGene_Name  = if (!is.null(rd$Gene_Symbol)) rd$Gene_Symbol else rd$gene,
  # UCSC_RefGene_Group = if (!is.null(rd$Feature_Type)) rd$Feature_Type else NA,
  UCSC_RefGene_Group = if (!is.null(rd$Feature_Type)) rd$Feature_Type else NA,
  chr = if (!is.null(rd$chrm_A)) rd$chrm_A else rd$seqnames,
  stringsAsFactors = FALSE)

# The GDC/sesame rowData carries gene symbols but NO region field (no
# UCSC_RefGene_Group / Feature_Type), so the fallback above leaves region = NA
# for every probe -- which makes the promoter/body split, and therefore the
# concordance term, score ZERO. The 450k annotation package supplies the region
# labels; they are per-probe and join by cg-ID, so the anno package being hg19
# while the data is hg38 does not matter for the LABELS.
if (requireNamespace("IlluminaHumanMethylation450kanno.ilmn12.hg19", quietly = TRUE) &&
    requireNamespace("minfi", quietly = TRUE)) {
  message("using IlluminaHumanMethylation450kanno for region annotation")
  # getAnnotation() -> updateObject() looks the package up on the search list, so
  # it must be ATTACHED with library(), not merely referenced with ::. The ::
  # form errors with "no item called package:IlluminaHumanMethylation450kanno...".
  suppressMessages({
    library(minfi)
    library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
  })
  ann450k <- minfi::getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)
  anno <- data.frame(
    probe = rownames(ann450k),
    UCSC_RefGene_Name  = ann450k$UCSC_RefGene_Name,
    UCSC_RefGene_Group = ann450k$UCSC_RefGene_Group,
    chr    = ann450k$chr,
    island = ann450k$Relation_to_UCSC_CpG_Island,
    stringsAsFactors = FALSE)
}

# Guard against the silent-zero failure: if no probe carries a region label, the
# promoter/body split cannot fire and methylation contributes nothing. Fail loud
# rather than producing an expression-only result that looks multi-modal.
if (all(is.na(anno$UCSC_RefGene_Group) | anno$UCSC_RefGene_Group == ""))
  stop("No probe region annotation (UCSC_RefGene_Group is entirely empty). ",
       "The GDC/sesame rowData has no region field; install ",
       "IlluminaHumanMethylation450kanno.ilmn12.hg19 + minfi so the promoter/",
       "body split -- and the concordance term -- can score.")

# Cross-reactive and SNP-overlapping probes should also be excluded here.
# Chen et al. 2013 and Zhou et al. 2017 publish the lists; load one and pass it
# as bad_probes. Skipping this inflates apparent signal at affected loci.
beta <- filter_probes(beta, anno, drop_sex = TRUE, bad_probes = character(0))
message(nrow(beta), " probes after filtering")

dm <- timeit("diff_methylation",
             diff_methylation(beta, mgrp, mpat,
                              paired = PAIRED && length(design_info$meth_pairs) >= 10))
dmr <- timeit("annotate_probes", annotate_probes(dm, anno))
# Warning message:
#   In annotate_probes(dm, anno) :
#   91063 probes dropped: gene and region fields differ in length.
message(nrow(dmr), " gene x region rows, ",
        sum(dmr$padj < 0.05 & abs(dmr$delta_beta) > 0.1), " passing thresholds")
# 231770 gene x region rows, 8755 passing thresholds
utils::write.table(dmr, file.path(OUT, "dmr_table.tsv"), sep = "\t",
                   row.names = FALSE, quote = FALSE)

# ---------------------------------------------------------------------------
# 4. Prizes and PCSF
# ---------------------------------------------------------------------------
message("\n--- prizes ---")

# Restrict the expression evidence to protein_coding BEFORE scoring, so
# rank-normalisation ranks each gene only against genes that can actually
# participate in the network -- non-coding genes would otherwise shift the
# ranks of the protein-coding genes that matter.
de_expr <- de
if (PROTEIN_CODING_ONLY) {
  keep_pc <- !is.na(de_expr$biotype) & de_expr$biotype == "protein_coding"
  message(sprintf("protein-coding restriction (expression): %d -> %d genes",
                  nrow(de_expr), sum(keep_pc)))
  de_expr <- de_expr[keep_pc, , drop = FALSE]
}

expr_s <- expression_scores(de_expr, padj_cutoff = DE_PADJ, lfc_cutoff = DE_LFC)
meth_s <- methylation_scores(dmr, padj_cutoff = METH_PADJ, min_delta = METH_DELTA)
if (nrow(meth_s) == 0L)
  warning("methylation_scores returned 0 rows: no probes matched the promoter/",
          "body region vocabulary, so methylation and the concordance term ",
          "contribute nothing. Check that dmr$region carries UCSC_RefGene_Group ",
          "labels (TSS200, Body, ...), not NA.")
message(sprintf("methylation scored: %d gene rows (promoter and/or body)", nrow(meth_s)))
pr_df  <- build_prizes(expr_s, meth_s, prize_config(beta = 1))
utils::write.table(pr_df, file.path(OUT, "prizes.tsv"), sep = "\t",
                   row.names = FALSE, quote = FALSE)

# Terminal count is the quantity that drives solver cost. Report it here, before
# touching the network, so the effect of the thresholds is immediate.
n_term <- sum(pr_df$prize > 0)
message(sprintf("prize terminals: %d nonzero of %d genes (DE_LFC=%.2f, DE_PADJ=%.3f)",
                n_term, nrow(pr_df), DE_LFC, DE_PADJ))

if (!file.exists("omnipath.tsv")) {
  message("\nomnipath.tsv not found. Create it with:")
  message('  library(OmnipathR)')
  message('  ia <- import_omnipath_interactions(organism = 9606L)')
  # message('  ia <- import_all_interactions()')
  message('  write.table(ia, "omnipath.tsv", sep="\\t", row.names=FALSE, quote=FALSE)')
  # quit(status = 0)
}

g  <- load_omnipath("omnipath.tsv")
pv <- prize_vector(pr_df)

# Terminal coverage is the check people skip and then wonder why the subnetwork
# is small. If most prized genes are absent from the network, the problem is
# identifier mismatch, not biology.
hit <- sum(names(pv) %in% V(g)$name)
message(sprintf("terminals: %d prized genes, %d present in network (%.0f%%)",
                length(pv), hit, 100 * hit / length(pv)))

# Coverage by biotype, over the FULL significant DEG list (pre-restriction), so
# the split between "structurally cannot map" (non-coding) and "protein-coding
# but missing" (the real identifier-mismatch signal) stays visible each run.
bcov <- biotype_coverage(de$gene[sig_mask], V(g)$name, biotype)
utils::write.table(bcov, file.path(OUT, "biotype_coverage.tsv"), sep = "\t",
                   row.names = FALSE, quote = FALSE)
message("biotype coverage of significant DEGs vs network:")
print(utils::head(bcov, 8), row.names = FALSE)
pc_cov <- bcov$coverage[bcov$biotype == "protein_coding"]
if (length(pc_cov) && pc_cov < 0.4)
  warning(sprintf("Only %.0f%% of protein-coding DEGs map to the network. ",
                  100 * pc_cov),
          "That is the identifier-mismatch tripwire -- check gene symbols ",
          "(Ensembl IDs, aliases, or a stale annotation build).")

om  <- suggest_omega(pv, 0.90)
res <- timeit("solve_pcsf", solve_pcsf(g, pv, omega = om, mu = 0.005))
cs  <- component_summary(res$subnetwork)
message(sprintf("subnetwork: %d nodes (%d terminal, %d Steiner), %d multi-node components",
                vcount(res$subnetwork), res$info$terminals, res$info$steiner,
                cs$n_multinode))

message("\n--- robustness and null ---")
ef <- timeit("robust_pcsf",
             robust_pcsf(g, pv, omega = om, mu = 0.005, n_runs = 100))
utils::write.table(ef, file.path(OUT, "edge_frequency.tsv"), sep = "\t",
                   row.names = FALSE, quote = FALSE)
message(sum(ef$frequency >= 0.9), " of ", nrow(ef), " edges stable at freq >= 0.9")

# Previously stalled overnight here: with ~19k terminals, each permutation
# relocates prizes onto fresh uncached nodes, so cost was terminals x perms
# full-graph Dijkstras. The fold-change floor above is the fix; the timing below
# confirms it. If this is still slow, drop n_perm or tighten DE_LFC further.
perm <- timeit("permutation_test",
               permutation_test(g, pv, omega = om, mu = 0.005, n_perm = 200))
message(sprintf("p_objective = %.4f (primary), p_compactness = %.4f",
                perm$p_objective, perm$p_compactness))

saveRDS(list(subnetwork = res$subnetwork, edge_freq = ef, perm = perm,
             params = list(beta = 1, omega = om, mu = 0.005,
                           DE_PADJ = DE_PADJ, DE_LFC = DE_LFC,
                           METH_PADJ = METH_PADJ, METH_DELTA = METH_DELTA,
                           n_terminals = n_term)),
        file.path(OUT, "pcsf_result.rds"))

# Timing summary -- persisted so a threshold sweep can compare runs on cost as
# well as on network quality.
tdf <- data.frame(stage = names(.timings),
                  seconds = unlist(.timings, use.names = FALSE),
                  row.names = NULL, stringsAsFactors = FALSE)
utils::write.table(tdf, file.path(OUT, "timings.tsv"), sep = "\t",
                   row.names = FALSE, quote = FALSE)
message("\n--- timings (s) ---")
print(tdf, row.names = FALSE)
message(sprintf("total wall clock: %.1f s", sum(tdf$seconds)))

message("\nwrote results to ", OUT, "/")
