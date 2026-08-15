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

data_dir=".data"
if (!dir.exists(data_dir)){dir.create(data_dir, recursive=T)}
message("project: ", PROJECT)

# ---------------------------------------------------------------------------
# 1. Fetch
# ---------------------------------------------------------------------------
options(timeout = 3000)
se_expr <- fetch_expression(PROJECT, cache_dir = data_dir)
se_meth <- fetch_methylation(PROJECT, cache_dir = data_dir)

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

de <- diff_expression(counts, grp, pat, symbols = symbols, paired = use_paired)
# Console message: "calcNormFactors has been renamed to normLibSizes"

message(nrow(de), " genes tested, ", sum(de$padj < 0.05), " at FDR < 0.05")
utils::write.table(de, file.path(OUT, "de_table.tsv"), sep = "\t",
                   row.names = FALSE, quote = FALSE)

# 29709 genes tested, 19440 at FDR < 0.05
# This is way too many DEGs, sus

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

if (requireNamespace("IlluminaHumanMethylation450kanno.ilmn12.hg19", quietly = TRUE)) {
  message("using IlluminaHumanMethylation450kanno for region annotation")
  ann450k <- minfi::getAnnotation(
    IlluminaHumanMethylation450kanno.ilmn12.hg19::IlluminaHumanMethylation450kanno.ilmn12.hg19)
  anno <- data.frame(
    probe = rownames(ann450k),
    UCSC_RefGene_Name  = ann450k$UCSC_RefGene_Name,
    UCSC_RefGene_Group = ann450k$UCSC_RefGene_Group,
    chr    = ann450k$chr,
    island = ann450k$Relation_to_UCSC_CpG_Island,
    stringsAsFactors = FALSE)
}

# Cross-reactive and SNP-overlapping probes should also be excluded here.
# Chen et al. 2013 and Zhou et al. 2017 publish the lists; load one and pass it
# as bad_probes. Skipping this inflates apparent signal at affected loci.
beta <- filter_probes(beta, anno, drop_sex = TRUE, bad_probes = character(0))
message(nrow(beta), " probes after filtering")

dm <- diff_methylation(beta, mgrp, mpat,
                       paired = PAIRED && length(design_info$meth_pairs) >= 10)
dmr <- annotate_probes(dm, anno)
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
expr_s <- expression_scores(de)
meth_s <- methylation_scores(dmr)
pr_df  <- build_prizes(expr_s, meth_s, prize_config(beta = 1))
utils::write.table(pr_df, file.path(OUT, "prizes.tsv"), sep = "\t",
                   row.names = FALSE, quote = FALSE)

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
if (hit / length(pv) < 0.5)
  warning("Under half the prized genes map to the network. Check gene symbols ",
          "-- Ensembl IDs, aliases, or a stale annotation build are the usual cause.")

om  <- suggest_omega(pv, 0.90)
res <- solve_pcsf(g, pv, omega = om, mu = 0.005)
cs  <- component_summary(res$subnetwork)
message(sprintf("subnetwork: %d nodes (%d terminal, %d Steiner), %d multi-node components",
                vcount(res$subnetwork), res$info$terminals, res$info$steiner,
                cs$n_multinode))

message("\n--- robustness and null ---")
ef <- robust_pcsf(g, pv, omega = om, mu = 0.005, n_runs = 100)
utils::write.table(ef, file.path(OUT, "edge_frequency.tsv"), sep = "\t",
                   row.names = FALSE, quote = FALSE)
message(sum(ef$frequency >= 0.9), " of ", nrow(ef), " edges stable at freq >= 0.9")

perm <- permutation_test(g, pv, omega = om, mu = 0.005, n_perm = 200)
message(sprintf("p_objective = %.4f (primary), p_compactness = %.4f",
                perm$p_objective, perm$p_compactness))

saveRDS(list(subnetwork = res$subnetwork, edge_freq = ef, perm = perm,
             params = list(beta = 1, omega = om, mu = 0.005)),
        file.path(OUT, "pcsf_result.rds"))
message("\nwrote results to ", OUT, "/")
