## Verify the analysis layer on simulated data shaped like TCGA output.
##
## This exercises everything except the GDC download: barcode parsing, aliquot
## deduplication, paired and unpaired designs, voom DE, M-value DM, the
## pairwise-aligned annotation unpacking, and the handoff into pcsfomics.
##
## Planted truth: 40 genes are downregulated in tumour AND promoter
## hypermethylated. A successful run recovers them with high concordance.

suppressMessages({library(limma); library(edgeR); library(igraph)})
for (f in list.files("R", full.names = TRUE)) source(f)

set.seed(11)

N_PATIENT <- 60
N_PAIRED  <- 15          # patients with an adjacent normal, as in real TCGA
N_GENE    <- 3000
N_PROBE   <- 6000
N_TRUE    <- 40

## --- barcodes ------------------------------------------------------------
pid <- sprintf("TCGA-XX-%04d", seq_len(N_PATIENT))
tumour_bc <- paste0(pid, "-01A-11R-A00Z-07")
normal_bc <- paste0(pid[seq_len(N_PAIRED)], "-11A-11R-A00Z-07")
# a duplicate aliquot for one patient, to exercise dedup
dup_bc <- paste0(pid[1], "-01B-11R-A00Z-07")
expr_bc <- c(tumour_bc, normal_bc, dup_bc)

cat("=== barcode parsing ===\n")
stopifnot(tcga_patient(tumour_bc[1]) == pid[1])
stopifnot(tcga_sample_type(normal_bc[1]) == "11")
g <- tcga_group(expr_bc)
cat("groups:", sum(g == "tumour"), "tumour,", sum(g == "normal"), "normal\n")
dd <- dedup_aliquots(expr_bc)
cat("dedup removed", length(expr_bc) - length(dd), "aliquot (expected 1)\n")
stopifnot(length(expr_bc) - length(dd) == 1)
expr_bc <- dd
g <- tcga_group(expr_bc); pt <- tcga_patient(expr_bc)

## --- simulated counts ----------------------------------------------------
symbols <- sprintf("GENE%04d", seq_len(N_GENE))
ens <- sprintf("ENSG%011d.%d", seq_len(N_GENE), sample(1:20, N_GENE, TRUE))
true_genes <- symbols[seq_len(N_TRUE)]

base_mu <- 2^runif(N_GENE, 4, 11)
counts <- matrix(rnbinom(N_GENE * length(expr_bc), mu = base_mu, size = 6),
                 nrow = N_GENE, dimnames = list(ens, expr_bc))
# planted downregulation in tumours
counts[seq_len(N_TRUE), g == "tumour"] <-
  matrix(rnbinom(N_TRUE * sum(g == "tumour"),
                 mu = base_mu[seq_len(N_TRUE)] * 0.25, size = 6),
         nrow = N_TRUE)

cat("\n=== differential expression ===\n")
de_unpaired <- diff_expression(counts, g, pt, symbols = symbols, paired = FALSE)
cat("unpaired: ", nrow(de_unpaired), " genes tested, ",
    sum(de_unpaired$padj < 0.05), " significant\n", sep = "")
rec <- sum(head(de_unpaired$gene[order(de_unpaired$padj)], N_TRUE) %in% true_genes)
cat("planted genes in top ", N_TRUE, ": ", rec, "/", N_TRUE, "\n", sep = "")
stopifnot(rec >= N_TRUE * 0.8)

de_paired <- diff_expression(counts, g, pt, symbols = symbols, paired = TRUE)
cat("paired  : ", nrow(de_paired), " genes tested, ",
    sum(de_paired$padj < 0.05), " significant (n=", N_PAIRED, " pairs)\n", sep = "")

## --- simulated methylation ----------------------------------------------
probes <- sprintf("cg%08d", seq_len(N_PROBE))
beta <- matrix(rbeta(N_PROBE * length(expr_bc), 2, 5), nrow = N_PROBE,
               dimnames = list(probes, expr_bc))

# probe annotation, including multi-gene probes with aligned fields
gene_field <- character(N_PROBE); region_field <- character(N_PROBE)
chr <- sample(c(paste0("chr", 1:22), "chrX"), N_PROBE, TRUE)

# first 2*N_TRUE probes are promoters of the planted genes
promoter_probes <- seq_len(2 * N_TRUE)
gene_field[promoter_probes] <- rep(true_genes, each = 2)
region_field[promoter_probes] <- rep(c("TSS200", "TSS1500"), N_TRUE)
chr[promoter_probes] <- "chr1"

rest <- setdiff(seq_len(N_PROBE), promoter_probes)
for (i in rest) {
  k <- sample(1:3, 1)
  gs <- sample(symbols, k)
  rs <- sample(c("TSS200", "TSS1500", "Body", "5'UTR", "3'UTR"), k, TRUE)
  gene_field[i] <- paste(gs, collapse = ";")
  region_field[i] <- paste(rs, collapse = ";")
}
anno <- data.frame(probe = probes, UCSC_RefGene_Name = gene_field,
                   UCSC_RefGene_Group = region_field, chr = chr,
                   stringsAsFactors = FALSE)

# planted hypermethylation at those promoters
beta[promoter_probes, g == "tumour"] <-
  pmin(beta[promoter_probes, g == "tumour"] + runif(1, .3, .35), 0.99)

cat("\n=== probe filtering ===\n")
beta_f <- filter_probes(beta, anno, drop_sex = TRUE)
cat("dropped ", nrow(beta) - nrow(beta_f), " sex-chromosome probes, ",
    nrow(beta_f), " remain\n", sep = "")
stopifnot(nrow(beta_f) < nrow(beta))

cat("\n=== differential methylation ===\n")
dm <- diff_methylation(beta_f, g, pt, paired = FALSE)
cat(nrow(dm), " probes tested, ", sum(dm$padj < 0.05 & abs(dm$delta_beta) > 0.1),
    " significant with |delta beta| > 0.1\n", sep = "")

cat("\n=== annotation unpacking ===\n")
dmr <- annotate_probes(dm, anno)
cat(nrow(dmr), " gene x region rows from ", nrow(dm), " probes\n", sep = "")
# verify pairwise alignment survived: planted promoter probes must keep TSS labels
chk <- dmr[dmr$probe %in% probes[promoter_probes], ]
cat("planted promoter probes: ", nrow(chk), " rows, all TSS regions: ",
    all(chk$region %in% c("TSS200", "TSS1500")), "\n", sep = "")
cat("planted genes correctly assigned: ",
    all(chk$gene %in% true_genes), "\n", sep = "")
stopifnot(all(chk$region %in% c("TSS200", "TSS1500")), all(chk$gene %in% true_genes))

# multi-gene probe spot check
mg <- anno[grepl(";", anno$UCSC_RefGene_Name), ][1, ]
expanded <- dmr[dmr$probe == mg$probe, c("gene", "region")]
cat("multi-gene probe ", mg$probe, ": ", mg$UCSC_RefGene_Name, " / ",
    mg$UCSC_RefGene_Group, "\n  -> expanded to ", nrow(expanded), " rows\n", sep = "")
print(expanded, row.names = FALSE)

## --- handoff to pcsfomics ------------------------------------------------
cat("\n=== prize construction ===\n")
expr_s <- expression_scores(de_unpaired)
meth_s <- methylation_scores(dmr)
prizes <- build_prizes(expr_s, meth_s)

cat("expression terminals : ", nrow(expr_s), "\n", sep = "")
cat("methylation terminals: ", nrow(meth_s), "\n", sep = "")
cat("mean concordance -- planted: ",
    sprintf("%.3f", mean(prizes$concordance[prizes$gene %in% true_genes])),
    " | others: ",
    sprintf("%.3f", mean(prizes$concordance[!prizes$gene %in% true_genes])), "\n", sep = "")
cat("planted genes in top ", N_TRUE, " by prize: ",
    sum(head(prizes$gene, N_TRUE) %in% true_genes), "/", N_TRUE, "\n", sep = "")

cat("\nAll analysis-layer checks passed.\n")
