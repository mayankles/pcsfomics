# CLAUDE.md

Context for an LLM assistant picking up this project. Read this before making
changes.

## What this is

`pcsfomics` is an R package that integrates **differential expression** and
**differential methylation** at the **cohort level** over a signalling network,
using a **Prize-Collecting Steiner Forest** (PCSF) to pull out connected
subnetworks where both modalities agree.

The scientific motivation is mechanistic: promoter hypermethylation silencing a
transcription factor, which shifts its regulon, is a directional and testable
hypothesis. The method is built to surface exactly that shape of signal, and to
find **connector genes that carry no differential signal of their own** but are
structurally necessary — those are recovered as Steiner nodes and are invisible
to any threshold-and-intersect approach.

A second motivation shapes several design choices: the user wants to **compare
findings across datasets**, not just analyse one TCGA cohort. That is why
prizes are rank-based (see below).

The user prefers using R for omics data, but is also fluent in Python.

## Design decisions already made

These were argued through and settled. Each has a reason that is easy to miss.

1. **Prizes are rank-normalised, never `-log10(p)`.** P-values scale with sample
   size, so `-log10(p)` prizes from a 30-sample cohort and a 500-sample TCGA
   cohort are not on the same scale and their subnetworks cannot be meaningfully
   intersected. Rank-normalisation is what makes cross-dataset comparison valid.
   **If you change the prize function, preserve this property.**

2. **Promoter and gene-body methylation are scored separately.** Their expected
   correlation with expression has opposite sign — promoter methylation is
   typically silencing, gene-body methylation is often positively correlated
   with expression. Pooling them cancels signal.

3. **A concordance term is the actual integration mechanism.** Adding expression
   and methylation evidence together treats "strong but incoherent" the same as
   "strong and mechanistically consistent". `concordance_bonus()` gives full
   credit for promoter-hyper + down, partial credit for body-hyper + up, and
   zero otherwise, scaled by the geometric mean of the two scores. Benchmark:
   0.857 for coherent genes, 0.000 for genes with identical effect magnitudes in
   an incoherent direction.

4. **Shortest paths were considered and rejected.** PPI networks have diameter
   ~4-6, so nearly all gene pairs are 2-3 hops apart and path length is not
   discriminative; paths funnel through literature-biased hubs (UBC, ELAVL1,
   TP53); and argmax paths are unstable under small weight changes. PCSF
   aggregates over structure instead. Do not reintroduce a shortest-path scoring
   scheme.

5. **The network is OmniPath, not STRING or plain PPI.** A PPI edge means
   physical binding; methylation's effect on another gene's expression
   propagates through *regulation*. Adding the `collectri` dataset for TF-target
   edges is recommended and not yet done by default.

6. **PCSF is undirected.** Direction and sign from OmniPath ride along as edge
   attributes for post-hoc interpretation only. A recovered edge asserts
   connection, not causal flow. Say so in any methods text.

## Repository layout

```
R/prizes.R        rank_normalise, expression_scores, methylation_scores,
                  concordance_bonus, prize_config, build_prizes, prize_vector
R/network.R       fetch_omnipath, omnipath_graph, omnipath_df_to_graph,
                  load_omnipath, suggest_omega, component_summary,
                  OMNIPATH_FIELDS, OMNIPATH_DATASETS
R/solver.R        solve_pcsf, new_sp_cache, .kmb_steiner, .prune_tree_dp
R/robustness.R    robust_pcsf, degree_matched_shuffle, permutation_test
R/tcga_input.R    tcga_* barcode helpers, fetch_expression, fetch_methylation,
                  diff_expression, diff_methylation, annotate_probes,
                  filter_probes, summarise_design, ensembl_to_symbol

inst/examples/demo_synthetic.R       planted-module benchmark + beta sweep
inst/examples/test_tcga_analysis.R   analysis layer on simulated TCGA data
inst/examples/run_tcga_pipeline.R    end-to-end driver (needs GDC access)
inst/examples/get_omnipath.R         network download + troubleshooting
```

Only hard dependency is **igraph**. limma/edgeR are needed for `tcga_input.R`,
TCGAbiolinks only for the fetch layer.

## Verification status — read this before trusting anything

| Component | Status |
|---|---|
| `solve_pcsf` | Tested. Three unit cases (cheap path connects, costly path splits into 2 trees, low-prize leaf pruned), matching a Python reference implementation. |
| Prize construction, concordance | Tested on synthetic benchmark with planted module. |
| `robust_pcsf`, `permutation_test` | Run end-to-end on synthetic data. |
| `diff_expression`, `diff_methylation`, `annotate_probes`, barcode helpers | Tested via `test_tcga_analysis.R` on simulated TCGA-shaped data with planted signal. Recovers 39/40 planted genes. |
| `fetch_expression`, `fetch_methylation` | **NOT TESTED.** No GDC access in the dev environment. |
| `fetch_omnipath`, `omnipath_graph` | **Parsing and validation tested; the live HTTP round-trip is NOT.** omnipathdb.org was unreachable from the dev sandbox. |

The fetch layer is the untested surface. If something breaks, look there first.

## Benchmark findings that should guide tuning

From `demo_synthetic.R` (20 planted concordant genes in two network-distant
clusters, joined only through a signal-free connector, among 300 decoys):

```
definition                          nodes     core   prec   conn    p_obj
sum of all evidence (w_conc=1.5)      240    20/20   0.08  FALSE    0.484
concordance-weighted (w_conc=6)       240    20/20   0.08  FALSE    0.484
concordance ONLY                       22    19/20   0.86   TRUE    0.032
baseline (top-40 intersected)          17     9/20      -  FALSE       -
```

1. **Prize definition dominates solver tuning.** Sweeping `beta` 1→20 changed
   subnetwork size 240→449 but left precision flat (0.08→0.04). Changing the
   prize *definition* moved precision 0.08→0.86. Spend effort on the prize
   function, not on parameter sweeps.

2. **`p_objective` is the primary test statistic, not `p_compactness`.** Node
   count is weak because decoy terminals inflate observed and null subnetworks
   alike; it stayed non-significant in every regime. Objective value separated
   good from bad prize definitions cleanly (0.484 vs 0.032).

3. **Restrictive prizes did not cost discovery on this draw** — concordance-only
   recovered the signal-free connector *and* hit 0.86 precision. This is
   seed-dependent; on a different synthetic draw the permissive setting found
   the connector and the restrictive one did not. Check on real data rather than
   assuming either way.

## Gotchas already hit and fixed — do not regress these

- **`igraph::bfs()`'s `father` field does not survive `as_ids()`.** It returns
  vertex names in index order rather than each vertex's parent, and indexing it
  directly throws on the root's `NA`. `.prune_tree_dp()` uses a hand-rolled BFS
  instead. Commented in `solver.R`.
- **`%--%` only works inside igraph's indexing context.** Use
  `igraph::get.edge.ids()` for edge lookup.
- **`graph_from_adjacency_matrix(mode = "undirected")` warns** on floating-point
  distance matrices it cannot verify as exactly symmetric. Use `mode = "max"`.
- **OmniPath rejects unknown parameters with HTTP 200** and a one-line error page
  reading "Something is not entirely good.", which `read.delim` parses as a
  header. `fetch_omnipath()` validates against `OMNIPATH_FIELDS` /
  `OMNIPATH_DATASETS` before sending and detects the error page explicitly.
  **`n_references` and `n_resources` are NOT server fields** — OmnipathR derives
  them client-side, and so do we. Values are case-sensitive (`collectri`).
- **OmnipathR's `import_omnipath_interactions()` fails with
  `html_table` / `xml_missing`.** The cause is organism-name-to-taxid validation
  scraping the Ensembl species HTML page, not the interaction data. Passing a
  numeric taxon ID over plain HTTP avoids the code path. This has hit
  OmnipathR's own Bioconductor build server.
- **450K annotation fields are positionally aligned.** `UCSC_RefGene_Name` and
  `UCSC_RefGene_Group` pack multiple semicolon-delimited entries per probe, and
  entry *i* of one belongs with entry *i* of the other. Splitting them
  independently scrambles the gene-to-region mapping the promoter/body split
  depends on. `annotate_probes()` unpacks pairwise; the test asserts it.

## Known gaps — good next tasks

1. **Cross-reactive and SNP-overlapping probe lists are not shipped.**
   `filter_probes()` accepts `bad_probes` but defaults to empty. Chen et al.
   (2013) and Zhou et al. (2017) publish the standard exclusion lists. Skipping
   this inflates apparent signal at affected loci. **Highest-value gap.**
2. **CpG island context is not used.** `Relation_to_UCSC_CpG_Island`
   (island/shore/shelf/open sea) is a stronger regulatory stratifier than the
   promoter/body split alone. Island-shore methylation is the more reliable
   signal.
3. **Cis-correlation is assumed, not measured.** The concordance term assumes
   the canonical direction. Computing the actual per-gene methylation-expression
   correlation across patients would make seeds much stronger than merely
   appearing on two differential lists.
4. **Cross-dataset intersection is unimplemented.** This is the stated
   motivation for the whole project. The infrastructure is there (rank-based
   prizes, stable edge-frequency output); what is missing is a function to
   intersect or compare subnetworks across cohorts.
5. **No formal test suite.** `inst/examples/*.R` are scripts with `stopifnot()`,
   not `testthat` tests. Converting them would be straightforward.
6. **Cross-solver check not run.** `solve_pcsf(backend = "pcsf")` delegates to
   `PCSF::PCSF()` (Goemans-Williamson) versus the builtin KMB heuristic.
   Agreement between two approximation schemes on real data is worth reporting.

## How to run things

```r
# iterate without installing
for (f in list.files("R", full.names = TRUE)) source(f)

# the benchmark (a few minutes, no network)
Rscript inst/examples/demo_synthetic.R

# the analysis layer on simulated TCGA data (no network)
Rscript inst/examples/test_tcga_analysis.R

# network download
Rscript inst/examples/get_omnipath.R

# full pipeline (needs TCGAbiolinks + GDC access, downloads GB)
Rscript inst/examples/run_tcga_pipeline.R TCGA-LUAD results/
```

Cache directory is `.data/` (gitignored). Note the `fetch_*` functions default
to `cache_dir = "cache"` — pass `.data` explicitly or they will create a second
cache directory.

## Working style the user prefers

- Test claims rather than asserting them. When something is untested, say so.
- Flag when a previous conclusion turns out to be wrong; the user tracks these.
- Comments should explain *why*, especially where a choice looks arbitrary but
  is not (rank normalisation, M-values for testing, degree-matched permutation).
- The user is a bioinformatician and does not need statistics explained from
  first principles, but does want reasoning made explicit.
- Prefer working code over description of code.

## Analysis caveats worth remembering

- The user's own cohort is roughly 30 samples. **Tune `beta`, `omega`, `mu` on a
  large TCGA cohort, freeze them, then apply.** Parameter selection at n=30 will
  overfit undetectably.
- `summarise_design()` reports patients with *both* assays. For a cross-modal
  method that overlap is the real sample size, not the per-assay totals.
- TCGA adjacent normals are scarce (often 10-60 per tumour type). Paired designs
  gain power but restrict to complete pairs; the driver falls back to unpaired
  below 10 pairs.
- Check terminal coverage — what fraction of prized genes actually appear in the
  network. Ensembl-versus-symbol mismatch is the most common silent failure and
  looks exactly like a biological null result.
