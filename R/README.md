# pcsfomics

Cohort-level Prize-Collecting Steiner Forest over a signalling network, with
prizes built from differential expression + differential methylation.

Tested on R 4.3.3 / igraph 1.6.0. Only hard dependency is **igraph**.

## Install

```r
# from the package directory
install.packages(".", repos = NULL, type = "source")

# or just source the files while iterating
for (f in list.files("R", full.names = TRUE)) source(f)
```

## Usage

```r
library(igraph)

g    <- load_omnipath("omnipath.tsv")
expr <- expression_scores(de_table)     # gene, log2FC, padj
meth <- methylation_scores(dmr_table)   # gene, region, delta_beta, padj

prizes <- prize_vector(build_prizes(expr, meth, prize_config(beta = 1)))
omega  <- suggest_omega(prizes, 0.90)

res <- solve_pcsf(g, prizes, omega = omega, mu = 0.005)
component_summary(res$subnetwork)

freq <- robust_pcsf(g, prizes, omega = omega, mu = 0.005, n_runs = 100)
perm <- permutation_test(g, prizes, omega = omega, mu = 0.005, n_perm = 200)
```

Run the benchmark:

```r
Rscript inst/examples/demo_synthetic.R
```

## Getting the network

```r
library(OmnipathR)
ia <- import_omnipath_interactions()
write.table(ia, "omnipath.tsv", sep = "\t", row.names = FALSE, quote = FALSE)
```

`load_omnipath()` derives edge cost from `curation_effort` and carries
`is_directed` / `is_stimulation` / `is_inhibition` through as edge attributes.
PCSF itself is undirected — a recovered edge asserts connection, not causal
flow. Use the sign attributes for post-hoc interpretation and say so in methods.

## Two solver backends

```r
solve_pcsf(g, prizes, backend = "builtin")  # default, pure R + igraph
solve_pcsf(g, prizes, backend = "pcsf")     # delegates to PCSF::PCSF()
```

`builtin` is a KMB Steiner heuristic followed by an exact bottom-up tree DP for
prize pruning. `pcsf` calls the Goemans–Williamson heuristic of Akhmedov et al.
(`remotes::install_github("murodzhon/PCSF")`; needs Boost, GitHub-only, docs
last built 2019).

Running both on identical inputs is a real robustness check. Two different
approximation schemes converging on the same core is worth a methods sentence;
sharp disagreement means your parameters sit in a regime where the solution
isn't well determined.

## Parameters

| Parameter | Meaning | Starting point |
|---|---|---|
| `beta` | prize scale vs. edge cost | 1; sweep 1–20 |
| `omega` | per-tree admission charge | `suggest_omega(prizes, 0.90)` |
| `mu` | hub penalty (`prize - mu * degree`) | 0.005 |

`mu` between 1e-4 and 5e-2 is the range the PCSF package authors recommend for
biological networks, selected by inspecting the resulting Steiner/terminal
ratio. `omega` only does useful work inside the prize distribution — far below
it every terminal is retained, far above it nothing is admitted.

## What the benchmark shows

Planted module: 20 concordant genes in two network-distant clusters, joined only
through a signal-free connector, among 300 decoys carrying single-modality or
direction-incoherent signal.

```
concordance separates coherent from incoherent genes:
  core (hyper + down)        mean concordance = 0.857
  discordant (hypo + down)   mean concordance = 0.000

prize definition comparison (beta held at 1):
definition                          nodes     core   prec   conn    p_obj
sum of all evidence (w_conc=1.5)      240    20/20   0.08  FALSE    0.484
concordance-weighted (w_conc=6)       240    20/20   0.08  FALSE    0.484
concordance ONLY                       22    19/20   0.86   TRUE    0.032

baseline (top-40 each, intersected):  17 genes, core 9/20, connector FALSE
```

1. **The concordance term works.** Genes with promoter hypermethylation and
   downregulation score 0.857; genes with the same effect magnitudes in an
   incoherent direction score 0.000. Summing the two modalities cannot make
   that distinction.

2. **Prize definition dominates solver tuning.** Sweeping `beta` from 1 to 20
   changed subnetwork size from 240 to 449 nodes and left precision flat
   (0.08 → 0.04). Changing the prize *definition* moved precision from 0.08 to
   0.86 and flipped the permutation test from p = 0.484 to p = 0.032. Spend
   your time on the prize function.

3. **Restrictive prizes did not cost discovery here.** The concordance-only
   definition recovered the signal-free connector as a Steiner node while the
   permissive definition missed it inside a 240-node blob. Fewer, better
   terminals gave the solver a cleaner problem. (This is draw-dependent — on a
   different synthetic seed the permissive setting found the connector and the
   restrictive one did not, so check it on your own data rather than assuming.)

4. **Node count is a weak test statistic; objective is better.** Decoy terminals
   inflate observed and null subnetworks alike, so `p_compactness` stayed
   non-significant in every regime. `p_objective` separated the good prize
   definition from the bad ones cleanly. `permutation_test()` returns both.

## Getting the inputs from TCGA

`R/tcga_input.R` covers download through to `de_table.tsv` / `dmr_table.tsv`.

```bash
Rscript inst/examples/run_tcga_pipeline.R TCGA-LUAD results/
```

Needs TCGAbiolinks and GDC access. Downloads are cached as `.rds`, so re-runs
skip the network. The fetch layer (`fetch_*`) and the analysis layer (`diff_*`,
`annotate_probes`) are kept separate so you can iterate on the analysis without
re-downloading.

Verify the analysis layer without touching the network:

```bash
Rscript inst/examples/test_tcga_analysis.R
```

That runs the whole chain on simulated TCGA-shaped data with a planted signal
(40 genes downregulated + promoter hypermethylated) and asserts recovery.

### Things this handles that are easy to get wrong

- **Aliquot duplicates.** TCGA often has several aliquots per patient and sample
  type. Leaving them in inflates n with technical replicates. `dedup_aliquots()`
  keeps one, deterministically.
- **Semicolon-delimited annotation is positionally aligned.**
  `UCSC_RefGene_Name` and `UCSC_RefGene_Group` pack multiple entries per probe,
  and entry *i* of one belongs with entry *i* of the other. Splitting them
  independently scrambles exactly the gene-to-region mapping the promoter/body
  split depends on. `annotate_probes()` unpacks them pairwise and the test
  asserts it.
- **M-values for testing, beta for effect size.** Beta values are bounded and
  heteroscedastic at the extremes, which breaks the linear model's
  constant-variance assumption (Du et al. 2010). Effect sizes are reported back
  on the beta scale so "30% methylation difference" stays interpretable.
- **Ensembl-to-symbol mapping.** STAR counts use versioned Ensembl IDs; OmniPath
  uses symbols. The driver reports what fraction of prized genes actually appear
  in the network and warns below 50%. This is the most common silent failure in
  network integration, and it looks like a biological null result.
- **Paired vs unpaired.** Blocking on patient is a real power gain, but TCGA
  normals are scarce and blocking restricts you to complete pairs.
  `summarise_design()` reports the pair counts before you commit; the driver
  falls back to unpaired below 10 pairs.

### What you still need to add

Cross-reactive and SNP-overlapping probe exclusion lists (Chen et al. 2013,
Zhou et al. 2017). `filter_probes()` takes them via `bad_probes` but ships
empty — the lists are external files. Skipping this inflates apparent signal at
affected loci.

## Notes

- Prizes are rank-based, not `-log10(p)`-based, specifically so subnetworks from
  cohorts of different sizes and platforms stay comparable. Preserve that
  property if you modify the prize function — it is what makes cross-dataset
  intersection meaningful.
- `robust_pcsf()` perturbs **prizes**; `PCSF::PCSF_rand()` perturbs **edge
  costs**. Prizes from a ~30-sample cohort are far less stable than OmniPath
  curation scores, so prize noise is usually the more honest perturbation.
  Reporting edges stable under both is stronger.
- Shortest-path distances are cached in `new_sp_cache()` and reused across
  resampling replicates. Pass one explicitly if you write your own loop.
- Tune `beta` and `omega` on a large TCGA cohort, freeze them, then apply to the
  small cohort. Parameter selection on ~30 samples will overfit.
