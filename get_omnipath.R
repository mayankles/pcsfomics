## get_omnipath.R -- work around the OmnipathR Ensembl-scraping failure.
##
## The error
##
##   no applicable method for 'html_table' applied to an object of class "xml_missing"
##
## is NOT a problem with the interaction data or with your setup. Look at where
## it happens in the traceback: omnipath_check_param -> ncbi_taxid ->
## get_db("organisms") -> ensembl_organisms_raw() -> html_table().
##
## OmnipathR validates the `organisms` parameter before sending any query, and
## that validation translates organism names to NCBI taxon IDs by SCRAPING
## https://www.ensembl.org/info/about/species.html with rvest. When that page's
## structure changes, html_element() returns xml_missing and html_table() dies.
## The OmniPath web service itself is unaffected -- the query never got sent.
##
## This has hit OmnipathR's own Bioconductor build server with the identical
## traceback, so it is a known upstream fragility, not something you configured
## wrong.
##
## Options below are ordered by how likely they are to just work.

source("R/network.R")
source("R/prizes.R")          # rank_normalise
suppressMessages(library(igraph))


## ---------------------------------------------------------------------------
## Option 1: skip OmnipathR entirely (recommended)
## ---------------------------------------------------------------------------
## The web service is a plain HTTP endpoint returning TSV. Passing the taxon ID
## as a NUMBER means nothing ever needs to look up an organism name.

g <- omnipath_graph(datasets = "omnipath", organism = 9606L,
                    destfile = "omnipath.tsv")

cat(sprintf("\nnetwork: %d nodes, %d edges\n", vcount(g), ecount(g)))
cat("saved to omnipath.tsv -- load_omnipath('omnipath.tsv') will read it back\n")

## Add TF-target edges. Worth doing for this project: methylation acts on
## transcription, and a promoter-hypermethylation-to-expression story propagates
## through regulatory edges, not through physical binding.
##
##   g <- omnipath_graph(datasets = c("omnipath", "collectri"), organism = 9606L)


## ---------------------------------------------------------------------------
## Option 2: clear the OmnipathR cache, then retry
## ---------------------------------------------------------------------------
## OmnipathR caches downloads including the Ensembl page. If a truncated or
## error page got cached, every retry replays the same failure. Wiping forces a
## fresh fetch and sometimes that is all it takes.
##
##   library(OmnipathR)
##   omnipath_cache_wipe()
##   # or inspect first:
##   getOption("omnipathr.cachedir")
##   omnipath_cache_search("ensembl")
##
## Then retry with a NUMERIC organism, which is worth doing regardless:
##
##   ia <- import_omnipath_interactions(organism = 9606L)


## ---------------------------------------------------------------------------
## Option 3: update OmnipathR
## ---------------------------------------------------------------------------
## Fixes for the scraping code usually land in the development version before
## they reach the release branch.
##
##   BiocManager::install("saezlab/OmnipathR")
##
## Enable logging to see the actual HTTP exchange if it still fails:
##
##   options(omnipathr.console_loglevel = "trace")


## ---------------------------------------------------------------------------
## Sanity check whatever you loaded
## ---------------------------------------------------------------------------
## Before spending a run on it, confirm the network looks right. A few hub genes
## should be present and well connected; if TP53 has three neighbours, something
## is wrong with the identifiers.

for (hub in c("TP53", "EGFR", "MYC", "CTNNB1")) {
  if (hub %in% V(g)$name) {
    cat(sprintf("%-8s degree %d\n", hub, degree(g, hub)))
  } else {
    cat(sprintf("%-8s ABSENT -- check gene symbol conventions\n", hub))
  }
}

comp <- components(g)
cat(sprintf("\nlargest connected component: %d of %d nodes (%.0f%%)\n",
            max(comp$csize), vcount(g), 100 * max(comp$csize) / vcount(g)))
cat("edge cost range:", sprintf("%.3f", range(E(g)$cost)), "\n")