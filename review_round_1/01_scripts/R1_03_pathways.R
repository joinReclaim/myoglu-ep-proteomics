# =====================================================================
# R1_03_pathways.R — ORA mot MÅLT proteom, separat opp/ned
# Adresserer R5 punkt 4, R2 concern 4
#
# Feilen i 006_ECV_pathway_analysis.R: bg_ids ble beregnet (linje 277-282)
# men aldri sendt som universe= til enrichGO/enrichKEGG, så all enrichment
# gikk mot helgenomet (BgRatio-nevner 18986). I 101_ECV_acute_chronic_
# overlap.R er universe eksplisitt utkommentert.
#
# Krever: R1_01_longterm.R
# =====================================================================
suppressPackageStartupMessages({
  library(data.table); library(clusterProfiler); library(org.Hs.eg.db)
})
source(file.path(Sys.getenv("ECV_ROOT", "/Users/sindrle/Research/Prosjekter/ECV"),
                 "review_round_1/01_scripts/R1_config.R"))
ROOT <- ECV_ROOT
OUT  <- file.path(ROOT, "review_round_1/02_output"); setwd(ROOT)
P <- readRDS(file.path(OUT, "R1_01_longterm.rds"))$PRIM

univ <- unique(P$Genes)                    # ALT som faktisk ble testet
sig  <- P[q < 0.05]
sets <- list(alle = unique(sig$Genes),
             opp  = unique(sig[Estimate > 0]$Genes),
             ned  = unique(sig[Estimate < 0]$Genes))
cat(sprintf("universe (målt proteom): %d | sig %d (opp %d / ned %d)\n",
            length(univ), length(sets$alle), length(sets$opp), length(sets$ned)))

to_entrez <- function(x) unique(na.omit(suppressMessages(
  bitr(x, "SYMBOL", "ENTREZID", org.Hs.eg.db))$ENTREZID))
U <- to_entrez(univ)
cat(sprintf("universe mappet til Entrez: %d\n\n", length(U)))

run_one <- function(genes, lab, ont) {
  G <- to_entrez(genes); if (!length(G)) return(NULL)
  e <- tryCatch({
    if (ont == "KEGG") enrichKEGG(G, universe = U, organism = "hsa",
                                  pAdjustMethod = "BH", pvalueCutoff = .05)
    else enrichGO(G, OrgDb = org.Hs.eg.db, keyType = "ENTREZID", ont = ont,
                  universe = U, pAdjustMethod = "BH", pvalueCutoff = .05,
                  readable = TRUE)
  }, error = function(e) NULL)
  if (is.null(e) || !nrow(as.data.frame(e))) return(NULL)
  d <- as.data.table(as.data.frame(e))[, `:=`(set = lab, ontology = ont)]
  d
}
all <- rbindlist(lapply(names(sets), function(s)
  rbindlist(lapply(c("BP","CC","KEGG"), function(o) run_one(sets[[s]], s, o)),
            fill = TRUE)), fill = TRUE)

if (!nrow(all)) {
  cat("INGEN signifikant enrichment mot målt proteom.\n")
} else {
  fwrite(all, file.path(OUT, "R1_03_pathways_measured_universe.csv"))
  for (s in names(sets)) for (o in c("BP","CC","KEGG")) {
    d <- all[set == s & ontology == o]
    if (!nrow(d)) { cat(sprintf("--- %s / %s: ingen ---\n", s, o)); next }
    cat(sprintf("\n--- %s / %s (%d termer) ---\n", s, o, nrow(d)))
    print(head(d[order(p.adjust)][, .(Description, GeneRatio, BgRatio,
                                      q = signif(p.adjust,3), Count)], 8))
  }
}
cat("\n=== Til sammenligning: HELGENOM-bakgrunn (som i manuskriptet) ===\n")
Uold <- U; U <<- NULL
wg <- rbindlist(lapply(c("BP","KEGG"), function(o) {
  G <- to_entrez(sets$alle)
  e <- tryCatch(if (o=="KEGG") enrichKEGG(G, organism="hsa", pvalueCutoff=.05)
                else enrichGO(G, OrgDb=org.Hs.eg.db, keyType="ENTREZID",
                              ont=o, pvalueCutoff=.05, readable=TRUE),
                error=function(e) NULL)
  if (is.null(e)||!nrow(as.data.frame(e))) return(NULL)
  as.data.table(as.data.frame(e))[, ontology := o]
}), fill=TRUE)
if (nrow(wg)) print(head(wg[order(p.adjust)][, .(ontology, Description, BgRatio,
                                                 q=signif(p.adjust,3), Count)], 10))
