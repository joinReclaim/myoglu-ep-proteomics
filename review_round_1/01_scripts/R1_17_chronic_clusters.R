# =====================================================================
# R1_17_chronic_clusters.R — korrelasjonsklynging av responderne (Fig 3C/3D)
# Erstatter 006_ECV_pathway_analysis.R linje 192-260
#
# Endringer:
#   - loess-normaliserte verdier og ny responderliste (ikke unormalisert/314)
#   - klyngingen vises MED OG UTEN deltakernivå-sentrering (R3 ber om dette)
#   - den induserte korrelasjonen av sentreringen kvantifiseres (R5 minor:
#     å trekke fra et deltakergjennomsnitt beregnet fra de samme proteinene
#     induserer negativ korrelasjon av orden 1/p — påstanden om at det
#     «cannot induce between-protein correlations» er ikke riktig)
#   - ORA per klynge mot MÅLT proteom (ikke helgenom)
# Krever: R1_01_longterm.R, R1_02_phenotypes.R
# =====================================================================
suppressPackageStartupMessages({library(data.table); library(pheatmap)
  library(cluster); library(clusterProfiler); library(org.Hs.eg.db)})
source(file.path(Sys.getenv("ECV_ROOT", "/Users/sindrle/Research/Prosjekter/ECV"),
                 "review_round_1/01_scripts/R1_config.R"))
ROOT <- ECV_ROOT
OUT  <- file.path(ROOT,"review_round_1/02_output"); setwd(ROOT); set.seed(1)
source(file.path(ROOT,"review_round_1/01_scripts/R1_theme.R"))
P  <- readRDS(file.path(OUT,"R1_01_longterm.rds"))$PRIM
PH <- readRDS(file.path(OUT,"R1_02_phenotypes.rds"))

for (THR in c(0.05, 0.25)) {
  sig <- intersect(rownames(PH$dp), unique(P[q < THR, Genes]))
  D <- PH$dp[sig, , drop = FALSE]                     # proteiner x deltakere
  D <- D[rowSums(is.finite(D)) >= 20, , drop = FALSE]
  cat(sprintf("\n########## RESPONDERE q<%.2f: %d proteiner x %d deltakere ##########\n",
              THR, nrow(D), ncol(D)))

  # deltakernivå-sentrering: trekk fra hver deltakers gjennomsnittsdelta
  gm <- colMeans(D, na.rm = TRUE)
  Dc <- sweep(D, 2, gm, "-")

  cm_raw <- cor(t(D),  use = "pairwise.complete.obs", method = "spearman")
  cm_ctr <- cor(t(Dc), use = "pairwise.complete.obs", method = "spearman")
  off <- function(m) mean(m[upper.tri(m)], na.rm = TRUE)
  cat(sprintf("  gj.snittlig korrelasjon mellom proteiner: usentrert %+.4f | sentrert %+.4f\n",
              off(cm_raw), off(cm_ctr)))
  cat(sprintf("  forventet indusert skift ved sentrering: -1/p = %+.4f (observert %+.4f)\n",
              -1/nrow(D), off(cm_ctr) - off(cm_raw)))
  cat("  => sentreringen INDUSERER en liten negativ korrelasjon. Neglisjerbar i\n")
  cat("     størrelse, men påstanden om at den ikke kan gjøre det er ikke riktig.\n")

  for (nm in c("usentrert","sentrert")) {
    cm <- if (nm == "usentrert") cm_raw else cm_ctr
    cm[!is.finite(cm)] <- 0; diag(cm) <- 1
    hc <- hclust(as.dist(1 - cm), method = "average")
    sil <- sapply(2:8, function(k) mean(silhouette(cutree(hc, k), as.dist(1-cm))[,3]))
    cat(sprintf("\n  --- %s --- silhouette k=2..8: %s\n", nm,
                paste(sprintf("%.2f", sil), collapse = " ")))
    K <- (2:8)[which.max(sil)]
    cl <- cutree(hc, K)
    cat(sprintf("      beste k=%d, størrelser: %s\n", K,
                paste(as.integer(table(cl)), collapse = ", ")))
    if (THR == 0.25 && nm == "sentrert") {
      ann <- data.frame(Klynge = factor(paste0("C", cl))); rownames(ann) <- rownames(cm)
      f <- file.path(FIGDIR, "R1_17_chronic_cluster_heatmap.pdf")
      cairo_pdf(f, width = 5, height = 4, bg = "transparent")
      pheatmap(cm, cluster_rows = hc, cluster_cols = hc, treeheight_row = 0,
               treeheight_col = 0, annotation_row = ann, annotation_col = ann,
               show_rownames = FALSE, show_colnames = FALSE,
               color = colorRampPalette(c("blue","white","red"))(100),
               breaks = seq(-1, 1, length.out = 101))
      dev.off(); message("figur: ", f)
      fwrite(data.table(Genes = names(cl), cluster = paste0("C", cl)),
             file.path(OUT, "R1_17_chronic_clusters.csv"))

      cat("\n  --- ORA per klynge, universe = MÅLT proteom ---\n")
      U <- unique(na.omit(suppressMessages(bitr(unique(P$Genes), "SYMBOL",
                          "ENTREZID", org.Hs.eg.db))$ENTREZID))
      for (k in sort(unique(cl))) {
        g <- names(cl)[cl == k]
        G <- unique(na.omit(suppressMessages(bitr(g, "SYMBOL", "ENTREZID",
                            org.Hs.eg.db))$ENTREZID))
        e <- tryCatch(enrichGO(G, OrgDb = org.Hs.eg.db, keyType = "ENTREZID",
                ont = "BP", universe = U, pAdjustMethod = "BH",
                pvalueCutoff = .05, readable = TRUE), error = function(e) NULL)
        d <- if (is.null(e)) NULL else as.data.table(as.data.frame(e))
        if (is.null(d) || !nrow(d)) { cat(sprintf("    C%d (n=%d): ingen signifikante termer\n", k, length(g))); next }
        cat(sprintf("    C%d (n=%d): %s\n", k, length(g),
            paste(head(d[order(p.adjust)]$Description, 3), collapse = " | ")))
        fwrite(d[, cluster := paste0("C",k)],
               file.path(OUT, sprintf("R1_17_ORA_C%d.csv", k))) }
    }
  }
}
