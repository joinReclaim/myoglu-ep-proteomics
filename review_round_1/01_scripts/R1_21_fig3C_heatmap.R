# =====================================================================
# R1_21_fig3C_heatmap.R — manuskriptets FIGUR 3, PANEL C (Figures/Figure1.tif)
# Erstatter 006_ECV_pathway_analysis.R:242-256
#
# Uendret: pheatmap, treeheight 0, annotasjonsfelt, ingen rad-/kolonnenavn,
# colorRampPalette(c("blue","white","red"))(100), breaks -1..1, 5 x 4 tommer,
# bg = "transparent".
# Endret: ny responderliste og loess-data. Genererer BÅDE sentrert og
# usentrert (R3 ber eksplisitt om å se begge).
# Krever: R1_01_longterm.R, R1_02_phenotypes.R
# =====================================================================
suppressPackageStartupMessages({library(data.table); library(pheatmap); library(cluster)})
source(file.path(Sys.getenv("ECV_ROOT", "/Users/sindrle/Research/Prosjekter/ECV"),
                 "review_round_1/01_scripts/R1_config.R"))
ROOT <- ECV_ROOT
OUT  <- file.path(ROOT,"review_round_1/02_output"); setwd(ROOT)
source(file.path(ROOT,"review_round_1/01_scripts/R1_theme.R"))
P  <- readRDS(file.path(OUT,"R1_01_longterm.rds"))$PRIM
PH <- readRDS(file.path(OUT,"R1_02_phenotypes.rds"))

THR <- 0.25    # sett 0.05 for primærlisten
K   <- 2       # silhouette-maks; originalen brukte k=7 uten begrunnelse

sig <- intersect(rownames(PH$dp), unique(P[q < THR, Genes]))
D <- PH$dp[sig, , drop = FALSE]; D <- D[rowSums(is.finite(D)) >= 20, , drop = FALSE]
Dc <- sweep(D, 2, colMeans(D, na.rm = TRUE), "-")

mk <- function(M, tag) {
  cm <- cor(t(M), use = "pairwise.complete.obs", method = "spearman")
  cm[!is.finite(cm)] <- 0; diag(cm) <- 1
  hc <- hclust(as.dist(1 - cm), method = "average")
  cl <- cutree(hc, K)
  sil <- mean(silhouette(cl, as.dist(1 - cm))[, 3])
  cat(sprintf("%-11s n=%d | silhouette k=%d: %.3f | blokker: %s\n",
              tag, nrow(M), K, sil, paste(as.integer(table(cl)), collapse = "/")))
  ann <- data.frame(Module = factor(paste0("M", cl))); rownames(ann) <- rownames(cm)
  f <- file.path(FIGDIR, sprintf("R1_21_fig3C_%s.pdf", tag))
  cairo_pdf(f, width = 5, height = 4, bg = "transparent")
  pheatmap(cm, cluster_rows = hc, cluster_cols = hc, treeheight_row = 0,
           treeheight_col = 0, annotation_row = ann, annotation_col = ann,
           show_rownames = FALSE, show_colnames = FALSE,
           color = PAL$heat, breaks = seq(-1, 1, length.out = 101))
  dev.off(); message("figur: ", f)
  invisible(list(cl = cl, cm = cm))
}
a <- mk(Dc, "centred")
b <- mk(D,  "uncentred")
cat(sprintf("\nmoduletilhørighet lik i begge versjoner: %.0f%% av proteinene\n",
            100 * max(mean(a$cl == b$cl), mean(a$cl != b$cl))))
fwrite(data.table(Genes = names(a$cl), module_centred = paste0("M", a$cl),
                  module_uncentred = paste0("M", b$cl)),
       file.path(OUT, "R1_21_fig3C_modules.csv"))
