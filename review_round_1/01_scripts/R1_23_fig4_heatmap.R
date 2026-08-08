# =====================================================================
# R1_23_fig4_heatmap.R — manuskriptets FIGUR 4 (Figures/Figure2.tif)
# Erstatter 005_ECV_phenotype_correlations.R:185-240
#
# Uendret fra originalen: pheatmap, clustering på rader og kolonner,
# show_colnames = F, labels_row = nice_names, blue-white-red,
# fontsize_row 10, legend_breaks -0.6/0/0.6, 8 x 4 TOMMER, bg transparent.
# Originalen satte IKKE breaks — fargeskalaen auto-skalerer til dataområdet.
# Beholdt slik.
# Endret: loess-normaliserte verdier, ny responderliste (q<0.25), INGEN
# mean-imputering av proteinverdier (005:178-182 gjorde det).
# Boks og pil rundt VO2max-raden settes i PowerPoint, som i originalen.
# Krever: R1_01_longterm.R, R1_02_phenotypes.R
# =====================================================================
suppressPackageStartupMessages({library(data.table); library(pheatmap)})
source(file.path(Sys.getenv("ECV_ROOT", "/Users/sindrle/Research/Prosjekter/ECV"),
                 "review_round_1/01_scripts/R1_config.R"))
ROOT <- ECV_ROOT
OUT  <- file.path(ROOT,"review_round_1/02_output"); setwd(ROOT)
source(file.path(ROOT,"review_round_1/01_scripts/R1_theme.R"))
PH <- readRDS(file.path(OUT,"R1_02_phenotypes.rds"))
P  <- readRDS(file.path(OUT,"R1_01_longterm.rds"))$PRIM
THR <- 0.25

nice_names <- c(
  Weight = "ΔBody weight (kg)", BMI = "ΔBMI (kg/m²)",
  FFM = "ΔFat-free mass (kg)", M_mg.kg.min = "ΔM-value (mg/kg/min)",
  VO2max.kg = "ΔVO₂max (mL/kg/min)", ChestPress = "ΔChest press (kg)",
  PullDown = "ΔPull-down (kg)", LegPress = "ΔLeg press (kg)",
  mri_VL = "ΔVastus lateralis volume", mri_IP = "ΔIntraperitoneal fat (AU)",
  mri_TAT = "ΔTotal adipose tissue (L)", mri_SAT = "ΔSubcutaneous adipose tissue (L)",
  mrs_Hep = "ΔLiver fat (AU)", mrs_Musc = "ΔMuscle fat (AU)")

sig  <- intersect(rownames(PH$dp), unique(P[q < THR, Genes]))
M    <- PH$dp[sig, , drop = FALSE]
cm <- sapply(PH$PHENOS, function(v)
  apply(M, 1, function(x) { ok <- is.finite(x) & is.finite(PH$dph[[v]])
    if (sum(ok) < 10) NA_real_ else
      suppressWarnings(cor(x[ok], PH$dph[[v]][ok], method = "spearman")) }))
cm <- t(cm)
cm <- cm[, colSums(is.na(cm)) == 0, drop = FALSE]
cat(sprintf("q<%.2f: %d respondere -> %d proteiner med komplett korrelasjon x %d fenotyper\n",
            THR, length(sig), ncol(cm), nrow(cm)))
cat(sprintf("rho-område: %.2f til %.2f\n", min(cm), max(cm)))

f <- file.path(FIGDIR, "R1_23_fig4_phenotype_heatmap.pdf")
cairo_pdf(f, width = 8, height = 4, bg = "transparent")
pheatmap(cm,
         cluster_rows = TRUE, cluster_cols = TRUE, show_colnames = FALSE,
         labels_row = nice_names[rownames(cm)],
         color = PAL$heat,
         fontsize_row = 10, fontsize_col = 6, legend = TRUE,
         legend_breaks = c(-0.6, 0, 0.6), legend_labels = c("-0.6", "0", "0.6"))
dev.off(); message("figur: ", f)
fwrite(data.table(phenotype = rownames(cm), as.data.table(cm)),
       file.path(OUT, "R1_23_fig4_matrix.csv"))
