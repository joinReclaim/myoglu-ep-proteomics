# =====================================================================
# R1_24_fig5AB.R — manuskriptets FIGUR 5, PANEL A (Figures/Figure3.tif)
#
# KOMPOSITTSCOREN ER TATT UT AV FIGUREN. Den fantes for å slippe fire XY-plott,
# men lollipop-panelet viser de samme fire korrelasjonene direkte og er den
# ærlige versjonen. Permutasjonsnullens median er +0.774 mot observert +0.724
# (R1_02), altså under forventning — panelet ville hatt en figurtekst som sier
# at det ikke viser noe. Beholdt i koden for sporbarhet, skrives til fil merket
# NOT_USED. Lollipop blir nytt panel A.
# Erstatter 005_ECV_phenotype_correlations.R:251-316
#
# Uendret: A = XY med geom_smooth(lm, se=FALSE, black, size 1.2), farger
# dysglyc=red / control=black, theme_classic, 8 x 8 cm.
# B = lollipop (geom_segment + geom_point size 4), coord_flip, 8 x 8 cm.
# Endret: loess-normaliserte data, ingen mean-imputering, ECV -> EP i
# aksetitler (jf. Q- to be fixed).
# MERK: kompositten i A er en ILLUSTRASJON for å slippe fire XY-plott.
# Permutasjonstest (R1_02): nullfordelingens median rho = +0.774, publisert
# kompositt +0.724 -> p = 0.77. Dette skal stå i legenden.
# Krever: R1_02_phenotypes.R
# =====================================================================
suppressPackageStartupMessages({library(data.table)})
source(file.path(Sys.getenv("ECV_ROOT", "/Users/sindrle/Research/Prosjekter/ECV"),
                 "review_round_1/01_scripts/R1_config.R"))
ROOT <- ECV_ROOT
OUT  <- file.path(ROOT,"review_round_1/02_output"); setwd(ROOT)
source(file.path(ROOT,"review_round_1/01_scripts/R1_theme.R"))
PH <- readRDS(file.path(OUT,"R1_02_phenotypes.rds"))
CAND <- c("LTBP4","VCL","MEP1B","HSP90B1")

M  <- PH$dp[intersect(CAND, rownames(PH$dp)), , drop = FALSE]
z  <- t(scale(t(M)))
sc <- colMeans(z, na.rm = TRUE)
d  <- data.table(score = sc, dVO2 = PH$dph$VO2max.kg, Group = PH$dph$Gr)
d  <- d[is.finite(score) & is.finite(dVO2)]
ct <- suppressWarnings(cor.test(d$score, d$dVO2, method = "spearman"))
cat(sprintf("panel A: n=%d | rho = %+.3f, p = %.4f\n", nrow(d), ct$estimate, ct$p.value))

pA <- ggplot(d, aes(x = score, y = dVO2)) +
  geom_point(aes(colour = Group), alpha = 1) +
  geom_smooth(aes(group = 1), method = "lm", se = FALSE, colour = "black", linewidth = 1.2) +
  scale_colour_manual(values = PAL$group_xy) +
  labs(x = expression(paste(Delta, "EP VO"[2], "max score, log"[2], "(int)")),
       y = expression(paste(Delta, "VO"[2], "max (mL/kg/min)"))) +
  theme_ecv() + theme(legend.position = "none")
save_ecv("R1_24_fig5_composite_NOT_USED.pdf", pA, width = 8, height = 8, unit = "cm")

rho <- data.table(protein = rownames(M), rho = apply(M, 1, function(x) {
  ok <- is.finite(x) & is.finite(PH$dph$VO2max.kg)
  suppressWarnings(cor(x[ok], PH$dph$VO2max.kg[ok], method = "spearman")) }))
cat("panel B:\n"); print(rho[order(-rho)])
pB <- ggplot(rho, aes(x = reorder(protein, rho), y = rho)) +
  geom_segment(aes(xend = protein, y = 0, yend = rho)) +
  geom_point(size = 4) + coord_flip() +
  labs(x = NULL, y = expression(paste("Spearman's rho with ", Delta, "VO"[2], "max"))) +
  theme_ecv() + theme(legend.position = "none",
                      plot.margin = margin(t = 10, r = 30, b = 10, l = 10))
save_ecv("R1_24_fig5A_lollipop.pdf", pB, width = 8, height = 8, unit = "cm")
fwrite(rho[order(-rho)], file.path(OUT, "R1_24_fig5A_rho.csv"))
