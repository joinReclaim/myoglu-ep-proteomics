# =====================================================================
# R1_15_partial_phenotype.R — gjensidig justerte fenotypeassosiasjoner
# Adresserer R5 punkt 3 («the two outcomes are also not independent in
# this cohort ... a protein associating with one is expected to associate
# with the other») — direkte, i stedet for å omgå det.
#
# Spørsmålet er ikke hvilke proteiner som treffer BEGGE utfall (forventet
# når utfallene selv korrelerer), men hvilke som treffer ETT uavhengig av
# det andre. Rangbasert partiell korrelasjon.
# Krever: R1_01_longterm.R, R1_02_phenotypes.R
# =====================================================================
suppressPackageStartupMessages({library(data.table)})
source(file.path(Sys.getenv("ECV_ROOT", "/Users/sindrle/Research/Prosjekter/ECV"),
                 "review_round_1/01_scripts/R1_config.R"))
ROOT <- ECV_ROOT
OUT  <- file.path(ROOT,"review_round_1/02_output"); setwd(ROOT); set.seed(1)
PH <- readRDS(file.path(OUT,"R1_02_phenotypes.rds"))
P  <- readRDS(file.path(OUT,"R1_01_longterm.rds"))$PRIM
dph <- PH$dph
M <- PH$dp[intersect(rownames(PH$dp), unique(P[q < 0.25, Genes])), , drop = FALSE]

ok <- is.finite(dph$VO2max.kg) & is.finite(dph$M_mg.kg.min)
ct <- suppressWarnings(cor.test(dph$VO2max.kg[ok], dph$M_mg.kg.min[ok], method = "spearman"))
cat(sprintf("ΔVO₂max mot ΔM-value: rho = %+.3f (p = %.3g, n = %d)\n",
            ct$estimate, ct$p.value, sum(ok)))
cat(sprintf("proteiner i analysen (q<0.25): %d | deltakere: %d\n\n", nrow(M), ncol(M)))

# rangbasert partiell: y ~ x + z, p-verdi for x
part <- function(y, z) apply(M, 1, function(x) {
  k <- is.finite(x) & is.finite(y) & is.finite(z)
  if (sum(k) < 12) return(NA_real_)
  summary(lm(rank(y[k]) ~ rank(x[k]) + rank(z[k])))$coefficients[2, 4] })
marg <- function(y) apply(M, 1, function(x) {
  k <- is.finite(x) & is.finite(y); if (sum(k) < 10) return(NA_real_)
  suppressWarnings(cor.test(x[k], y[k], method = "spearman"))$p.value })

r <- data.table(Genes = rownames(M),
                p_VO2_marg = marg(dph$VO2max.kg),
                p_M_marg   = marg(dph$M_mg.kg.min),
                p_VO2_adj  = part(dph$VO2max.kg, dph$M_mg.kg.min),
                p_M_adj    = part(dph$M_mg.kg.min, dph$VO2max.kg))
exp_n <- 0.05 * nrow(M)
cat("=== NOMINELT SIGNIFIKANTE (forventet ved tilfeldighet: ", round(exp_n), ") ===\n", sep = "")
cat(sprintf("  ΔVO₂max, ujustert           : %3d\n", sum(r$p_VO2_marg < .05, na.rm = TRUE)))
cat(sprintf("  ΔVO₂max, justert for ΔM     : %3d  <- signal på kondisjonsaksen\n", sum(r$p_VO2_adj < .05, na.rm = TRUE)))
cat(sprintf("  ΔM-value, ujustert          : %3d\n", sum(r$p_M_marg < .05, na.rm = TRUE)))
cat(sprintf("  ΔM-value, justert for ΔVO₂max: %3d  <- intet signal på metabolsk akse\n", sum(r$p_M_adj < .05, na.rm = TRUE)))

cat("\n=== KANDIDATENE ===\n")
print(r[Genes %in% c("LTBP4","VCL","MEP1B","HSP90B1","TPI1"),
        .(Genes, p_VO2_marg = signif(p_VO2_marg,3), p_VO2_adj = signif(p_VO2_adj,3),
          p_M_marg = signif(p_M_marg,3), p_M_adj = signif(p_M_adj,3))])
fwrite(r[order(p_VO2_adj)], file.path(OUT, "R1_15_partial_phenotype.csv"))
cat("\nskrevet: R1_15_partial_phenotype.csv\n")
