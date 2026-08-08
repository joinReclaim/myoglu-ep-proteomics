# =====================================================================
# R1_02_phenotypes.R — fenotypeassosiasjoner med multiplisitetskontroll
# Adresserer R5 punkt 3, R1 punkt 4, R3 punkt 2, R4 (multiple testing)
#
#   FDR over HELE protein x fenotype-matrisen | antall fenotyper oppgitt
#   | partielle korrelasjoner justert for gruppe | permutasjonstest for
#   kompositten (som fanger seleksjonseffekten)
#
# Merk: den opprinnelige analysen (005_ECV_phenotype_correlations.R:219)
# beregnet ALDRI p-verdier — kandidatene ble rangert på korrelasjons-
# koeffisient. Ved fast n er |rho| og p monotont relatert, så rangering på
# rho ER rangering på p; multiplisiteten ligger i seleksjonen.
#
# Krever: R1_00_setup.R
# =====================================================================
suppressPackageStartupMessages({library(data.table)})
source(file.path(Sys.getenv("ECV_ROOT", "/Users/sindrle/Research/Prosjekter/ECV"),
                 "review_round_1/01_scripts/R1_config.R"))
ROOT <- ECV_ROOT
OUT  <- file.path(ROOT, "review_round_1/02_output"); setwd(ROOT)
S <- readRDS(file.path(OUT, "R1_00_normalised.rds"))
set.seed(1)

PHENOS <- c("Weight","BMI","FFM","M_mg.kg.min","VO2max.kg","ChestPress",
            "PullDown","LegPress","mri_VL","mri_IP","mri_TAT","mri_SAT",
            "mrs_Hep","mrs_Musc")
CAND <- c("LTBP4","VCL","MEP1B","HSP90B1")

# ---- fenotypedeltaer --------------------------------------------------
ph <- as.data.table(S$pheno)
setnames(ph, "timepoint", "Time", skip_absent = TRUE)
pa <- ph[Time == "A1"][order(ID2)]; pb <- ph[Time == "B1"][order(ID2)]
ids <- intersect(pa$ID2, pb$ID2)
pa <- pa[match(ids, ID2)]; pb <- pb[match(ids, ID2)]
dph <- as.data.table(lapply(PHENOS, function(v) pb[[v]] - pa[[v]]))
setnames(dph, PHENOS); dph[, ID2 := ids][, Gr := pa$group]

# ---- proteindeltaer (loess, >=70% completeness) -----------------------
idx <- which(S$si$SampleCode %in% c("A1","B1"))
M <- S$Xloe[, idx]; s <- S$si[idx]
nA <- rowSums(!is.na(M[, s$SampleCode == "A1"])); nB <- rowSums(!is.na(M[, s$SampleCode == "B1"]))
M <- M[nA >= .7*26 & nB >= .7*26, ]; g <- S$genes[nA >= .7*26 & nB >= .7*26]
dp <- sapply(ids, function(i) {
  a <- M[, s$SampleCode == "A1" & s$ID2 == i, drop = FALSE]
  b <- M[, s$SampleCode == "B1" & s$ID2 == i, drop = FALSE]
  if (!ncol(a) || !ncol(b)) return(rep(NA_real_, nrow(M)))
  b[,1] - a[,1]
})
rownames(dp) <- make.unique(g)
dp <- dp[rowSums(!is.na(dp)) >= 15, ]
cat(sprintf("proteiner %d x fenotyper %d = %d tester | n deltakere %d\n",
            nrow(dp), length(PHENOS), nrow(dp)*length(PHENOS), length(ids)))

# ---- korrelasjonsmatrise med p -----------------------------------------
res <- rbindlist(lapply(PHENOS, function(v) {
  y <- dph[[v]]
  rbindlist(lapply(seq_len(nrow(dp)), function(i) {
    x <- dp[i, ]; ok <- !is.na(x) & !is.na(y)
    if (sum(ok) < 10) return(NULL)
    ct <- suppressWarnings(cor.test(x[ok], y[ok], method = "spearman"))
    data.table(phenotype = v, Genes = rownames(dp)[i], n = sum(ok),
               rho = unname(ct$estimate), p = ct$p.value)
  }))
}))
res[, q_global := p.adjust(p, "BH")]                       # over HELE matrisen
res[, q_within := p.adjust(p, "BH"), by = phenotype]       # per fenotype
setorder(res, p)
fwrite(res, file.path(OUT, "R1_02_phenotype_matrix.csv"))

cat("\n=== MULTIPLISITETSKONTROLL ===\n")
cat(sprintf("nominell p<0.05 : %d (forventet ved tilfeldighet ~%.0f)\n",
            sum(res$p < .05), .05*nrow(res)))
cat(sprintf("FDR<0.05 per fenotype : %d\n", sum(res$q_within < .05)))
cat(sprintf("FDR<0.05 hele matrisen: %d\n", sum(res$q_global < .05)))

cat("\n=== dVO2max ===\n")
v <- res[phenotype == "VO2max.kg"][order(p)]
cat(sprintf("nominell p<0.05: %d av %d (forventet ~%.0f)\n",
            sum(v$p < .05), nrow(v), .05*nrow(v)))
print(head(v[, .(Genes, rho = round(rho,3), p = signif(p,3),
                 q_within = signif(q_within,3), q_global = signif(q_global,3))], 15))
cat("\nKandidatene mot dVO2max:\n")
print(v[Genes %in% CAND, .(Genes, rho = round(rho,3), p = signif(p,3),
                           q_within = signif(q_within,3), q_global = signif(q_global,3))])
cat("\nTPI1 mot alle fenotyper (nominell p<0.05):\n")
print(res[Genes == "TPI1" & p < .05,
          .(phenotype, rho = round(rho,3), p = signif(p,3),
            q_within = signif(q_within,3), q_global = signif(q_global,3))])

# ---- partiell korrelasjon justert for gruppe --------------------------
cat("\n=== PARTIELL KORRELASJON justert for gruppe (rangbasert) ===\n")
grp <- as.numeric(factor(dph$Gr))
pc <- rbindlist(lapply(CAND, function(gn) {
  if (!gn %in% rownames(dp)) return(NULL)
  x <- rank(dp[gn, ]); y <- rank(dph$VO2max.kg); ok <- !is.na(dp[gn, ])
  m <- summary(lm(y[ok] ~ x[ok] + grp[ok]))
  data.table(Genes = gn, beta = coef(m)[2,1], p_adj_group = coef(m)[2,4],
             p_unadj = cor.test(dp[gn,ok], dph$VO2max.kg[ok],
                                method = "spearman")$p.value)
}))
print(pc[, .(Genes, beta = round(beta,3), p_unadj = signif(p_unadj,3),
             p_adj_group = signif(p_adj_group,3))])

# ---- permutasjonstest for komposittscoren -----------------------------
# Nullen må inkludere SELEKSJONEN: for hver permutasjon velges de fire
# beste på nytt mot den permuterte utfallsvariabelen.
cat("\n=== PERMUTASJONSTEST, firer-komposittscore mot dVO2max ===\n")
# NB: z-scorene må fortegnsjusteres før de midles, ellers kansellerer
# positive og negative korrelater hverandre og nullfordelingen blir ugyldig.
build <- function(y) {
  r <- apply(dp, 1, function(x) { ok <- !is.na(x) & !is.na(y)
    if (sum(ok) < 10) return(NA_real_)
    suppressWarnings(cor(x[ok], y[ok], method = "spearman")) })
  top <- names(sort(abs(r), decreasing = TRUE))[1:4]
  z <- t(scale(t(dp[top, , drop = FALSE]))) * sign(r[top])
  sc <- colMeans(z, na.rm = TRUE)
  ok <- !is.na(sc) & !is.na(y)
  list(top = top, rho = suppressWarnings(cor(sc[ok], y[ok], method = "spearman")))
}
obs <- build(dph$VO2max.kg)
cat("valgt på ekte data:", paste(obs$top, collapse = ", "), "\n")
cat(sprintf("observert rho = %+.3f\n", obs$rho))
NP <- 2000
null <- replicate(NP, build(sample(dph$VO2max.kg))$rho)
cat(sprintf("permutasjons-p = %.4f  (nullfordeling: median %+.3f, 95%%-kvantil %+.3f)\n",
            (1 + sum(abs(null) >= abs(obs$rho))) / (NP + 1),
            median(null, na.rm = TRUE), quantile(abs(null), .95, na.rm = TRUE)))
cat("\nSammenligning: den PUBLISERTE kompositten (LTBP4/VCL/MEP1B/HSP90B1)\n")
zz <- t(scale(t(dp[intersect(CAND, rownames(dp)), , drop = FALSE])))
sc <- colMeans(zz, na.rm = TRUE); ok <- !is.na(sc)
ctp <- suppressWarnings(cor.test(sc[ok], dph$VO2max.kg[ok], method = "spearman"))
cat(sprintf("  rho = %+.3f, nominell p = %.4f, mot samme permutasjonsnull p = %.4f\n",
            ctp$estimate, ctp$p.value,
            (1 + sum(abs(null) >= abs(ctp$estimate))) / (NP + 1)))

saveRDS(list(res = res, dp = dp, dph = dph, obs = obs, null = null, PHENOS = PHENOS),
        file.path(OUT, "R1_02_phenotypes.rds"))
