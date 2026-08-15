# =====================================================================
# R1_41_within_between_peptide.R — riktig håndtering av batch/peptid
#
# R1_39 la peptidkonsentrasjon inn som ÉN kovariat. Det blander to
# størrelser med ulik betydning, og gir en koeffisient som er et
# kompromiss mellom dem:
#
#   pep_between = deltakerens gjennomsnitt over injeksjoner
#   pep_within  = injeksjonens avvik fra eget gjennomsnitt
#
# Strukturen i dataene (verifisert):
#   - batch er KONSTANT innad i deltaker for 24 av 26 (kun ID 4 og 17
#     krysser). Batch kan derfor ikke konfundere en innen-deltaker-
#     kontrast; random intercept absorberer den allerede.
#   - peptid varierer derimot innad i deltaker (median spenn 0.080 mg/ml,
#     mellom-deltaker SD 0.081), og innen-delen skiller seg på tidspunkt
#     (+0.040 mg/ml, p = 0.024). DET er den som kan konfundere.
#   - batch henger sammen med pep_between (p = 0.059), ikke med
#     pep_within (p = 0.59). Kollineariteten ligger altså i den delen som
#     ikke påvirker tidspunktestimatet.
#
# Konklusjon: splitt peptid, behold batch. Da estimeres tidspunkteffekten
# justert for den lastvariasjonen som faktisk kan forstyrre den, uten at
# batch og peptid konkurrerer om samme varians.
#
# Krever: R1_00_setup.R, R1_01_longterm.R
# =====================================================================
suppressPackageStartupMessages({
  library(data.table); library(lme4); library(lmerTest); library(car)
})
source(file.path(Sys.getenv("ECV_ROOT", "/Users/sindrle/Research/Prosjekter/ECV"),
                 "review_round_1/01_scripts/R1_config.R"))
ROOT <- ECV_ROOT
OUT  <- file.path(ROOT, "review_round_1/02_output"); setwd(ROOT)
S <- readRDS(file.path(OUT, "R1_00_normalised.rds"))
THR <- 0.70

idx <- which(S$si$SampleCode %in% c("A1","B1"))
s <- S$si[idx]
mu <- tapply(s$peptide_mg_ml, s$ID2, mean)
s[, pep_between := as.numeric(mu[as.character(ID2)])]
s[, pep_within  := peptide_mg_ml - pep_between]

fit_all <- function(form) {
  m <- S$Xloe[, idx, drop = FALSE]
  nA <- rowSums(!is.na(m[, s$SampleCode == "A1", drop = FALSE]))
  nB <- rowSums(!is.na(m[, s$SampleCode == "B1", drop = FALSE]))
  ok <- nA >= THR*26 & nB >= THR*26
  m <- m[ok, , drop = FALSE]; g <- S$genes[ok]
  d0 <- data.table(SampleCode = factor(s$SampleCode, c("A1","B1")), Gr = s$Gr,
                   batch = s$batch, ID2 = factor(s$ID2),
                   pep_within = s$pep_within, pep_between = s$pep_between)
  rbindlist(lapply(seq_len(nrow(m)), function(i) {
    d <- copy(d0); d[, y := m[i, ]]; d <- d[!is.na(y)]
    if (uniqueN(d$SampleCode) < 2 || uniqueN(d$ID2) < 5) return(NULL)
    if (uniqueN(d$batch) < 2) return(NULL)
    f <- tryCatch(suppressWarnings(suppressMessages(
           lmer(form, data = d, REML = FALSE))), error = function(e) NULL)
    if (is.null(f)) return(NULL)
    co <- summary(f)$coefficients
    if (!"SampleCodeB1" %in% rownames(co)) return(NULL)
    a <- tryCatch(car::Anova(f, type = 2), error = function(e) NULL)
    if (is.null(a)) return(NULL)
    data.table(Genes = g[i], Estimate = co["SampleCodeB1","Estimate"],
               p = a["SampleCode","Pr(>Chisq)"])
  }))[, q := p.adjust(p, "BH")][order(p)]
}

mods <- list(
  "primaer  (batch)"                 = y ~ SampleCode + Gr + batch + (1|ID2),
  "R1_39    (batch + pep samlet)"    = y ~ SampleCode + Gr + batch + pep_within + pep_between + (1|ID2),
  "kun innen-deltaker peptid"        = y ~ SampleCode + Gr + batch + pep_within + (1|ID2),
  "uten batch, med splittet peptid"  = y ~ SampleCode + Gr + pep_within + pep_between + (1|ID2))

res <- list()
for (nm in names(mods)) {
  r <- fit_all(mods[[nm]]); res[[nm]] <- r
  cat(sprintf("%-32s testet %4d | FDR<0.05: %3d (opp %3d / ned %3d)\n",
              nm, nrow(r), sum(r$q<.05), sum(r$q<.05 & r$Estimate>0), sum(r$q<.05 & r$Estimate<0)))
}
base <- res[[1]][q<.05, Genes]
cat("\noverlapp med primaeranalysen:\n")
for (nm in names(res)[-1]) {
  g <- res[[nm]][q<.05, Genes]
  m <- merge(res[[1]][, .(Genes, e0=Estimate)], res[[nm]][, .(Genes, e1=Estimate)], by="Genes")
  cat(sprintf("  %-32s %3d delt av %3d | r(estimat) = %.4f\n",
              nm, length(intersect(base,g)), length(g), cor(m$e0, m$e1)))
}
saveRDS(res, file.path(OUT, "R1_41_within_between.rds"))
fwrite(rbindlist(lapply(names(res), function(nm)
  data.table(model=nm, n_tested=nrow(res[[nm]]), n_fdr=sum(res[[nm]]$q<.05),
             n_up=sum(res[[nm]]$q<.05 & res[[nm]]$Estimate>0),
             n_down=sum(res[[nm]]$q<.05 & res[[nm]]$Estimate<0)))),
  file.path(OUT, "R1_41_within_between.csv"))
cat("\nskrevet: R1_41_within_between.csv\n")
