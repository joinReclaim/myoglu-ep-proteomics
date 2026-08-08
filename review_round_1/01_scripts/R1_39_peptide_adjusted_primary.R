# =====================================================================
# R1_39_peptide_adjusted_primary.R — holder responderlisten når vi
# justerer for målt peptidmengde per injeksjon?
#
# Bakgrunn: verifiseringen av R1_33 viste at peptidlast IKKE forklarer
# A1–B1-forskjellen i identifikasjonsdybde (18 % av gapet; tidspunkt-
# effekten overlever justering, p = 0.038), og at lastforskjellen mellom
# A1 og B1 ikke er signifikant (Wilcoxon p = 0.12). Vi kan derfor ikke
# påstå en lastmekanisme. Det vi KAN gjøre er å vise at konklusjonene
# ikke avhenger av last i det hele tatt — en direkte sensitivitetsanalyse
# framfor en mekanismepåstand.
#
# Krever: R1_00_setup.R
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

fit_all <- function(M, use_pep) {
  idx <- which(S$si$SampleCode %in% c("A1","B1"))
  m <- M[, idx, drop = FALSE]; s <- S$si[idx]
  nA <- rowSums(!is.na(m[, s$SampleCode == "A1", drop = FALSE]))
  nB <- rowSums(!is.na(m[, s$SampleCode == "B1", drop = FALSE]))
  ok <- nA >= THR*26 & nB >= THR*26
  m <- m[ok, , drop = FALSE]; g <- S$genes[ok]
  d0 <- data.table(SampleCode = factor(s$SampleCode, c("A1","B1")),
                   Gr = s$Gr, batch = s$batch, ID2 = factor(s$ID2),
                   pep = as.numeric(s$peptide_mg_ml))
  form <- if (use_pep) y ~ SampleCode + Gr + batch + pep + (1 | ID2)
          else         y ~ SampleCode + Gr + batch + (1 | ID2)
  res <- rbindlist(lapply(seq_len(nrow(m)), function(i) {
    d <- copy(d0); d[, y := m[i, ]]; d <- d[!is.na(y)]
    if (uniqueN(d$SampleCode) < 2 || uniqueN(d$ID2) < 5) return(NULL)
    if (uniqueN(d$batch) < 2) return(NULL)
    # NB: fang KUN error. lmer advarer rutinemessig; warning-fangst dreper alt.
    f <- tryCatch(suppressWarnings(suppressMessages(
           lmer(form, data = d, REML = FALSE))), error = function(e) NULL)
    if (is.null(f)) return(NULL)
    co <- summary(f)$coefficients
    if (!"SampleCodeB1" %in% rownames(co)) return(NULL)
    a <- tryCatch(car::Anova(f, type = 2), error = function(e) NULL)
    if (is.null(a)) return(NULL)
    data.table(Genes = g[i], Estimate = co["SampleCodeB1","Estimate"],
               SE = co["SampleCodeB1","Std. Error"],
               p = a["SampleCode","Pr(>Chisq)"])
  }))
  res[, q := p.adjust(p, "BH")][order(p)]
}

cat("=== uten peptidjustering (= primæranalysen i R1_01) ===\n")
r0 <- fit_all(S$Xloe, FALSE)
cat(sprintf("testet %d | FDR<0.05: %d (opp %d / ned %d)\n",
            nrow(r0), sum(r0$q<.05), sum(r0$q<.05 & r0$Estimate>0), sum(r0$q<.05 & r0$Estimate<0)))

cat("\n=== MED peptidkonsentrasjon som kovariat ===\n")
r1 <- fit_all(S$Xloe, TRUE)
cat(sprintf("testet %d | FDR<0.05: %d (opp %d / ned %d)\n",
            nrow(r1), sum(r1$q<.05), sum(r1$q<.05 & r1$Estimate>0), sum(r1$q<.05 & r1$Estimate<0)))

s0 <- r0[q<.05, Genes]; s1 <- r1[q<.05, Genes]
cat(sprintf("\noverlapp: %d | kun ujustert: %d | kun justert: %d | Jaccard %.3f\n",
            length(intersect(s0,s1)), length(setdiff(s0,s1)), length(setdiff(s1,s0)),
            length(intersect(s0,s1))/length(union(s0,s1))))

m <- merge(r0[, .(Genes, e0=Estimate, q0=q)], r1[, .(Genes, e1=Estimate, q1=q)], by="Genes")
cat(sprintf("effektestimater: Pearson r = %.4f | Spearman rho = %.4f | samme fortegn %.1f%%\n",
            cor(m$e0, m$e1), cor(m$e0, m$e1, method="spearman"),
            100*mean(sign(m$e0)==sign(m$e1))))
cat(sprintf("median |endring i estimat| = %.4f (median |estimat| = %.4f)\n",
            median(abs(m$e1-m$e0)), median(abs(m$e0))))

fwrite(m, file.path(OUT, "R1_39_peptide_adjusted.csv"))
saveRDS(list(r0=r0, r1=r1), file.path(OUT, "R1_39_peptide_adjusted.rds"))
cat("\nskrevet: R1_39_peptide_adjusted.csv\n")
