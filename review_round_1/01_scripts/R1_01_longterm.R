# =====================================================================
# R1_01_longterm.R — 12-ukers responderliste
# Adresserer R5 punkt 1 og 9, R3 (statistisk rammeverk), R4 (Methods)
#
#   ingen imputering | eksplisitt completeness-filter | cyclic loess
#   | batch som MODELLTERM (ikke forhåndsfjernet) | lmer, BH-FDR
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

THR_PRIMARY <- 0.70   # completeness innen hvert tidspunkt

fit_all <- function(M, thr, use_batch = TRUE) {
  idx <- which(S$si$SampleCode %in% c("A1", "B1"))
  m <- M[, idx, drop = FALSE]; s <- S$si[idx]
  nA <- rowSums(!is.na(m[, s$SampleCode == "A1", drop = FALSE]))
  nB <- rowSums(!is.na(m[, s$SampleCode == "B1", drop = FALSE]))
  ok <- nA >= thr * 26 & nB >= thr * 26
  m <- m[ok, , drop = FALSE]; g <- S$genes[ok]
  d0 <- data.table(SampleCode = factor(s$SampleCode, c("A1", "B1")),
                   Gr = s$Gr, batch = s$batch, ID2 = factor(s$ID2))
  form <- if (use_batch) y ~ SampleCode + Gr + batch + (1 | ID2)
          else           y ~ SampleCode + Gr + (1 | ID2)
  res <- rbindlist(lapply(seq_len(nrow(m)), function(i) {
    d <- copy(d0); d[, y := m[i, ]]; d <- d[!is.na(y)]
    if (uniqueN(d$SampleCode) < 2 || uniqueN(d$ID2) < 5) return(NULL)
    if (use_batch && uniqueN(d$batch) < 2) return(NULL)
    f <- tryCatch(suppressWarnings(suppressMessages(
           lmer(form, data = d, REML = FALSE))), error = function(e) NULL)
    if (is.null(f)) return(NULL)
    co <- summary(f)$coefficients
    if (!"SampleCodeB1" %in% rownames(co)) return(NULL)
    a <- tryCatch(car::Anova(f, type = 2), error = function(e) NULL)
    if (is.null(a)) return(NULL)
    data.table(Protein.Group = rownames(m)[i], Genes = g[i],
               Estimate = co["SampleCodeB1","Estimate"],
               SE = co["SampleCodeB1","Std. Error"],
               p = a["SampleCode","Pr(>Chisq)"], n_ID = uniqueN(d$ID2))
  }))
  if (nrow(res)) res[, q := p.adjust(p, "BH")][order(p)] else res
}
summ <- function(r, lab) {
  s <- r[q < 0.05]
  cat(sprintf("%-42s testet %5d | FDR<0.05 %4d (opp %3d/ned %3d) | median est %+0.4f\n",
              lab, nrow(r), nrow(s), sum(s$Estimate > 0), sum(s$Estimate < 0),
              median(r$Estimate)))
  invisible(r)
}

cat("=== SENSITIVITET: normalisering x completeness (batch i modellen) ===\n")
sens <- list()
for (nm in c("rå","median","loess")) {
  M <- switch(nm, "rå" = S$X, "median" = S$Xmed, "loess" = S$Xloe)
  for (thr in c(.50,.70,.90,1.00)) {
    r <- fit_all(M, thr); sens[[paste(nm,thr)]] <- r
    summ(r, sprintf("%-6s | completeness >=%3.0f%%", nm, 100*thr))
  }
  cat("\n")
}
cat("=== BATCH-TERMENS BETYDNING (loess, >=70%) ===\n")
for (ub in c(TRUE, FALSE))
  summ(fit_all(S$Xloe, THR_PRIMARY, ub),
       sprintf("batch %s", ifelse(ub, "i modellen", "utelatt")))

PRIM <- sens[[paste("loess", THR_PRIMARY)]]
cat("\n=== PRIMÆRANALYSE ===\n"); summ(PRIM, "loess + >=70% + batch")
fwrite(PRIM, file.path(OUT, "R1_01_longterm_primary.csv"))

cat("\nKandidater:\n")
print(PRIM[Genes %in% c("LTBP4","VCL","MEP1B","HSP90B1","TPI1"),
           .(Genes, Estimate = round(Estimate,3), SE = round(SE,3),
             p = signif(p,3), q = signif(q,3), n_ID)])
cat("\nTopp 25:\n")
print(head(PRIM[, .(Genes, Estimate = round(Estimate,3),
                    p = signif(p,3), q = signif(q,3), n_ID)], 25))

old <- fread("A1_vs_B1_all_carAnova_type2.csv")
oldsig <- unique(old[q < 0.05 & p > 0 & !is.na(Genes), Genes])
newsig <- unique(PRIM[q < 0.05, Genes])
cat(sprintf("\npublisert %d | ny %d | overlapp %d | publiserte testbare her %d\n",
            length(oldsig), length(newsig), length(intersect(oldsig,newsig)),
            length(intersect(oldsig, unique(PRIM$Genes)))))
saveRDS(list(PRIM = PRIM, sens = sens, THR = THR_PRIMARY),
        file.path(OUT, "R1_01_longterm.rds"))
