# =====================================================================
# R1_00_setup.R  —  felles innlesing og normalisering
#
# Kjøres FØRST. Alle andre R1-skript leser cachen denne lager.
# Skriver kun til review_round_1/02_output/. Rører ingen originalfiler.
#
# Kilder (leses, endres ikke):
#   pg.report_Frode_plasma.csv        rå DIA-NN, 2971 pg x 104 prøver, 30.2% NA
#   Frode_Norheim_sample_info.csv     prøve -> ID/tidspunkt/plate/brønn/RunID
#   MyoGlu_pheno_2019-23-08.csv       fenotyper og gruppe
#
# Normalisering: cyclic loess mot radvis referanse, tilpasset på observerte
# par (NA-sikker), kjørt på alle 104 injeksjoner fordi den intensitets-
# avhengige forvrengningen er en egenskap ved injeksjonen, ikke ved kontrasten.
#
# ANTAKELSE som må stå i Methods: som all mellom-prøve-normalisering
# forutsetter loess at flertallet av proteiner ikke endrer seg. Absolutt
# mengde EP-materiale skal leses fra peptidkonsentrasjon/ELISA, ikke herfra.
# =====================================================================

suppressPackageStartupMessages({library(data.table)})
source(file.path(Sys.getenv("ECV_ROOT", "/Users/sindrle/Research/Prosjekter/ECV"),
                 "review_round_1/01_scripts/R1_config.R"))
ROOT <- ECV_ROOT
OUT  <- file.path(ROOT, "review_round_1/02_output")
setwd(ROOT); dir.create(OUT, showWarnings = FALSE, recursive = TRUE)


# ---- prøvemetadata ----------------------------------------------------
si <- fread("Frode_Norheim_sample_info.csv")
pg <- fread("pg.report_Frode_plasma.csv")
pg <- pg[, which(names(pg) != ""), with = FALSE]
samp <- intersect(si$original_colname, names(pg))
si   <- si[match(samp, original_colname)]

pheno <- read.table(PHENO_FILE, sep = ";", header = TRUE, row.names = 1)
pheno$ID2 <- as.numeric(as.character(pheno$ID))
# ID-transformasjon som i 002_ECV_vulcano.R:70 — første siffer droppes
si[, ID2 := as.numeric(substr(as.character(ID), 2, nchar(as.character(ID))))]
si[, Gr  := factor(pheno$group[match(ID2, pheno$ID2)], c("dysglyc", "control"))]
stopifnot(!any(is.na(si$Gr)))

# batch: nanoLC-havari deler injeksjonssekvensen rent i to.
# Utledet av differansen mellom Perseus-arkene 't-tests' og 'Before batch
# correction'; splittet samsvarer eksakt med RunID uten overlapp.
si[, batch := factor(ifelse(as.numeric(RunID) <= 18371, "A", "B"))]

# peptidkonsentrasjon (målt retrospektivt av kjernefasiliteten, aug 2026)
pepA <- rbind(
  c(0.38,0.54,0.34,0.48,0.31,0.27,0.35,0.29,0.30,0.25,0.26,0.29),
  c(0.34,0.30,0.47,0.46,0.32,0.24,0.37,0.40,0.15,0.23,0.24,0.21),
  c(0.29,0.32,0.27,0.42,0.26,0.27,0.32,0.38,0.27,0.38,0.16,0.24),
  c(0.28,0.24,0.28,0.35,0.34,0.35,0.43,0.40,0.32,0.23,0.51,0.18),
  c(0.38,0.55,0.29,0.31,0.35,0.44,0.36,0.36,0.24,0.38,0.34,0.37),
  c(0.40,0.30,0.44,0.51,0.37,0.37,0.37,0.44,0.35,0.54,0.51,0.46),
  c(0.37,0.39,0.42,0.23,0.37,0.37,0.25,0.41,0.26,0.10,0.08,0.19))
rownames(pepA) <- LETTERS[1:7]
pepH  <- c(0.39,0.26,0.32,0.27,0.33,0.31,0.33,0.32,0.40,0.35,0.52,0.47) # plate A rad H
pepS2 <- c(0.39,0.33,0.40,0.34,0.40,0.52,0.54,0.48)                     # plate S2 A1-A8
pep <- function(plate, well) {
  r <- substr(well, 1, 1); cc <- as.integer(substr(well, 2, nchar(well)))
  if (plate == "S2") return(pepS2[cc])
  if (r == "H") return(pepH[cc])
  pepA[r, cc]
}
si[, peptide_mg_ml := mapply(pep, Plate, Well)]
stopifnot(!any(is.na(si$peptide_mg_ml)))

# ---- log2-matrise -----------------------------------------------------
X <- as.matrix(pg[, ..samp]); storage.mode(X) <- "double"
X[X <= 0] <- NA; X <- log2(X)
rownames(X) <- pg$Protein.Group
genes <- pg$Genes
keep <- !is.na(genes) & genes != ""
X <- X[keep, ]; genes <- genes[keep]

# ---- normalisering ----------------------------------------------------
loess_normalise <- function(M, iter = 3, span = 0.7) {
  for (k in seq_len(iter)) {
    ref <- apply(M, 1, median, na.rm = TRUE)
    for (j in seq_len(ncol(M))) {
      ok <- !is.na(M[, j]) & !is.na(ref)
      if (sum(ok) < 100) next
      a <- (M[ok, j] + ref[ok]) / 2; m <- M[ok, j] - ref[ok]
      fit <- tryCatch(loess(m ~ a, span = span, degree = 1,
                            family = "symmetric", surface = "direct"),
                      error = function(e) NULL)
      if (!is.null(fit)) M[ok, j] <- M[ok, j] - predict(fit, data.frame(a = a))
    }
  }
  M
}
cm   <- apply(X, 2, median, na.rm = TRUE)
Xmed <- sweep(X, 2, cm - mean(cm))
message("cyclic loess ...")
Xloe <- loess_normalise(X)

saveRDS(list(X = X, Xmed = Xmed, Xloe = Xloe, genes = genes, si = si,
             pheno = pheno, PHENO_FILE = PHENO_FILE),
        file.path(OUT, "R1_00_normalised.rds"))

cat(sprintf("prøver %d | proteingrupper %d | NA %.1f%%\n",
            nrow(si), nrow(X), 100 * mean(is.na(X))))
print(table(si$SampleCode, si$batch))
cat(sprintf("peptidkons A1 %.3f vs B1 %.3f mg/ml\n",
            mean(si[SampleCode == "A1"]$peptide_mg_ml),
            mean(si[SampleCode == "B1"]$peptide_mg_ml)))
cat("skrevet:", file.path(OUT, "R1_00_normalised.rds"), "\n")
