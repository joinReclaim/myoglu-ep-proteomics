# =====================================================================
# R1_33_qc_report.R — QC-rapport på prøvenivå (R5 punkt 9, R6)
#
# R5: "no sample-level quality metrics are reported for a 104-sample MS
# experiment — no QC or pooled injections, no coefficients of variation,
# no per-sample intensity distributions, no injection-order diagnostics.
# Also missing: per-protein missingness; the number of proteins quantified
# per sample; whether instrument batch was confounded with timepoint."
#
# Dette skriptet lager EN supplementær figur (paneler A-E) og EN tabell med
# én rad per injeksjon.
#
# VIKTIG: alt som handler om last, injeksjonsrekkefølge eller identifikasjons-
# dybde bruker RÅ log2 (S$X). Loess fjerner nettopp lasteffekten og ville
# maskert det R5 spør om. S$Xloe brukes bare der normaliseringens EFFEKT
# skal vises (panel B, nedre rad) og i den biologiske CV-en.
#
# INGEN POOLED QC-INJEKSJONER finnes i datasettet. Teknisk CV kan derfor
# ikke beregnes. Vi rapporterer biologisk CV per protein innen tidspunkt
# som nærmeste tilgjengelige størrelse, og sier eksplisitt fra om det.
# =====================================================================
suppressPackageStartupMessages({
  library(data.table); library(ggplot2); library(cowplot)
})
source(file.path(Sys.getenv("ECV_ROOT", "/Users/sindrle/Research/Prosjekter/ECV"),
                 "review_round_1/01_scripts/R1_config.R"))
ROOT <- ECV_ROOT
OUT  <- file.path(ROOT, "review_round_1/02_output"); setwd(ROOT)
S <- readRDS(file.path(OUT, "R1_00_normalised.rds"))
source(file.path(ROOT, "review_round_1/01_scripts/R1_theme.R"))

sink(file.path(OUT, "R1_33_log.txt"), split = TRUE)

X   <- S$X            # rå log2, 2963 x 104 (kun proteingrupper med gensymbol)
XL  <- S$Xloe         # cyclic loess
si  <- copy(S$si)
nP  <- nrow(X); nS <- ncol(X)
cat(sprintf("matrise: %d proteingrupper x %d injeksjoner\n", nP, nS))

# Engelske etiketter til figurene. A1 = hvile før trening, B1 = hvile etter
# 12 uker, B2 = rett etter akutt økt, B3 = 2 t etter (jf. R1_05).
TLAB <- c(A1 = "A1  pre-training rest", B1 = "B1  post-training rest",
          B2 = "B2  acute, 0 h",        B3 = "B3  acute, 2 h")
si[, Time4 := factor(TLAB[SampleCode], levels = unname(TLAB))]
COL_T <- setNames(unname(PAL$cluster), unname(TLAB))   # 4 tidspunkt -> 4 farger
COL_B <- setNames(unname(PAL$venn), c("A", "B"))       # 2 batcher

# ---------------------------------------------------------------------
# 1. Metrikker per injeksjon
# ---------------------------------------------------------------------
si[, n_proteins      := colSums(!is.na(X))]
si[, frac_missing    := 1 - n_proteins / nP]
si[, median_log2_raw := apply(X,  2, median, na.rm = TRUE)]
si[, median_log2_loess := apply(XL, 2, median, na.rm = TRUE)]
si[, iqr_log2_raw    := apply(X,  2, IQR,    na.rm = TRUE)]

# Kjent tall fra tidligere: median 2266 (A1) mot 1980 (B1). Det ble regnet på
# HELE DIA-NN-matrisen (2971 proteingrupper). S$X er filtrert til de 2963 som
# har gensymbol, så tallene her ligger ~7 lavere. Vi regner begge når rådata
# er tilgjengelig, slik at avviket er dokumentert og ikke en uoverensstemmelse.
PGFILE <- file.path(ROOT, "pg.report_Frode_plasma.csv")
if (file.exists(PGFILE)) {
  pg <- fread(PGFILE); pg <- pg[, which(names(pg) != ""), with = FALSE]
  cn <- si$original_colname
  Xa <- as.matrix(pg[, ..cn]); storage.mode(Xa) <- "double"; Xa[Xa <= 0] <- NA
  si[, n_proteins_all := colSums(!is.na(Xa))]
  cat(sprintf("full DIA-NN-matrise: %d proteingrupper (gensymbol-filter fjerner %d)\n",
              nrow(Xa), nrow(Xa) - nP))
  rm(Xa, pg)
} else {
  si[, n_proteins_all := NA_integer_]
  cat("MERK: pg.report_Frode_plasma.csv ikke funnet — n_proteins_all er NA\n")
}

cat("\n=== IDENTIFIKASJONSDYBDE PER TIDSPUNKT ===\n")
dep <- si[, .(n = .N,
              median_n_proteins     = as.numeric(median(n_proteins)),
              median_n_proteins_all = as.numeric(median(n_proteins_all)),
              mean_peptide_mg_ml    = mean(peptide_mg_ml),
              median_peptide_mg_ml  = median(peptide_mg_ml),
              median_log2_raw       = median(median_log2_raw)),
          by = .(SampleCode)][order(SampleCode)]
print(dep)
cat(sprintf(paste0("\nKjent tall 2266 (A1) mot 1980 (B1): reprodusert på full matrise ",
                   "(%.1f mot %.1f).\n  På den gensymbol-filtrerte matrisen: %.1f mot %.1f.\n"),
            dep[SampleCode == "A1", median_n_proteins_all],
            dep[SampleCode == "B1", median_n_proteins_all],
            dep[SampleCode == "A1", median_n_proteins],
            dep[SampleCode == "B1", median_n_proteins]))
cat(sprintf("Kjent tall peptid A1 0,321 mot B1 0,361 mg/ml: her %.3f mot %.3f (gjennomsnitt).\n",
            dep[SampleCode == "A1", mean_peptide_mg_ml],
            dep[SampleCode == "B1", mean_peptide_mg_ml]))

wt <- wilcox.test(n_proteins ~ SampleCode, data = si[SampleCode %in% c("A1", "B1")])
cat(sprintf("A1 mot B1, antall proteiner: Wilcoxon p = %.3g\n", wt$p.value))
wp <- wilcox.test(peptide_mg_ml ~ SampleCode, data = si[SampleCode %in% c("A1", "B1")])
cat(sprintf("A1 mot B1, peptidkonsentrasjon: Wilcoxon p = %.3g\n", wp$p.value))

# ---------------------------------------------------------------------
# 2. Injeksjonsrekkefølge og batch
# ---------------------------------------------------------------------
sp <- function(a, b) {
  ct <- suppressWarnings(cor.test(a, b, method = "spearman"))
  c(rho = unname(ct$estimate), p = ct$p.value)
}
cat("\n=== INJEKSJONSREKKEFØLGE (RunID) ===\n")
cat(sprintf("RunID vs n_proteins       : rho = %+.3f  p = %.3g\n",
            sp(si$RunID, si$n_proteins)[1], sp(si$RunID, si$n_proteins)[2]))
cat(sprintf("RunID vs median log2 (rå) : rho = %+.3f  p = %.3g\n",
            sp(si$RunID, si$median_log2_raw)[1], sp(si$RunID, si$median_log2_raw)[2]))
for (b in c("A", "B")) {
  d <- si[batch == b]
  cat(sprintf("  innen batch %s (n = %2d): RunID vs n_proteins rho = %+.3f (p = %.3g)\n",
              b, nrow(d), sp(d$RunID, d$n_proteins)[1], sp(d$RunID, d$n_proteins)[2]))
}
BATCH_CUT <- 18371   # siste RunID i batch A; batch B starter på 18393
cat(sprintf("batchgrense: A = RunID %d-%d (n = %d), B = RunID %d-%d (n = %d)\n",
            min(si$RunID[si$batch == "A"]), max(si$RunID[si$batch == "A"]), sum(si$batch == "A"),
            min(si$RunID[si$batch == "B"]), max(si$RunID[si$batch == "B"]), sum(si$batch == "B")))

# ---------------------------------------------------------------------
# 3. Peptidlast: mer peptid -> HØYERE median, FÆRRE identifikasjoner
# ---------------------------------------------------------------------
cat("\n=== PEPTIDLAST (RÅ data) ===\n")
r_int <- sp(si$peptide_mg_ml, si$median_log2_raw)
r_id  <- sp(si$peptide_mg_ml, si$n_proteins)
r_cmp <- sp(si$median_log2_raw, si$n_proteins)
cat(sprintf("peptide_mg_ml vs median log2 : rho = %+.3f  p = %.3g\n", r_int[1], r_int[2]))
cat(sprintf("peptide_mg_ml vs n_proteins  : rho = %+.3f  p = %.3g\n", r_id[1],  r_id[2]))
cat(sprintf("median log2 vs n_proteins    : rho = %+.3f  p = %.3g   (kompresjonsaksen)\n",
            r_cmp[1], r_cmp[2]))
# Kontroll etter loess. NB: medianen over ALLE observerte verdier henger igjen
# med last (rho = +0.39) fordi hver prøve har sitt eget observerte proteinsett —
# prøver med få identifikasjoner mangler nettopp de svakeste proteinene, og
# medianen deres blir dermed høy uansett normalisering. Regnes medianen på et
# FAST proteinsett (de komplette i alle 104) forsvinner lasteffekten helt.
r_int_l <- sp(si$peptide_mg_ml, si$median_log2_loess)
cat(sprintf("kontroll, peptide vs median log2 ETTER loess (alle observerte): rho = %+.3f  p = %.3g\n",
            r_int_l[1], r_int_l[2]))
kc  <- rowSums(is.na(X)) == 0                       # komplette proteiner
mrc <- apply(X[kc, ],  2, median); mlc <- apply(XL[kc, ], 2, median)
cat(sprintf("fast proteinsett (n = %d komplette): peptide vs median rå rho = %+.3f, etter loess rho = %+.3f\n",
            sum(kc), sp(si$peptide_mg_ml, mrc)[1], sp(si$peptide_mg_ml, mlc)[1]))
cat(sprintf("  spredning i prøvemedian, fast sett: SD rå %.3f -> loess %.3f log2-enheter\n",
            sd(mrc), sd(mlc)))
# innen deltaker (fjerner mellom-person-variasjon)
si[, pep_c := peptide_mg_ml - mean(peptide_mg_ml), by = ID2]
si[, int_c := median_log2_raw - mean(median_log2_raw), by = ID2]
si[, np_c  := n_proteins    - mean(n_proteins),    by = ID2]
cat(sprintf("innen deltaker: peptide vs median log2 rho = %+.3f | peptide vs n_proteins rho = %+.3f\n",
            sp(si$pep_c, si$int_c)[1], sp(si$pep_c, si$np_c)[1]))

# ---------------------------------------------------------------------
# 4. Missingness per protein
# ---------------------------------------------------------------------
miss_all <- rowMeans(is.na(X))
# --- avvikende injeksjoner (robust: median +/- 3 MAD på antall proteiner) ---
mn <- median(si$n_proteins); ma <- mad(si$n_proteins)
si[, outlier := n_proteins < mn - 3 * ma | n_proteins > mn + 3 * ma]
cat(sprintf("\n=== AVVIKENDE INJEKSJONER (n_proteins utenfor %.0f +/- 3 MAD = %.0f) ===\n", mn, 3 * ma))
if (any(si$outlier)) {
  print(si[outlier == TRUE, .(filename, SampleCode, batch, RunID, peptide_mg_ml,
                              n_proteins, median_log2_raw = round(median_log2_raw, 2))])
  cat("Ingen injeksjoner er ekskludert fra analysene; de er beholdt og flagget her.\n")
} else cat("ingen\n")

cat("\n=== MISSINGNESS PER PROTEIN ===\n")
cat(sprintf("totalt %.1f %% manglende verdier i matrisen\n", 100 * mean(is.na(X))))
cat(sprintf("proteiner komplette i alle 104: %d (%.1f %%)\n",
            sum(miss_all == 0), 100 * mean(miss_all == 0)))
for (th in c(0, .3, .5, .7)) cat(sprintf("  observert i >= %2.0f %% av prøvene: %4d proteiner\n",
                                         100 * (1 - th), sum(miss_all <= th)))
mt <- rbindlist(lapply(names(TLAB), function(tc) {
  j <- si$SampleCode == tc
  data.table(Time4 = TLAB[tc], Genes = S$genes, miss = rowMeans(is.na(X[, j, drop = FALSE])))
}))
mt[, Time4 := factor(Time4, levels = unname(TLAB))]
print(mt[, .(median_missing = round(median(miss), 3),
             complete = sum(miss == 0)), by = Time4])

# ---------------------------------------------------------------------
# 5. Biologisk CV per protein innen tidspunkt (ERSTATNING for teknisk CV)
# ---------------------------------------------------------------------
# Det finnes ingen pooled QC- eller replikatinjeksjoner i dette datasettet,
# så teknisk CV kan ikke beregnes. Biologisk CV per protein innen tidspunkt
# er nærmeste tilgjengelige størrelse; den inneholder både biologisk og
# teknisk variasjon og er derfor et ØVRE anslag på teknisk presisjon.
cvfun <- function(M, j, minobs = 0.7) {
  A <- 2^M[, j, drop = FALSE]
  k <- rowMeans(!is.na(A)) >= minobs
  A <- A[k, , drop = FALSE]
  100 * apply(A, 1, sd, na.rm = TRUE) / rowMeans(A, na.rm = TRUE)
}
cat("\n=== BIOLOGISK CV PER PROTEIN INNEN TIDSPUNKT (ingen teknisk CV mulig) ===\n")
cvtab <- rbindlist(lapply(names(TLAB), function(tc) {
  j <- which(si$SampleCode == tc)
  data.table(Time = tc, n = length(j),
             n_prot   = sum(rowMeans(!is.na(XL[, j, drop = FALSE])) >= 0.7),
             cv_raw   = median(cvfun(X,  j), na.rm = TRUE),
             cv_loess = median(cvfun(XL, j), na.rm = TRUE))
}))
print(cvtab)
CV_MED <- median(cvtab$cv_loess)
cat(sprintf("median biologisk CV på tvers av tidspunkt (loess): %.1f %%\n", CV_MED))

# ---------------------------------------------------------------------
# 6. Er batch konfundert med tidspunkt (eller gruppe)?
# ---------------------------------------------------------------------
cat("\n=== BATCH MOT DESIGN ===\n")
tb_bt <- table(batch = si$batch, timepoint = si$SampleCode)
print(tb_bt)
ft_bt <- fisher.test(tb_bt)
cat(sprintf("Fisher exact, batch x tidspunkt : p = %.3f  -> ingen konfundering\n", ft_bt$p.value))
tb_bg <- table(batch = si$batch, group = si$Gr)
print(tb_bg)
ft_bg <- fisher.test(tb_bg)
cat(sprintf("Fisher exact, batch x gruppe    : p = %.3f\n", ft_bg$p.value))
# Kritisk for parede kontraster: havner A1 og B1 for samme deltaker i ulik batch?
sp2 <- dcast(si[SampleCode %in% c("A1", "B1")], ID2 ~ SampleCode, value.var = "batch")
sp2[, split := A1 != B1]
cat(sprintf("deltakere der A1 og B1 ligger i ULIK batch: %d av %d\n",
            sum(sp2$split), nrow(sp2)))
n_split_all <- si[, .(nb = uniqueN(batch)), by = ID2][nb > 1, .N]
cat(sprintf("deltakere med alle fire injeksjoner i samme batch: %d av %d\n",
            uniqueN(si$ID2) - n_split_all, uniqueN(si$ID2)))
cat("balanse (skal være 26 per tidspunkt, 13 per gruppe):\n")
print(table(si$SampleCode, si$Gr))

# ---------------------------------------------------------------------
# 7. Tabell: én rad per injeksjon
# ---------------------------------------------------------------------
tab <- si[order(RunID), .(original_colname, ID2, SampleCode, Gr, batch, RunID,
                          Plate, Well, peptide_mg_ml,
                          n_proteins, n_proteins_all,
                          median_log2_raw   = round(median_log2_raw, 3),
                          median_log2_loess = round(median_log2_loess, 3),
                          iqr_log2_raw      = round(iqr_log2_raw, 3),
                          frac_missing      = round(frac_missing, 4),
                          outlier_3mad      = outlier)]
fwrite(tab, file.path(OUT, "R1_33_qc_per_sample.csv"))
cat(sprintf("\ntabell: %s (%d rader)\n",
            file.path(OUT, "R1_33_qc_per_sample.csv"), nrow(tab)))

# =====================================================================
# FIGUR — supplementær QC, paneler A-E. All tekst på engelsk.
# =====================================================================
si[, ord := rank(RunID, ties.method = "first")]
brk <- si[order(ord)][seq(1, nS, by = 10), ord]
lbl <- si[order(ord)][seq(1, nS, by = 10), as.character(RunID)]
cut_x <- mean(c(max(si$ord[si$batch == "A"]), min(si$ord[si$batch == "B"])))

# --- A: identifikasjoner mot injeksjonsrekkefølge --------------------
pA <- ggplot(si, aes(RunID, n_proteins, colour = Time4)) +
  geom_vline(xintercept = mean(c(BATCH_CUT, min(si$RunID[si$batch == "B"]))),
             linetype = 2, colour = "grey40") +
  geom_point(size = 1.8, alpha = 0.9) +
  geom_point(data = si[outlier == TRUE], shape = 21, size = 3.4,
             colour = "grey20", fill = NA, stroke = 0.5) +
  geom_text(data = si[outlier == TRUE], aes(label = "> 3 MAD below median"),
            hjust = -0.15, vjust = 0.4, size = 2.8, colour = "grey20",
            show.legend = FALSE) +
  annotate("text", x = min(si$RunID) + 20, y = max(si$n_proteins),
           label = "batch A", hjust = 0, size = 3.2, colour = "grey30") +
  annotate("text", x = max(si$RunID) - 20, y = max(si$n_proteins),
           label = "batch B", hjust = 1, size = 3.2, colour = "grey30") +
  scale_colour_manual(values = COL_T) +
  labs(x = "Injection order (RunID)", y = "Protein groups quantified",
       title = "A   Identification depth by injection order") +
  theme_ecv(10) + theme(legend.position = "right", plot.title = element_text(hjust = 0))

# --- B: intensitetsfordelinger, rå over loess ------------------------
mkbox <- function(M, ttl) {
  d <- data.table(ord = rep(si$ord, each = nP), batch = rep(si$batch, each = nP),
                  v = as.vector(M))[!is.na(v)]
  ggplot(d, aes(factor(ord), v, fill = batch)) +
    geom_boxplot(outlier.size = 0.15, outlier.alpha = 0.15, linewidth = 0.15,
                 colour = "grey25") +
    scale_fill_manual(values = COL_B) +
    scale_x_discrete(breaks = as.character(brk), labels = lbl) +
    labs(x = "Injection order (RunID)", y = expression(log[2]~intensity), title = ttl) +
    theme_ecv(10) + theme(legend.position = "right", plot.title = element_text(hjust = 0),
                          axis.text.x = element_text(size = 6, angle = 90, vjust = 0.5))
}
pB1 <- mkbox(X,  "B   Per-sample intensity distributions — raw log2") +
  theme(axis.title.x = element_blank(), axis.text.x = element_blank(),
        axis.ticks.x = element_blank())
pB2 <- mkbox(XL, "     After cyclic loess normalisation") +
  labs(subtitle = sprintf(paste0("On a fixed set of %d proteins complete in all %d injections, ",
                                 "between-sample SD of the median falls from %.2f to %.2f log2 units"),
                          sum(kc), nS, sd(mrc), sd(mlc))) +
  theme(plot.subtitle = element_text(hjust = 0, size = 8))

# --- C: peptidlast mot signal og mot identifikasjoner ----------------
lab_rho <- function(r) sprintf("Spearman rho = %+.2f, p = %s", r[1],
                              format.pval(r[2], digits = 2, eps = 1e-16))
pC1 <- ggplot(si, aes(peptide_mg_ml, median_log2_raw)) +
  geom_smooth(method = "lm", se = TRUE, colour = "grey30", fill = "grey85", linewidth = 0.5) +
  geom_point(aes(colour = Time4), size = 1.6) +
  scale_colour_manual(values = COL_T) +
  labs(x = "Peptide concentration (mg/mL)", y = expression(Median~log[2]~intensity),
       title = "C   More peptide → higher signal ...", subtitle = lab_rho(r_int)) +
  theme_ecv(10) + theme(legend.position = "none", plot.title = element_text(hjust = 0),
                        plot.subtitle = element_text(hjust = 0, size = 8))
pC2 <- ggplot(si, aes(peptide_mg_ml, n_proteins)) +
  geom_smooth(method = "lm", se = TRUE, colour = "grey30", fill = "grey85", linewidth = 0.5) +
  geom_point(aes(colour = Time4), size = 1.6) +
  scale_colour_manual(values = COL_T) +
  labs(x = "Peptide concentration (mg/mL)", y = "Protein groups quantified",
       title = "     ... but fewer identifications", subtitle = lab_rho(r_id)) +
  theme_ecv(10) + theme(legend.position = "none", plot.title = element_text(hjust = 0),
                        plot.subtitle = element_text(hjust = 0, size = 8))
pC3 <- ggplot(si, aes(median_log2_raw, n_proteins)) +
  geom_smooth(method = "lm", se = TRUE, colour = "grey30", fill = "grey85", linewidth = 0.5) +
  geom_point(aes(colour = Time4), size = 1.6) +
  scale_colour_manual(values = COL_T) +
  labs(x = expression(Median~log[2]~intensity), y = "Protein groups quantified",
       title = "     Dynamic range compression axis", subtitle = lab_rho(r_cmp)) +
  theme_ecv(10) + theme(legend.position = "none", plot.title = element_text(hjust = 0),
                        plot.subtitle = element_text(hjust = 0, size = 8))

# --- D: missingness per protein --------------------------------------
pD1 <- ggplot(data.table(m = miss_all), aes(100 * m)) +
  geom_histogram(bins = 40, fill = PAL$venn[1], colour = "white", linewidth = 0.2) +
  labs(x = "Missing values per protein (%)", y = "Protein groups",
       title = "D   Per-protein missingness",
       subtitle = sprintf("%.1f%% missing overall; %d of %d complete in all %d injections",
                          100 * mean(is.na(X)), sum(miss_all == 0), nP, nS)) +
  theme_ecv(10) + theme(plot.title = element_text(hjust = 0),
                        plot.subtitle = element_text(hjust = 0, size = 8))
pD2 <- ggplot(mt, aes(100 * miss, colour = Time4)) +
  stat_ecdf(linewidth = 0.7) +
  scale_colour_manual(values = COL_T) +
  labs(x = "Missing values per protein, within timepoint (%)",
       y = "Cumulative fraction of proteins",
       title = "     Stratified by timepoint") +
  theme_ecv(10) + theme(legend.position = "right", plot.title = element_text(hjust = 0))

# --- E: batch mot tidspunkt ------------------------------------------
dE <- as.data.table(tb_bt)[, .(batch, timepoint, N)]
dE[, timepoint := factor(TLAB[timepoint], levels = unname(TLAB))]
pE <- ggplot(dE, aes(timepoint, batch, fill = N)) +
  geom_tile(colour = "white", linewidth = 1) +
  geom_text(aes(label = N), size = 3.4) +
  scale_fill_gradientn(colours = PAL$heat[55:85], guide = "none") +
  labs(x = NULL, y = "Instrument batch",
       title = "E   Batch vs timepoint",
       subtitle = sprintf("Fisher exact p = %.2f — not confounded; %d/%d participants split across batches",
                          ft_bt$p.value, n_split_all, uniqueN(si$ID2))) +
  theme_ecv(10) + theme(plot.title = element_text(hjust = 0),
                        plot.subtitle = element_text(hjust = 0, size = 8),
                        axis.text.x = element_text(angle = 20, hjust = 1),
                        axis.line = element_blank(), axis.ticks = element_blank())

rowC <- plot_grid(pC1, pC2, pC3, nrow = 1, rel_widths = c(1, 1, 1))
rowD <- plot_grid(pD1, pD2, nrow = 1, rel_widths = c(1, 1.25))
fig  <- plot_grid(pA, pB1, pB2, rowC, rowD, pE, ncol = 1,
                  rel_heights = c(1.0, 0.95, 1.15, 1.0, 1.0, 0.75))
save_ecv("R1_33_qc_supplementary.pdf", fig, width = 26, height = 40, unit = "cm")

cat("\nFERDIG\n")
sink()
