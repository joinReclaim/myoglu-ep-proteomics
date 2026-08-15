# =====================================================================
# R1_36_overlap_recompute.R — overlapp akutt vs kronisk, regnet om
#
# Bakgrunn: det innsendte manuskriptet oppga overlappen som 21,6 %.
# Det tallet var overlapp / UNION. Referentene påpekte at nevneren var
# feil for spørsmålet som ble stilt; rettelsen i tilsvaret var
# 137/314 = 43,6 %. BEGGE tallene er nå utdaterte — reanalysen
# (loess, batch som modellterm, >=70 % completeness, ingen imputering)
# ga nye responderlister.
#
# Dette skriptet regner overlappen på alle fire rimelige måter og
# skriver teller/nevner eksplisitt, slik at ingen nevner kan forveksles
# igjen. I tillegg: hypergeometrisk test mot det universet som faktisk
# ble TESTET i begge analysene, og retningskonkordans blant de
# overlappende.
#
# Krever: R1_01_longterm.R (R1_01_longterm_primary.csv)
#         R1_05_acute_chronic.R (R1_05_acute.csv)
# =====================================================================
suppressPackageStartupMessages({
  library(data.table); library(ggplot2); library(ggforce)
})
source(file.path(Sys.getenv("ECV_ROOT", "/Users/sindrle/Research/Prosjekter/ECV"),
                 "review_round_1/01_scripts/R1_config.R"))
ROOT <- ECV_ROOT
OUT  <- file.path(ROOT, "review_round_1/02_output"); setwd(ROOT)
S <- readRDS(file.path(OUT, "R1_00_normalised.rds"))
source(file.path(ROOT, "review_round_1/01_scripts/R1_theme.R"))

FDR <- 0.05

# --- 1. les inn de to komplette resultattabellene ------------------
# Begge filene inneholder ALLE testede proteiner, ikke bare de
# signifikante — nødvendig for å definere universet i punkt 4.
CHR <- fread(file.path(OUT, "R1_01_longterm_primary.csv"))   # A1 -> B1, 12 uker
ACU <- fread(file.path(OUT, "R1_05_acute.csv"))              # B1 -> B2 -> B3, akutt

CHR <- CHR[!is.na(Genes) & Genes != ""]
ACU <- ACU[!is.na(Genes) & Genes != ""]

# Én rad per gen. I R1_01 er Protein.Group og Genes 1:1 i praksis, men
# vi sikrer oss: behold raden med lavest q per gen.
setorder(CHR, q); CHR1 <- unique(CHR, by = "Genes")
setorder(ACU, q_time); ACU1 <- unique(ACU, by = "Genes")

tested_chr <- CHR1$Genes
tested_acu <- ACU1$Genes

# --- 2. universet: testbart i BEGGE analyser -----------------------
# IKKE alle 2971 proteingrupper. Completeness-filteret (>=70 % innen
# hvert tidspunkt) og lmer-konvergens fjerner ulike proteiner i de to
# analysene, så universet er snittet av de to testede settene.
UNIV <- intersect(tested_chr, tested_acu)
n_univ <- length(UNIV)

# --- 3. responderlistene ------------------------------------------
# To varianter: "alle" (som listene faktisk er) og "univers-begrenset"
# (bare gener som var testbare i begge). Den siste er den som hører
# sammen med den hypergeometriske testen, fordi et gen som aldri ble
# testet akutt umulig kan havne i snittet.
C_all <- CHR1[q      < FDR]$Genes
A_all <- ACU1[q_time < FDR]$Genes
C_u   <- intersect(C_all, UNIV)
A_u   <- intersect(A_all, UNIV)

ov_all <- intersect(C_all, A_all)   # identisk med ov_u per konstruksjon
ov_u   <- intersect(C_u,   A_u)
un_all <- union(C_all, A_all)
un_u   <- union(C_u,   A_u)

# --- 4. alle fire brøker, side om side ----------------------------
mk <- function(scope, metric, num, den, num_lab, den_lab, note) {
  data.table(scope = scope, metric = metric,
             numerator = num, denominator = den,
             numerator_what = num_lab, denominator_what = den_lab,
             value = num / den, percent = round(100 * num / den, 1),
             note = note)
}
RES <- rbindlist(list(
  mk("all_significant", "overlap / chronic", length(ov_all), length(C_all),
     "genes significant in BOTH", "genes significant for 12-week training",
     "Answers: what fraction of long-term responders also respond acutely."),
  mk("all_significant", "overlap / acute", length(ov_all), length(A_all),
     "genes significant in BOTH", "genes significant for the acute bout",
     "Answers: what fraction of acute responders also change over 12 weeks."),
  mk("all_significant", "overlap / union", length(ov_all), length(un_all),
     "genes significant in BOTH", "genes significant in EITHER analysis",
     "THIS IS THE DENOMINATOR USED IN THE SUBMITTED MANUSCRIPT (reported as 21.6%). It is not a proportion of responders; it equals the Jaccard index."),
  mk("all_significant", "Jaccard", length(ov_all), length(un_all),
     "intersection", "union",
     "Identical to overlap/union by definition; reported separately only because the manuscript labelled it as an overlap percentage."),
  mk("common_universe", "overlap / chronic", length(ov_u), length(C_u),
     "genes significant in BOTH", "chronic responders testable in both analyses",
     "PRIMARY. Denominator excludes chronic responders that were not testable in the acute model and therefore could never enter the intersection."),
  mk("common_universe", "overlap / acute", length(ov_u), length(A_u),
     "genes significant in BOTH", "acute responders testable in both analyses",
     "Reciprocal view of the same intersection."),
  mk("common_universe", "overlap / union", length(ov_u), length(un_u),
     "genes significant in BOTH", "genes significant in EITHER (universe-restricted)",
     "Shown only to document what the erroneous denominator would be after re-analysis."),
  mk("common_universe", "Jaccard", length(ov_u), length(un_u),
     "intersection", "union", "Identical to overlap/union.")
))

# --- 5. hypergeometrisk test mot riktig univers --------------------
# Urne: n_univ gener testbare i begge. "Hvite kuler" = kroniske
# respondere i universet. Trekk = akutte respondere i universet.
# P = Pr(X >= observert snitt).
q_hit <- length(ov_u); m_wh <- length(C_u); n_bl <- n_univ - m_wh; k_dr <- length(A_u)
p_hyper  <- phyper(q_hit - 1, m_wh, n_bl, k_dr, lower.tail = FALSE)
expected <- k_dr * m_wh / n_univ
fold_enr <- q_hit / expected
# Fishers eksakte på samme 2x2 gir samme ensidige p — tas med som kontroll.
tab <- matrix(c(q_hit, k_dr - q_hit, m_wh - q_hit,
                n_univ - m_wh - k_dr + q_hit), nrow = 2)
ft  <- fisher.test(tab, alternative = "greater")

cat(sprintf(paste0("\n=== UNIVERS ===\ntestet kronisk %d | testet akutt %d | ",
                   "testbart i BEGGE %d\n"),
            length(tested_chr), length(tested_acu), n_univ))
cat(sprintf("kroniske respondere %d (%d i universet) | akutte %d (%d i universet)\n",
            length(C_all), m_wh, length(A_all), k_dr))
cat(sprintf("\n=== HYPERGEOMETRISK ===\nsnitt %d | forventet %.2f | fold enrichment %.2f | p = %.3g (Fisher OR %.2f, p = %.3g)\n",
            q_hit, expected, fold_enr, p_hyper, unname(ft$estimate), ft$p.value))
cat("\n  ⚠ P-VERDIEN ER ANTIKONSERVATIV OG SKAL IKKE SITERES.\n",
    "  Kronisk kontrast = B1-A1. Den akutte omnibus-testen gaar over B1/B2/B3.\n",
    "  B1 inngaar i BEGGE teststatistikkene, saa et protein med ekstrem B1 faar\n",
    "  oppblaast statistikk i begge analysene. De to signifikanslistene er\n",
    "  dermed ikke uavhengig trukket, og hypergeometrisk test forutsetter det.\n",
    "  Maalt avhengighet: Spearman(-log10 p_kronisk, -log10 p_akutt) = +0.39.\n\n",
    "  Overlappet rapporteres DESKRIPTIVT. Den formelle sammenligningen av\n",
    "  akutt mot kronisk gjoeres terskelfritt med split-half-designet i\n",
    "  R1_05_acute_chronic.R (disjunkte deltakere, rho = -0.41), som ikke\n",
    "  deler noen maaling mellom de to estimatene.\n")

# --- 6. retningskonkordans blant de overlappende ------------------
# Kronisk estimat = B1 - A1. Akutt: est_B2 = B2 - B1 (rett etter økt),
# est_B3 = B3 - B1 (2 t etter). B3-B2 er den eneste akutte kontrasten
# som ikke deler B1 med det kroniske estimatet.
OVD <- merge(CHR1[Genes %in% ov_u, .(Genes, chronic_est = Estimate, chronic_q = q)],
             ACU1[Genes %in% ov_u, .(Genes, acute_est_B2 = est_B2,
                                     acute_est_B3 = est_B3, acute_q = q_time)],
             by = "Genes")
OVD[, acute_est_B3_B2 := acute_est_B3 - acute_est_B2]
OVD[, `:=`(same_sign_B2    = sign(chronic_est) == sign(acute_est_B2),
           same_sign_B3    = sign(chronic_est) == sign(acute_est_B3),
           same_sign_B3_B2 = sign(chronic_est) == sign(acute_est_B3_B2))]
setorder(OVD, chronic_q)

conc <- data.table(
  contrast = c("B2-B1 (immediately post)", "B3-B1 (2 h post)", "B3-B2 (2 h vs post)"),
  n_same_sign = c(sum(OVD$same_sign_B2), sum(OVD$same_sign_B3), sum(OVD$same_sign_B3_B2)),
  n_overlap   = nrow(OVD))
conc[, percent := round(100 * n_same_sign / n_overlap, 1)]

# ---- RIKTIG NULLHYPOTESE ------------------------------------------
# En binomialtest mot 50 % er FEIL her. Fortegnssamsvar mellom kronisk og
# B3-B2 er en egenskap ved hele proteomet (r = +0.56 over alle testede
# gener), ikke ved de overlappende. Mot 50 % blir p feil med ~3 størrelses-
# ordener. Nullen må være den EMPIRISKE bakgrunnsraten, målt blant kroniske
# respondere som IKKE er akutt-signifikante.
CHR1 <- unique(CHR, by = "Genes")
both <- merge(CHR1[, .(Genes, chronic_est = Estimate, chronic_q = q)],
              ACU1[, .(Genes, est_B2, est_B3, q_time)], by = "Genes")
both[, est_B3_minus_B2 := est_B3 - est_B2]
ref <- both[chronic_q < FDR & q_time >= FDR]     # kronisk responder, ikke akutt
for (i in seq_len(nrow(conc))) {
  col <- c("est_B2", "est_B3", "est_B3_minus_B2")[i]
  rate <- mean(sign(ref$chronic_est) == sign(ref[[col]]))
  conc[i, `:=`(bg_rate = round(rate, 3), bg_n = nrow(ref),
               p_vs_bg = binom.test(n_same_sign, n_overlap, rate, "greater")$p.value,
               p_binom_50 = binom.test(n_same_sign, n_overlap, 0.5, "greater")$p.value)]
}
cat("\nMERK: p_binom_50 beholdes kun for å vise hvor misvisende den er.\n",
    "Tallet som gjelder er p_vs_bg — mot den empiriske bakgrunnsraten.\n")

cat("\n=== RETNINGSKONKORDANS BLANT DE OVERLAPPENDE ===\n"); print(conc)
cat("\n=== OVERLAPPENDE GENER ===\n")
print(OVD[, .(Genes, chronic = round(chronic_est, 3), B2 = round(acute_est_B2, 3),
              B3 = round(acute_est_B3, 3), same_B2 = same_sign_B2)])

# --- 7. skriv ut ---------------------------------------------------
# Nøkkeltallene legges nederst i samme fil som brøkene, slik at hele
# regnestykket kan siteres fra ett sted.
EXTRA <- rbindlist(list(
  data.table(scope = "universe", metric = "tested_chronic", numerator = length(tested_chr),
             denominator = NA_integer_, numerator_what = "genes tested in the 12-week model",
             denominator_what = NA_character_, value = NA_real_, percent = NA_real_,
             note = "R1_01_longterm_primary.csv, all rows"),
  data.table(scope = "universe", metric = "tested_acute", numerator = length(tested_acu),
             denominator = NA_integer_, numerator_what = "genes tested in the acute model",
             denominator_what = NA_character_, value = NA_real_, percent = NA_real_,
             note = "R1_05_acute.csv, all rows"),
  data.table(scope = "universe", metric = "tested_both", numerator = n_univ,
             denominator = NA_integer_, numerator_what = "genes testable in BOTH models",
             denominator_what = NA_character_, value = NA_real_, percent = NA_real_,
             note = "Hypergeometric universe. NOT the 2971 protein groups quantified."),
  data.table(scope = "enrichment", metric = "expected_overlap", numerator = NA_integer_,
             denominator = NA_integer_, numerator_what = "k*m/N under independence",
             denominator_what = NA_character_, value = expected, percent = NA_real_,
             note = sprintf("k=%d acute, m=%d chronic, N=%d", k_dr, m_wh, n_univ)),
  data.table(scope = "enrichment", metric = "fold_enrichment", numerator = q_hit,
             denominator = NA_integer_, numerator_what = "observed overlap",
             denominator_what = "expected overlap", value = fold_enr, percent = NA_real_,
             note = "observed / expected"),
  data.table(scope = "enrichment", metric = "hypergeometric_p", numerator = NA_integer_,
             denominator = NA_integer_, numerator_what = "Pr(X >= observed)",
             denominator_what = NA_character_, value = p_hyper, percent = NA_real_,
             note = sprintf("phyper(%d, %d, %d, %d, lower.tail=FALSE)",
                            q_hit - 1, m_wh, n_bl, k_dr)),
  data.table(scope = "concordance", metric = "same_sign_B2", numerator = conc$n_same_sign[1],
             denominator = conc$n_overlap[1], numerator_what = "same sign, chronic vs B2-B1",
             denominator_what = "overlapping genes", value = conc$n_same_sign[1]/conc$n_overlap[1],
             percent = conc$percent[1],
             note = sprintf("p = %.3g mot empirisk bakgrunnsrate %.3f (IKKE mot 50%%)",
                            conc$p_vs_bg[1], conc$bg_rate[1])),
  data.table(scope = "concordance", metric = "same_sign_B3", numerator = conc$n_same_sign[2],
             denominator = conc$n_overlap[2], numerator_what = "same sign, chronic vs B3-B1",
             denominator_what = "overlapping genes", value = conc$n_same_sign[2]/conc$n_overlap[2],
             percent = conc$percent[2],
             note = sprintf("p = %.3g mot empirisk bakgrunnsrate %.3f (IKKE mot 50%%)",
                            conc$p_vs_bg[2], conc$bg_rate[2])),
  data.table(scope = "concordance", metric = "same_sign_B3_B2", numerator = conc$n_same_sign[3],
             denominator = conc$n_overlap[3], numerator_what = "same sign, chronic vs B3-B2",
             denominator_what = "overlapping genes", value = conc$n_same_sign[3]/conc$n_overlap[3],
             percent = conc$percent[3],
             note = sprintf(paste0("p = %.3g mot empirisk bakgrunnsrate %.3f. B3-B2 deler ",
                                   "ingen maaling med kronisk, men maaler RECOVERY, ikke ",
                                   "akuttrespons - se R1_05 split-half (rho = -0.41)"),
                            conc$p_vs_bg[3], conc$bg_rate[3]))
))
OUTTAB <- rbind(RES, EXTRA)
fwrite(OUTTAB, file.path(OUT, "R1_36_overlap.csv"))
fwrite(OVD,    file.path(OUT, "R1_36_overlap_genes.csv"))

cat("\n=== ALLE BRØKER SIDE OM SIDE ===\n")
print(RES[, .(scope, metric, numerator, denominator, percent)])

# --- 8. Venn (to sirkler, PAL$venn) -------------------------------
# Gir merverdi fordi den viser at nevnerne 125, 35 og 148 er tre ulike
# tall i samme figur.
circ <- data.table(x = c(-0.55, 0.55), y = c(0, 0), r = c(1.15, 0.85),
                   set = c("12-week training", "Acute exercise"))
labs <- data.table(
  x = c(-1.05, 1.05, 0.12), y = c(0, 0, 0),
  lab = c(sprintf("%d", m_wh - q_hit), sprintf("%d", k_dr - q_hit), sprintf("%d", q_hit)))
hdr <- data.table(x = c(-1.05, 1.05), y = c(1.45, 1.15),
                  lab = c(sprintf("12-week training\n(n = %d)", m_wh),
                          sprintf("Acute exercise\n(n = %d)", k_dr)))
p <- ggplot() +
  geom_circle(data = circ, aes(x0 = x, y0 = y, r = r, fill = set),
              alpha = 0.45, colour = NA) +
  geom_text(data = labs, aes(x, y, label = lab), size = 4.2, family = fnt) +
  geom_text(data = hdr, aes(x, y, label = lab), size = 3.6, family = fnt,
            lineheight = 0.95) +
  annotate("text", x = 0, y = -1.55, family = fnt, size = 3.2,
           label = sprintf("%d / %d chronic responders (%.1f%%) also respond acutely\nhypergeometric p = %.3g, %.1f-fold enrichment, universe = %d proteins tested in both",
                           q_hit, m_wh, 100 * q_hit / m_wh, p_hyper, fold_enr, n_univ)) +
  scale_fill_manual(values = unname(PAL$venn)) +
  coord_fixed(clip = "off") + theme_void() +
  theme(legend.position = "none", text = element_text(family = fnt))
save_ecv("R1_36_venn_acute_chronic.pdf", p, width = 12, height = 10, unit = "cm")

cat(sprintf("\nTIL ABSTRACT: %d av %d (%.1f%%) — overlapp / kroniske respondere testbare i begge.\n",
            q_hit, m_wh, 100 * q_hit / m_wh))
cat(sprintf("IKKE bruk %d/%d = %.1f%% (overlapp/union) — det var feilen i det innsendte manuset.\n",
            q_hit, length(un_u), 100 * q_hit / length(un_u)))
