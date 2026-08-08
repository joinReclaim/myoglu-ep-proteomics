# =====================================================================
# R1_34_participants.R — deltakerkarakteristikk og complete-case n
# Adresserer R5 punkt 9: «Please add a brief participant characteristics
# table and per-analysis complete-case numbers so the paper stands alone.»
#
# Del 1: karakteristikktabell fra S$pheno ved A1 (baseline) og B1 (etter
#        12 uker), delt på gruppe, med n / gjennomsnitt / SD og
#        gruppeforskjell. I tillegg rad «delta» (B1 - A1 per deltaker),
#        der gruppeforskjellen ER group x time-interaksjonen.
#
# Del 2: complete-case n for hver analyse i revisjonen. Tallene HENTES FRA
#        de lagrede resultatfilene i 02_output/ der det er mulig, slik at
#        tabellen beskriver det som faktisk ble kjørt — ikke en ny kjøring
#        som kan avvike. Der et tall ikke finnes lagret (f.eks. antall
#        proteinmålinger som gikk inn i modellene) rekonstrueres FILTERET
#        eksakt som i originalskriptet og telles på nytt.
#
# NB: ingen figur. Ingen eksisterende fil røres.
# Krever: R1_00_setup.R og alle R1_01..R1_17
# =====================================================================
suppressPackageStartupMessages({library(data.table); library(stringr)})
source(file.path(Sys.getenv("ECV_ROOT", "/Users/sindrle/Research/Prosjekter/ECV"),
                 "review_round_1/01_scripts/R1_config.R"))
ROOT <- ECV_ROOT
OUT  <- file.path(ROOT, "review_round_1/02_output"); setwd(ROOT)
S <- readRDS(file.path(OUT, "R1_00_normalised.rds"))

N_PART <- 26      # deltakere med MS-data
THR    <- 0.70    # completeness-terskel, som i R1_01

# =====================================================================
# DEL 1 — DELTAKERKARAKTERISTIKK
# =====================================================================
VARS <- c("Age","Weight","BMI","FFM","M_mg.kg.min","VO2max.kg","ChestPress",
          "PullDown","LegPress","mri_VL","mri_IP","mri_TAT","mri_SAT",
          "mrs_Hep","mrs_Musc")

ph <- as.data.table(S$pheno)
have    <- VARS[VARS %in% names(ph)]
mangler <- setdiff(VARS, names(ph))
cat("=== DEL 1: DELTAKERKARAKTERISTIKK ===\n")
cat(sprintf("variabler funnet i pheno : %d av %d (%s)\n", length(have), length(VARS),
            paste(have, collapse = ", ")))
if (length(mangler)) cat(sprintf("IKKE i pheno og utelatt  : %s\n", paste(mangler, collapse = ", ")))

# begrens til de 26 deltakerne som faktisk har MS-data — tabellen skal
# beskrive analysepopulasjonen, ikke hele MyoGlu-kohorten
ms_ids <- sort(unique(S$si$ID2))
ph <- ph[ID2 %in% ms_ids]
cat(sprintf("deltakere i pheno med MS-data: %d (av %d unike i pheno-fila)\n",
            uniqueN(ph$ID2), uniqueN(as.data.table(S$pheno)$ID2)))
cat("gruppefordeling per tidspunkt:\n"); print(table(ph$timepoint, ph$group))

GRP <- c("dysglyc", "control")
pa <- ph[timepoint == "A1"][order(ID2)]
pb <- ph[timepoint == "B1"][order(ID2)]
ids <- intersect(pa$ID2, pb$ID2)
pa <- pa[match(ids, ID2)]; pb <- pb[match(ids, ID2)]
stopifnot(identical(pa$group, pb$group))

# hjelpere — alt NA-sikkert, og testene faller til NA hvis n er for liten
msd  <- function(x) c(n = sum(is.finite(x)), mean = mean(x, na.rm = TRUE),
                      sd = sd(x, na.rm = TRUE))
p_t  <- function(a, b) tryCatch(t.test(a, b)$p.value, error = function(e) NA_real_)
p_w  <- function(a, b) tryCatch(suppressWarnings(wilcox.test(a, b)$p.value),
                                error = function(e) NA_real_)

row_for <- function(v, tp) {
  d <- switch(tp, "A1" = pa[[v]], "B1" = pb[[v]], "delta" = pb[[v]] - pa[[v]])
  g <- pa$group
  a <- d[g == "dysglyc"]; b <- d[g == "control"]
  sa <- msd(a); sb <- msd(b)
  data.table(variable = v, timepoint = tp,
             n_dysglyc = sa[["n"]], mean_dysglyc = sa[["mean"]], sd_dysglyc = sa[["sd"]],
             n_control = sb[["n"]], mean_control = sb[["mean"]], sd_control = sb[["sd"]],
             diff_dysglyc_minus_control = sa[["mean"]] - sb[["mean"]],
             p_group_welch    = if (sa[["n"]] >= 3 && sb[["n"]] >= 3) p_t(a, b) else NA_real_,
             p_group_wilcoxon = if (sa[["n"]] >= 3 && sb[["n"]] >= 3) p_w(a, b) else NA_real_)
}

# Age er konstant over de 12 ukene — kun baseline-raden er meningsfull
tab <- rbindlist(lapply(have, function(v) {
  tps <- if (v == "Age") "A1" else c("A1", "B1", "delta")
  rbindlist(lapply(tps, function(tp) row_for(v, tp)))
}))

# egen rad øverst med antall deltakere
hdr <- data.table(variable = "n participants", timepoint = "A1",
                  n_dysglyc = sum(pa$group == "dysglyc"), mean_dysglyc = NA_real_,
                  sd_dysglyc = NA_real_, n_control = sum(pa$group == "control"),
                  mean_control = NA_real_, sd_control = NA_real_,
                  diff_dysglyc_minus_control = NA_real_,
                  p_group_welch = NA_real_, p_group_wilcoxon = NA_real_)
tab <- rbind(hdr, tab)

# p-verdier MÅ holdes utenfor round(, 4) — flere gruppeforskjeller er ~1e-06 og
# ble skrevet som eksakt 0 i tabellen. «p = 0» er ikke publiserbart.
pcols <- intersect(c("p_group_welch", "p_group_wilcoxon"), names(tab))
num <- setdiff(names(tab), c("variable", "timepoint", pcols))
tab[, (num)   := lapply(.SD, function(x) round(as.numeric(x), 4)),  .SDcols = num]
tab[, (pcols) := lapply(.SD, function(x) signif(as.numeric(x), 3)), .SDcols = pcols]
fwrite(tab, file.path(OUT, "R1_34_participant_characteristics.csv"))

cat("\n--- baseline (A1) ---\n")
print(tab[timepoint == "A1", .(variable,
      dysglyc = ifelse(is.na(mean_dysglyc), sprintf("n=%d", n_dysglyc),
                       sprintf("%.1f (%.1f) n=%d", mean_dysglyc, sd_dysglyc, n_dysglyc)),
      control = ifelse(is.na(mean_control), sprintf("n=%d", n_control),
                       sprintf("%.1f (%.1f) n=%d", mean_control, sd_control, n_control)),
      p = signif(p_group_welch, 3))])
cat("\n--- etter 12 uker (B1) ---\n")
print(tab[timepoint == "B1", .(variable,
      dysglyc = sprintf("%.1f (%.1f) n=%d", mean_dysglyc, sd_dysglyc, n_dysglyc),
      control = sprintf("%.1f (%.1f) n=%d", mean_control, sd_control, n_control),
      p = signif(p_group_welch, 3))])
cat("\n--- endring B1-A1 (gruppeforskjell = group x time) ---\n")
print(tab[timepoint == "delta", .(variable,
      dysglyc = sprintf("%+.2f (%.2f) n=%d", mean_dysglyc, sd_dysglyc, n_dysglyc),
      control = sprintf("%+.2f (%.2f) n=%d", mean_control, sd_control, n_control),
      p = signif(p_group_welch, 3))])
cat(sprintf("\nskrevet: %s\n", file.path(OUT, "R1_34_participant_characteristics.csv")))

# =====================================================================
# DEL 2 — COMPLETE-CASE N PER ANALYSE
# =====================================================================
cat("\n\n=== DEL 2: COMPLETE-CASE N PER ANALYSE ===\n")
rows <- list()
add <- function(analysis, n_participants, n_proteins, n_observations, note)
  rows[[length(rows) + 1]] <<- data.table(analysis = analysis,
    n_participants = n_participants, n_proteins = n_proteins,
    n_observations = n_observations, note = note)

# ---- 0. kohorten som gikk inn i MS-kjøringen ------------------------
inj <- table(S$si$SampleCode)
cat(sprintf("injeksjoner per tidspunkt: %s | deltakere med alle fire: %d\n",
            paste(sprintf("%s %d", names(inj), as.integer(inj)), collapse = ", "),
            sum(rowSums(table(S$si$ID2, S$si$SampleCode) > 0) == 4)))
add("Cohort (LC-MS/MS, all injections)", uniqueN(S$si$ID2), nrow(S$X), ncol(S$X),
    sprintf("%d participants x 4 sampling points (A1, B1, B2, B3); all 26 complete at the injection level; %d protein groups with a gene symbol; %.1f%% missing values at the protein x injection level, no imputation",
            uniqueN(S$si$ID2), nrow(S$X), 100 * mean(is.na(S$X))))

# ---- R1_01: 12-ukers respondere -------------------------------------
# filteret rekonstrueres eksakt som i R1_01 for å telle målingene
i01 <- which(S$si$SampleCode %in% c("A1", "B1"))
M01 <- S$Xloe[, i01]; s01 <- S$si[i01]
nA <- rowSums(!is.na(M01[, s01$SampleCode == "A1", drop = FALSE]))
nB <- rowSums(!is.na(M01[, s01$SampleCode == "B1", drop = FALSE]))
ok01 <- nA >= THR * N_PART & nB >= THR * N_PART
P01 <- readRDS(file.path(OUT, "R1_01_longterm.rds"))$PRIM
add("R1_01 12-week response (A1 vs B1, lmer)", N_PART, nrow(P01), sum(!is.na(M01[ok01, ])),
    sprintf("%d protein groups passed the >=70%% completeness filter within each time point and %d yielded a converged model; per-protein n ranged %d-%d participants (median %d); %d significant at FDR<0.05; observations counted as non-missing protein x injection quantifications across the %d A1/B1 injections",
            sum(ok01), nrow(P01), min(P01$n_ID), max(P01$n_ID),
            as.integer(median(P01$n_ID)), sum(P01$q < 0.05), length(i01)))

# ---- R1_02: fenotypematrise -----------------------------------------
P02 <- readRDS(file.path(OUT, "R1_02_phenotypes.rds"))
r02 <- P02$res
add("R1_02 protein x phenotype matrix (Spearman)", ncol(P02$dp), nrow(P02$dp), nrow(r02),
    sprintf("delta protein (B1-A1) for %d participants x %d proteins correlated against %d phenotype deltas; %d of the %d possible protein x phenotype tests were evaluable (>=10 complete pairs); per-test n %d-%d (median %d); BH-FDR applied across the whole matrix",
            ncol(P02$dp), nrow(P02$dp), length(P02$PHENOS), nrow(r02),
            nrow(P02$dp) * length(P02$PHENOS), min(r02$n), max(r02$n),
            as.integer(median(r02$n))))

# ---- R1_03: pathways -------------------------------------------------
univ03 <- uniqueN(P01$Genes); sig03 <- uniqueN(P01[q < 0.05, Genes])
add("R1_03 over-representation analysis", N_PART, univ03, NA_integer_,
    sprintf("gene-set level analysis, no participant-level data; universe = the %d gene symbols actually tested in R1_01 (1522 mapped to Entrez), query = the %d responders at FDR<0.05 (%d up / %d down)",
            univ03, sig03, uniqueN(P01[q < 0.05 & Estimate > 0, Genes]),
            uniqueN(P01[q < 0.05 & Estimate < 0, Genes])))

# ---- R1_04: group x time --------------------------------------------
r04 <- fread(file.path(OUT, "R1_04_group_time.csv"))
add("R1_04 group x time interaction (lmer)", N_PART, nrow(r04), sum(!is.na(M01[ok01, ])),
    sprintf("13 dysglycaemic vs 13 control, same A1/B1 samples and same completeness filter as R1_01; %d proteins yielded both nested models; per-protein n %d-%d participants (median %d); %d significant at FDR<0.05",
            nrow(r04), min(r04$n_ID), max(r04$n_ID), as.integer(median(r04$n_ID)),
            sum(r04$q_int < 0.05)))

# ---- R1_05: akutt ----------------------------------------------------
i05 <- which(S$si$SampleCode %in% c("B1", "B2", "B3"))
M05 <- S$Xloe[, i05]; s05 <- S$si[i05]
ok05 <- Reduce(`&`, lapply(c("B1", "B2", "B3"), function(t)
  rowSums(!is.na(M05[, s05$SampleCode == t, drop = FALSE])) >= THR * N_PART))
r05 <- fread(file.path(OUT, "R1_05_acute.csv"))
add("R1_05 acute response (B1/B2/B3, lmer)", N_PART, nrow(r05), sum(!is.na(M05[ok05, ])),
    sprintf("%d protein groups passed >=70%% completeness at all three post-training sampling points and %d yielded a converged model; per-protein n %d-%d participants (median %d); %d significant at FDR<0.05; observations counted across the %d B1/B2/B3 injections",
            sum(ok05), nrow(r05), min(r05$n_ID), max(r05$n_ID),
            as.integer(median(r05$n_ID)), sum(r05$q_time < 0.05), length(i05)))

# ---- R1_11: serum/EP -------------------------------------------------
r11 <- fread(file.path(OUT, "R1_11_serum_ep_shared.csv"))
n_ol <- tryCatch({
  ol <- suppressWarnings(fread(file.path(OLINK_DIR, "MyoGlu_Olink_A1B1_NPX.csv"), sep = ";"))
  uniqueN(ol$ID) }, error = function(e) NA_integer_)
add("R1_11 EP vs serum (Olink) overlap", N_PART, nrow(r11), nrow(r11),
    sprintf("comparison is at the protein level between two summary tables from the same %d participants (serum Olink n = %s participants); %d gene symbols quantified on BOTH platforms form the denominator; %d of these were EP-significant and %d serum-significant at FDR<0.05",
            N_PART, ifelse(is.na(n_ol), "NA", n_ol), nrow(r11),
            sum(r11$q_ep < 0.05), sum(r11$q_se < 0.05)))

# ---- R1_12: mRNA-konkordans -----------------------------------------
r12 <- fread(file.path(OUT, "R1_12_mrna_protein_concordance.csv"))
n_skm <- NA_integer_; n_at <- NA_integer_
if (file.exists(RNASEQ_RDA)) {
  load(RNASEQ_RDA)
  # samme ID-konvensjon som R1_12 / 202_MixOmics_DIABLO.R
  paired_n <- function(mat) {
    sid <- rownames(mat)
    tm  <- str_extract(sid, "A1|B1"); raw <- as.numeric(str_extract(sid, "\\d+"))
    id  <- ifelse(tm == "A1" & raw >= 100, raw - 100,
           ifelse(tm == "B1" & raw >= 200, raw - 200, raw))
    keep <- !is.na(tm) & !is.na(id) & rowSums(!is.na(mat)) > 0
    tb <- table(id[keep], tm[keep])
    sum(rowSums(tb > 0) == 2 & rownames(tb) %in% as.character(ms_ids))
  }
  n_skm <- paired_n(fpkm.sm.fix); n_at <- paired_n(fpkm.at.fix)
}
add("R1_12 delta protein vs delta mRNA concordance", min(n_skm, n_at, na.rm = TRUE),
    nrow(r12), nrow(r12),
    sprintf("%d genes had both a filtered EP protein delta and a measured transcript; skeletal muscle: %d genes, %s participants with paired A1/B1 RNA-seq; adipose tissue: %d genes, %s participants (2 participants lack adipose samples entirely and were not imputed); per-gene correlations required >=10 complete pairs",
            nrow(r12), sum(!is.na(r12$rho_SkM)), ifelse(is.na(n_skm), "NA", n_skm),
            sum(!is.na(r12$rho_AT)), ifelse(is.na(n_at), "NA", n_at)))

# ---- R1_13: DIABLO ---------------------------------------------------
D13 <- readRDS(file.path(OUT, "R1_13_diablo.rds"))
n_comp_ep <- sum(rowSums(!is.na(M01)) == ncol(M01))   # komplette rader i EP-blokken
add("R1_13 DIABLO (4-block integration)", uniqueN(D13$s$ID2), ncol(D13$X$EP), length(D13$key),
    sprintf("%d of 26 participants had complete EP, serum, muscle AND adipose data at both A1 and B1, giving %d samples (%s); 2 participants were dropped because their adipose samples are entirely missing and were NOT imputed; each block reduced to the top %d features by variance (EP pool = %d protein groups with no missing values across all A1/B1 injections)",
            uniqueN(D13$s$ID2), length(D13$key),
            paste(sprintf("%s %d", names(table(D13$s$SampleCode)),
                          as.integer(table(D13$s$SampleCode))), collapse = ", "),
            ncol(D13$X$EP), n_comp_ep))

# ---- R1_15: partiell ------------------------------------------------
r15 <- fread(file.path(OUT, "R1_15_partial_phenotype.csv"))
add("R1_15 mutually adjusted phenotype associations", ncol(P02$dp), nrow(r15), ncol(P02$dp),
    sprintf("rank-based partial correlations on %d participant-level protein deltas; restricted to the %d proteins with q<0.25 in R1_01 that also survived the R1_02 delta filter; each partial model required >=12 complete cases",
            ncol(P02$dp), nrow(r15)))

# ---- R1_16: akutte klynger ------------------------------------------
r16 <- fread(file.path(OUT, "R1_16_acute_clusters.csv"))
ax16 <- fread(file.path(OUT, "R1_16_acute_axis_vs_gain.csv"))
add("R1_16 acute response clusters", N_PART, nrow(r16),
    ax16[utfall == "VO2max.kg", n][1],
    sprintf("%d acute responders (q<0.25 in R1_05) clustered on z-scored mean trajectories across B1/B2/B3; cluster scores computed per participant and related to training gain in %d participants with complete phenotype data",
            nrow(r16), ax16[utfall == "VO2max.kg", n][1]))

# ---- R1_17: kroniske klynger ----------------------------------------
r17 <- fread(file.path(OUT, "R1_17_chronic_clusters.csv"))
add("R1_17 chronic response clusters", ncol(P02$dp), nrow(r17), ncol(P02$dp),
    sprintf("%d chronic responders (q<0.25 in R1_01) with >=20 of %d participants having a finite delta; correlation clustering on participant-level deltas, shown both centred and uncentred",
            nrow(r17), ncol(P02$dp)))

cc <- rbindlist(rows)
fwrite(cc, file.path(OUT, "R1_34_complete_case_n.csv"))
print(cc[, .(analysis, n_participants, n_proteins, n_observations)])
cat(sprintf("\nskrevet: %s\n", file.path(OUT, "R1_34_complete_case_n.csv")))
