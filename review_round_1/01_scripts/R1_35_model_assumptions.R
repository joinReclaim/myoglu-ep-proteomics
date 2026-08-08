# =====================================================================
# R1_35_model_assumptions.R — modellantakelser og random slopes
# Adresserer R3: «begrunn modellspesifikasjonen; ble random slopes
# vurdert?» og «modellantakelser (normalitet, homoskedastisitet) ble
# ikke vurdert».
#
# Primærmodellen (R1_01_longterm.R) er
#     y ~ SampleCode + Gr + batch + (1 | ID2),  REML = FALSE
# på cyclic loess-normaliserte data, completeness >= 70 % innen hvert
# tidspunkt, A1 mot B1.
#
# DEL 1 — antakelser for ALLE proteiner i primæranalysen:
#   * Shapiro-Wilk på betingede residualer (og marginale, som kontroll)
#   * Breusch-Pagan (Koenker-studentisert: n*R^2 fra e^2 ~ fitted)
#   * Spearman |residual| mot fitted
#   * Levene: lik residualvarians ved A1 og B1
#   * Shapiro-Wilk på de estimerte tilfeldige konstantleddene
#
# DEL 2 — random slopes:
#     y ~ SampleCode + Gr + batch + (1 + SampleCode | ID2)
#   LRT mot primærmodellen per protein, og — viktigst — hva som skjer
#   med SampleCodeB1-estimatet og FDR-listen.
#
# MERK om identifiserbarhet: med nøyaktig to tidspunkter og én måling
# per deltaker per tidspunkt gir random slope-modellen tre varians-
# parametre (konstantledd, stigning, korrelasjon) pluss residualvarians
# fra to observasjoner per deltaker. Stigningsvariansen og residual-
# variansen er da ikke separat identifiserbare. Singulære tilpasninger
# er derfor forventet og er selve BEGRUNNELSEN for modellspesifikasjonen,
# ikke et beregningsproblem.
#
# Krever: R1_00_setup.R (cache), R1_01_longterm.R (primærliste)
# =====================================================================
suppressPackageStartupMessages({
  library(data.table); library(lme4); library(lmerTest); library(car)
  library(ggplot2); library(cowplot)
})
source(file.path(Sys.getenv("ECV_ROOT", "/Users/sindrle/Research/Prosjekter/ECV"),
                 "review_round_1/01_scripts/R1_config.R"))
ROOT <- ECV_ROOT
OUT  <- file.path(ROOT, "review_round_1/02_output"); setwd(ROOT)
S <- readRDS(file.path(OUT, "R1_00_normalised.rds"))
source(file.path(ROOT, "review_round_1/01_scripts/R1_theme.R"))
set.seed(1)

THR <- 0.70   # samme completeness-filter som primæranalysen

# ---------------------------------------------------------------------
# Datasett: nøyaktig samme utvalg og filter som R1_01_longterm.R
# ---------------------------------------------------------------------
idx <- which(S$si$SampleCode %in% c("A1", "B1"))
M   <- S$Xloe[, idx, drop = FALSE]; s <- S$si[idx]
nA  <- rowSums(!is.na(M[, s$SampleCode == "A1", drop = FALSE]))
nB  <- rowSums(!is.na(M[, s$SampleCode == "B1", drop = FALSE]))
ok  <- nA >= THR * 26 & nB >= THR * 26
M   <- M[ok, , drop = FALSE]; g <- S$genes[ok]
pg  <- rownames(M)

d0 <- data.table(SampleCode = factor(s$SampleCode, c("A1", "B1")),
                 Gr = s$Gr, batch = s$batch, ID2 = factor(s$ID2))

cat(sprintf("Deltakere: %d | injeksjoner A1+B1: %d | proteiner etter >=%.0f%%: %d\n",
            uniqueN(s$ID2), nrow(d0), 100 * THR, nrow(M)))
cat("Balanse SampleCode x batch:\n"); print(table(d0$SampleCode, d0$batch))
cat("Observasjoner per deltaker (A1/B1, før NA-fjerning per protein):\n")
print(table(table(d0$ID2)))

form_ri <- y ~ SampleCode + Gr + batch + (1 | ID2)
form_rs <- y ~ SampleCode + Gr + batch + (1 + SampleCode | ID2)

# --- Breusch-Pagan, Koenker-studentisert: n*R^2 fra lm(e^2 ~ fitted) ---
bp_test <- function(e, fit) {
  if (length(e) < 8 || sd(fit) == 0) return(c(NA_real_, NA_real_))
  u <- e^2
  r2 <- suppressWarnings(summary(lm(u ~ fit))$r.squared)
  if (!is.finite(r2)) return(c(NA_real_, NA_real_))
  st <- length(e) * r2
  c(st, pchisq(st, df = 1, lower.tail = FALSE))
}
safe_shapiro <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) < 5 || sd(x) == 0) return(NA_real_)
  tryCatch(shapiro.test(x)$p.value, error = function(e) NA_real_)
}

# ---------------------------------------------------------------------
# Hovedløkke: primærmodell + random slope-modell per protein
# ALDRI tryCatch(warning=) rundt lmer — kun error.
# ---------------------------------------------------------------------
resid_store <- vector("list", nrow(M))
rows <- vector("list", nrow(M))

t0 <- Sys.time()
for (i in seq_len(nrow(M))) {
  d <- copy(d0); d[, y := M[i, ]]; d <- d[!is.na(y)]
  if (uniqueN(d$SampleCode) < 2 || uniqueN(d$ID2) < 5) next
  if (uniqueN(d$batch) < 2) next

  f0 <- tryCatch(suppressWarnings(suppressMessages(
          lmer(form_ri, data = d, REML = FALSE))), error = function(e) NULL)
  if (is.null(f0)) next
  co0 <- summary(f0)$coefficients
  if (!"SampleCodeB1" %in% rownames(co0)) next
  a0 <- tryCatch(car::Anova(f0, type = 2), error = function(e) NULL)
  if (is.null(a0)) next

  # ---- DEL 1: antakelser på primærmodellen ----
  e_c <- as.numeric(residuals(f0))                       # betingede residualer
  fitc <- as.numeric(fitted(f0))
  e_m <- as.numeric(d$y - predict(f0, re.form = NA))     # marginale residualer
  fitm <- as.numeric(predict(f0, re.form = NA))
  re   <- as.numeric(ranef(f0)$ID2[, 1])

  bpc <- bp_test(e_c, fitc); bpm <- bp_test(e_m, fitm)
  sp  <- suppressWarnings(cor.test(abs(e_c), fitc, method = "spearman"))
  lev <- tryCatch(car::leveneTest(e_c, group = d$SampleCode,
                                  center = median)[1, "Pr(>F)"],
                  error = function(e) NA_real_)

  resid_store[[i]] <- e_c

  # ---- DEL 2: random slope ----
  # Først med lme4s standardkontroller. Med to observasjoner per deltaker og
  # to tilfeldige effekter per deltaker er antall tilfeldige effekter lik
  # antall observasjoner, og lme4 nekter å tilpasse modellen. Vi lagrer
  # feilmeldingen — den ER svaret til R3.
  f1d <- tryCatch(suppressWarnings(suppressMessages(
           lmer(form_rs, data = d, REML = FALSE))),
         error = function(e) conditionMessage(e))
  rs_default_ok <- inherits(f1d, "merMod")
  errmsg <- if (rs_default_ok) NA_character_ else gsub("\\s+", " ", f1d)

  # Deretter TVUNGET tilpasning med identifiserbarhetssjekkene slått av, for
  # å kunne kvantifisere hva som faktisk skjer (LRT, estimatendring, FDR).
  f1 <- tryCatch(suppressWarnings(suppressMessages(
          lmer(form_rs, data = d, REML = FALSE,
               control = lmerControl(check.nobs.vs.nRE = "ignore",
                                     check.nobs.vs.nlev = "ignore",
                                     check.conv.singular = "ignore")))),
        error = function(e) NULL)
  fit_rs <- !is.null(f1)
  sing <- NA; est1 <- se1 <- p1 <- ll1 <- lrt <- lrtp <- sdsl <- rcor <- NA_real_
  sig1 <- NA_real_
  conv_ok <- NA
  if (fit_rs) {
    sig1 <- sigma(f1)
    sing <- isSingular(f1, tol = 1e-4)
    # konvergenskoder fra optimereren
    conv_ok <- is.null(f1@optinfo$conv$lme4$code) ||
               length(f1@optinfo$conv$lme4$code) == 0
    co1 <- summary(f1)$coefficients
    if ("SampleCodeB1" %in% rownames(co1)) {
      est1 <- co1["SampleCodeB1", "Estimate"]
      se1  <- co1["SampleCodeB1", "Std. Error"]
      a1 <- tryCatch(car::Anova(f1, type = 2), error = function(e) NULL)
      if (!is.null(a1)) p1 <- a1["SampleCode", "Pr(>Chisq)"]
    }
    ll1 <- as.numeric(logLik(f1))
    lrt <- 2 * (ll1 - as.numeric(logLik(f0)))
    lrt <- max(lrt, 0)
    # df = 2 (stigningsvarians + korrelasjon). Testen ligger på grensen av
    # parameterrommet, så p er konservativ — det svekker ikke konklusjonen.
    lrtp <- pchisq(lrt, df = 2, lower.tail = FALSE)
    vc <- as.data.frame(VarCorr(f1))
    sdsl <- vc$sdcor[vc$grp == "ID2" & !is.na(vc$var1) &
                     vc$var1 == "SampleCodeB1" & is.na(vc$var2)][1]
    rcor <- vc$sdcor[vc$grp == "ID2" & !is.na(vc$var2)][1]
  }

  rows[[i]] <- data.table(
    Protein.Group = pg[i], Genes = g[i],
    n_obs = nrow(d), n_ID = uniqueN(d$ID2),
    est_ri = co0["SampleCodeB1", "Estimate"],
    se_ri  = co0["SampleCodeB1", "Std. Error"],
    p_ri   = a0["SampleCode", "Pr(>Chisq)"],
    logLik_ri = as.numeric(logLik(f0)),
    shapiro_p = safe_shapiro(e_c),
    shapiro_p_marginal = safe_shapiro(e_m),
    shapiro_p_ranef = safe_shapiro(re),
    bp_stat = bpc[1], bp_p = bpc[2],
    bp_p_marginal = bpm[2],
    spearman_rho_absres_fitted = as.numeric(sp$estimate),
    spearman_p = sp$p.value,
    levene_p = lev,
    rs_default_ok = rs_default_ok, rs_error = errmsg,
    rs_fitted = fit_rs, rs_converged = conv_ok, rs_singular = sing,
    est_rs = est1, se_rs = se1, p_rs = p1,
    logLik_rs = ll1, LRT_chisq = lrt, LRT_p = lrtp,
    sd_slope = sdsl, cor_int_slope = rcor,
    sigma_ri = sigma(f0), sigma_rs = sig1)
  if (i %% 250 == 0) cat(sprintf("  ... %d/%d\n", i, nrow(M)))
}
R <- rbindlist(rows)
cat(sprintf("Løkke ferdig på %.1f min | modeller: %d\n",
            as.numeric(difftime(Sys.time(), t0, units = "mins")), nrow(R)))

R[, q_ri := p.adjust(p_ri, "BH")]

# =====================================================================
# DEL 1 — rapport: antakelser
# =====================================================================
cat("\n=== DEL 1: MODELLANTAKELSER (primærmodellen) ===\n")
n <- nrow(R)
pr <- function(x) sprintf("%d (%.1f %%)", sum(x, na.rm = TRUE),
                          100 * mean(x, na.rm = TRUE))
cat(sprintf("Proteiner vurdert: %d\n", n))
cat(sprintf("Shapiro-Wilk p < 0.05 (betingede residualer): %s\n",
            pr(R$shapiro_p < 0.05)))
cat(sprintf("Shapiro-Wilk p < 0.05 etter BH over proteiner:  %s\n",
            pr(p.adjust(R$shapiro_p, "BH") < 0.05)))
cat(sprintf("Shapiro-Wilk p < 0.05 (marginale residualer):   %s\n",
            pr(R$shapiro_p_marginal < 0.05)))
cat(sprintf("Shapiro-Wilk p < 0.05 (tilfeldige konstantledd):%s\n",
            pr(R$shapiro_p_ranef < 0.05)))
cat(sprintf("Breusch-Pagan p < 0.05 (betinget):              %s\n",
            pr(R$bp_p < 0.05)))
cat(sprintf("Breusch-Pagan p < 0.05 (marginalt):             %s\n",
            pr(R$bp_p_marginal < 0.05)))
cat(sprintf("Breusch-Pagan p < 0.05 etter BH:                %s\n",
            pr(p.adjust(R$bp_p, "BH") < 0.05)))
cat(sprintf("Spearman |resid| vs fitted, p < 0.05:           %s\n",
            pr(R$spearman_p < 0.05)))
cat(sprintf("Levene (A1 vs B1) p < 0.05:                     %s\n",
            pr(R$levene_p < 0.05)))
cat(sprintf("BRUDD PÅ BEGGE (Shapiro < 0.05 OG BP < 0.05):   %s\n",
            pr(R$shapiro_p < 0.05 & R$bp_p < 0.05)))
cat(sprintf("  forventet ved uavhengighet under nullen:      %.1f %%\n",
            100 * mean(R$shapiro_p < 0.05, na.rm = TRUE) *
                  mean(R$bp_p < 0.05, na.rm = TRUE)))
cat(sprintf("Median Shapiro-p %.3f | median BP-p %.3f | median |rho| %.3f\n",
            median(R$shapiro_p, na.rm = TRUE), median(R$bp_p, na.rm = TRUE),
            median(abs(R$spearman_rho_absres_fitted), na.rm = TRUE)))

# Er bruddene konsentrert blant funnene? (avgjør om de truer konklusjonen)
sigset <- R[q_ri < 0.05]
cat(sprintf("Blant FDR<0.05-proteiner (n=%d): Shapiro<0.05 %s | BP<0.05 %s\n",
            nrow(sigset), pr(sigset$shapiro_p < 0.05), pr(sigset$bp_p < 0.05)))

# =====================================================================
# DEL 2 — rapport: random slopes
# =====================================================================
cat("\n=== DEL 2: RANDOM SLOPES (1 + SampleCode | ID2) ===\n")
cat("--- Med lme4s standardkontroller (det man ville gjort i praksis) ---\n")
cat(sprintf("Tilpasset uten feil:                   %s\n", pr(R$rs_default_ok)))
cat("Feilmelding (tall fjernet, så meldingene kan grupperes):\n")
R[, rs_error_type := gsub("[0-9]+", "N", rs_error)]
print(R[!is.na(rs_error_type), .N, by = rs_error_type][order(-N)])
cat(sprintf("\nn = %d deltakere x 2 tidspunkter = %d observasjoner; ",
            uniqueN(d0$ID2), nrow(d0)))
cat(sprintf("(1 + SampleCode | ID2) krever %d tilfeldige effekter.\n",
            2 * uniqueN(d0$ID2)))
cat("=> Antall tilfeldige effekter er lik antall observasjoner. Stignings-\n")
cat("   variansen er ikke identifiserbar sammen med residualvariansen.\n")

cat("\n--- Tvunget tilpasning (identifiserbarhetssjekker slått av) ---\n")
cat(sprintf("Tilpasset uten feil:                   %s\n", pr(R$rs_fitted)))
cat(sprintf("Optimereren rapporterte konvergens:    %s\n", pr(R$rs_converged)))
cat(sprintf("SINGULÆR tilpasning (isSingular):      %s\n", pr(R$rs_singular)))
cat(sprintf("Ikke-singulær OG konvergerte:          %s\n",
            pr(R$rs_fitted & R$rs_converged & !R$rs_singular)))
cat(sprintf("|korrelasjon konstantledd-stigning| > 0.99: %s\n",
            pr(abs(R$cor_int_slope) > 0.99)))
cat(sprintf("Estimert SD(stigning) = 0 (< 1e-4):    %s\n",
            pr(R$sd_slope < 1e-4)))
# Det tydeligste tegnet på ikke-identifiserbarhet: residualvariansen kollapser
# fordi stigningsvariansen suger den til seg. isSingular fanger ikke dette.
cat(sprintf("Residual-SD: median %.4f med (1|ID2) -> %.4f med random slope (ratio %.3f)\n",
            median(R$sigma_ri, na.rm = TRUE), median(R$sigma_rs, na.rm = TRUE),
            median(R$sigma_rs / R$sigma_ri, na.rm = TRUE)))
cat(sprintf("Residual-SD falt til under 1 %% av (1|ID2)-verdien: %s\n",
            pr(R$sigma_rs / R$sigma_ri < 0.01)))
cat(sprintf("SD(stigning) >= residual-SD i random slope-modellen: %s\n",
            pr(R$sd_slope >= R$sigma_rs)))

FIT <- R[rs_fitted == TRUE & is.finite(LRT_p) & is.finite(p_rs)]
cat(sprintf("\nLRT tilgjengelig for %d proteiner\n", nrow(FIT)))
cat(sprintf("LRT p < 0.05 (nominelt):               %s\n", pr(FIT$LRT_p < 0.05)))
cat(sprintf("LRT p < 0.05 etter BH:                 %s\n",
            pr(p.adjust(FIT$LRT_p, "BH") < 0.05)))
NS <- FIT[rs_singular == FALSE]
cat(sprintf("Blant ikke-singulære (n=%d): LRT p<0.05 %s\n",
            nrow(NS), if (nrow(NS)) pr(NS$LRT_p < 0.05) else "0 (NA)"))
cat(sprintf("Median LRT-chisq %.3f (df=2 => forventet 2 hvis modellen trengtes)\n",
            median(FIT$LRT_chisq, na.rm = TRUE)))

# ---- Effektestimat: endrer det seg? ----
if (nrow(FIT) < 3) {
  cat("\n--- SampleCodeB1-estimatet: for få tilpassede modeller til sammenligning ---\n")
  FIT[, `:=`(q_rs = NA_real_, q_ri_sub = NA_real_)]
} else {
cat("\n--- SampleCodeB1-estimatet ---\n")
cat(sprintf("Pearson  r(est_ri, est_rs) = %.5f\n",
            cor(FIT$est_ri, FIT$est_rs, use = "complete.obs")))
cat(sprintf("Spearman rho               = %.5f\n",
            cor(FIT$est_ri, FIT$est_rs, method = "spearman", use = "complete.obs")))
cat(sprintf("Median |est_rs - est_ri| = %.5f | maks = %.5f | median |est_ri| = %.4f\n",
            median(abs(FIT$est_rs - FIT$est_ri), na.rm = TRUE),
            max(abs(FIT$est_rs - FIT$est_ri), na.rm = TRUE),
            median(abs(FIT$est_ri), na.rm = TRUE)))
cat(sprintf("Samme fortegn: %.1f %%\n", 100 * mean(FIT$est_ri * FIT$est_rs > 0, na.rm = TRUE)))
cat(sprintf("Median SE-ratio (rs/ri) = %.4f\n", median(FIT$se_rs / FIT$se_ri, na.rm = TRUE)))

# ---- FDR-listen: robust? (BH beregnes på nytt innen samme proteinsett) ----
FIT[, q_rs := p.adjust(p_rs, "BH")]
FIT[, q_ri_sub := p.adjust(p_ri, "BH")]
a <- FIT[q_ri_sub < 0.05, Protein.Group]; b <- FIT[q_rs < 0.05, Protein.Group]
jac <- length(intersect(a, b)) / length(union(a, b))
cat("\n--- FDR-listen ---\n")
cat(sprintf("FDR<0.05 med (1|ID2): %d | med random slope: %d | felles: %d | Jaccard %.3f\n",
            length(a), length(b), length(intersect(a, b)), jac))
cat(sprintf("Tapt fra responderlisten: %d | nye: %d\n",
            length(setdiff(a, b)), length(setdiff(b, a))))
cat(sprintf("Spearman rho mellom p-verdiene: %.4f\n",
            cor(FIT$p_ri, FIT$p_rs, method = "spearman", use = "complete.obs")))
prim_hits <- R[q_ri < 0.05, Protein.Group]
cat(sprintf("Av de %d primære FDR-treffene beholdes %d (%.1f %%) i random slope-modellen\n",
            length(prim_hits), length(intersect(prim_hits, b)),
            100 * length(intersect(prim_hits, b)) / length(prim_hits)))
if (length(setdiff(a, b)))
  cat("Mistet:", paste(head(FIT[Protein.Group %in% setdiff(a, b), Genes], 20), collapse = ", "), "\n")
if (length(setdiff(b, a)))
  cat("Nye:   ", paste(head(FIT[Protein.Group %in% setdiff(b, a), Genes], 20), collapse = ", "), "\n")
}

# =====================================================================
# FIGUR — engelsk tekst, PAL-farger
# =====================================================================
# Representativt utvalg: 12 proteiner jevnt fordelt over Shapiro-p-rangen
ord <- R[is.finite(shapiro_p)][order(shapiro_p)]
pick <- ord[round(seq(1, nrow(ord), length.out = 12)), .(Protein.Group, Genes, shapiro_p)]
qq <- rbindlist(lapply(seq_len(nrow(pick)), function(k) {
  e <- resid_store[[match(pick$Protein.Group[k], pg)]]
  e <- e[is.finite(e)]; z <- (e - mean(e)) / sd(e)
  data.table(lab = sprintf("%s (SW p = %s)", pick$Genes[k], signif(pick$shapiro_p[k], 2)),
             ordr = k, theo = qnorm(ppoints(length(z)))[order(order(z))], samp = z)
}))
qq[, lab := factor(lab, levels = unique(qq[order(ordr)]$lab))]

pA <- ggplot(qq, aes(theo, samp)) +
  geom_abline(slope = 1, intercept = 0, colour = PAL$venn[2], linewidth = 0.5) +
  geom_point(size = 0.6, colour = "black", alpha = 0.7) +
  facet_wrap(~ lab, ncol = 4, scales = "fixed") +
  labs(x = "Theoretical quantiles", y = "Standardised conditional residuals",
       title = "A  Normal Q-Q plots, 12 proteins spanning the Shapiro-Wilk range") +
  theme_ecv(9) + theme(plot.title = element_text(hjust = 0, size = 10),
                       strip.text = element_text(size = 7))

hist_pan <- function(x, xlab, ttl) {
  dd <- data.table(p = x[is.finite(x)])
  dd[, sig := factor(ifelse(p < 0.05, "p < 0.05", "p >= 0.05"),
                     levels = c("p < 0.05", "p >= 0.05"))]
  ggplot(dd, aes(p, fill = sig)) +
    geom_histogram(breaks = seq(0, 1, 0.05), colour = "white", linewidth = 0.2) +
    geom_hline(yintercept = nrow(dd) / 20, linetype = 2, linewidth = 0.4) +
    scale_fill_manual(values = c("p < 0.05" = unname(PAL$venn[2]),
                                 "p >= 0.05" = unname(PAL$venn[1]))) +
    scale_x_continuous(expand = expansion(mult = c(0.01, 0.01))) +
    labs(x = xlab, y = "Proteins", title = ttl) +
    theme_ecv(9) + theme(plot.title = element_text(hjust = 0, size = 10),
                         legend.position = "inside",
                         legend.position.inside = c(0.75, 0.85),
                         legend.key.size = unit(0.35, "cm"),
                         legend.text = element_text(size = 7))
}
pB <- hist_pan(R$shapiro_p, "Shapiro-Wilk p (residual normality)",
               sprintf("B  Normality  (%.1f %% below 0.05)",
                       100 * mean(R$shapiro_p < 0.05, na.rm = TRUE)))
pC <- hist_pan(R$bp_p, "Breusch-Pagan p (homoscedasticity)",
               sprintf("C  Homoscedasticity  (%.1f %% below 0.05)",
                       100 * mean(R$bp_p < 0.05, na.rm = TRUE)))

fig <- plot_grid(pA, plot_grid(pB, pC, nrow = 1, rel_widths = c(1, 1)),
                 ncol = 1, rel_heights = c(1.5, 1))
save_ecv("R1_35_model_assumptions.pdf", fig, width = 20, height = 22, unit = "cm")

# =====================================================================
# UTFILER
# =====================================================================
setcolorder(R, c("Protein.Group", "Genes", "n_obs", "n_ID"))
fwrite(R[, .(Protein.Group, Genes, n_obs, n_ID, est_ri, se_ri, p_ri, q_ri,
             shapiro_p, shapiro_p_marginal, shapiro_p_ranef,
             bp_stat, bp_p, bp_p_marginal,
             spearman_rho_absres_fitted, spearman_p, levene_p)][order(shapiro_p)],
       file.path(OUT, "R1_35_assumptions.csv"))
fwrite(merge(R[, .(Protein.Group, Genes, n_obs, n_ID, est_ri, se_ri, p_ri, q_ri,
                   logLik_ri, rs_default_ok, rs_error,
                   rs_fitted, rs_converged, rs_singular,
                   est_rs, se_rs, p_rs, logLik_rs, LRT_chisq, LRT_p,
                   sd_slope, cor_int_slope, sigma_ri, sigma_rs)],
             FIT[, .(Protein.Group, q_rs, q_ri_sub)], by = "Protein.Group",
             all.x = TRUE, sort = FALSE)[order(p_ri)],
       file.path(OUT, "R1_35_random_slopes.csv"))
cat("\nSkrevet: R1_35_assumptions.csv, R1_35_random_slopes.csv\n")
