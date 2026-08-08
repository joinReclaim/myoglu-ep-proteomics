# =====================================================================
# R1_25_fig5CD.R — manuskriptets FIGUR 5, PANEL B og C (Figures/Figure3.tif)
#   (tidligere panel C og D — kompositten er tatt ut, panelene rykker opp)
# Erstatter 005_ECV_phenotype_correlations.R:320-425 og 201_Olink_vs_ECV.R
#
# PANEL B (tidl. C), endret:
#   - «EP enrichment»-formene (ΔEP − Δserum) er FJERNET. Størrelsen er ikke
#     en log-ratio (R5), og ble brukt som prediktor i lm(). Kun Serum og EP.
#   - originalen beregnet PEARSON men merket aksen «Spearman's rho»
#     (005:355 mot 005:419). Rettet til Spearman.
#   - ECV -> EP i legenden (Q- to be fixed).
# PANEL C (tidl. D), endret:
#   - Venn begrenset til de 505 genene som faktisk er MÅLBARE PÅ BEGGE
#     plattformer. Uten det har prosenten ingen tolkbar nevner (R5 punkt 2).
# Uendret: ggvenn med fill_color c("#4C78A8","#F58518"), 3.5 x 3.5 tommer;
#   panel C som dotplot med form per kompartment, 16 x 8 cm.
# Krever: R1_00_setup.R, R1_01_longterm.R, R1_02_phenotypes.R
# =====================================================================
suppressPackageStartupMessages({library(data.table); library(ggvenn)})
source(file.path(Sys.getenv("ECV_ROOT", "/Users/sindrle/Research/Prosjekter/ECV"),
                 "review_round_1/01_scripts/R1_config.R"))
ROOT <- ECV_ROOT
OUT  <- file.path(ROOT,"review_round_1/02_output"); setwd(ROOT)
source(file.path(ROOT,"review_round_1/01_scripts/R1_theme.R"))
PH <- readRDS(file.path(OUT,"R1_02_phenotypes.rds"))

# ---- serum (Olink) for MEP1B og HSP90B1 ----
ol <- fread(file.path(OLINK_DIR,"MyoGlu_Olink_A1B1_NPX.csv"), sep = ";")
an <- fread(file.path(OLINK_DIR,"Annotation.csv"), sep = ";", quote = "")
uid <- c(MEP1B = "Q16820", HSP90B1 = "P14625")
oid <- sapply(uid, function(u) an[UniProtID == u, OlinkID][1])
cat("Olink-ID:", paste(names(oid), oid, sep = "="), "\n")
sw <- dcast(melt(ol[, c("ID","Time", oid), with = FALSE], id.vars = c("ID","Time")),
            ID ~ variable + Time, value.var = "value")
ids <- PH$dph$ID2
ser <- sapply(names(oid), function(g)
  sw[[paste0(oid[g], "_B1")]][match(ids, sw$ID)] - sw[[paste0(oid[g], "_A1")]][match(ids, sw$ID)])

d <- rbindlist(lapply(names(oid), function(g) {
  ep <- if (g %in% rownames(PH$dp)) PH$dp[g, ] else rep(NA_real_, length(ids))
  rbindlist(list(
    data.table(protein = g, compartment = "Serum", x = ser[, g]),
    data.table(protein = g, compartment = "EP",    x = ep)))[
    , { ok <- is.finite(x) & is.finite(PH$dph$VO2max.kg)
        ct <- suppressWarnings(cor.test(x[ok], PH$dph$VO2max.kg[ok], method = "spearman"))
        .(rho = unname(ct$estimate), p = ct$p.value, n = sum(ok)) },
    by = .(protein, compartment)] }))
d[, compartment := factor(compartment, levels = c("Serum","EP"))]
d[, p_label := fifelse(p < 0.001, "p < 0.001", paste0("p = ", signif(p, 2)))]
cat("\npanel C:\n"); print(d)

pC <- ggplot(d, aes(x = rho, y = protein, shape = compartment)) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_point(size = 4, position = position_dodge(width = 0.55)) +
  geom_text(aes(label = p_label), position = position_dodge(width = 0.55),
            hjust = -0.15, size = 3.5, family = fnt) +
  labs(x = expression(paste("Spearman's rho with ", Delta, "VO"[2], "max")),
       y = NULL, shape = NULL) +
  xlim(min(d$rho) - 0.1, max(d$rho) + 0.35) +
  theme_ecv()
save_ecv("R1_25_fig5B_compartment.pdf", pC, width = 16, height = 8, unit = "cm")

# ---- panel D: Venn med ærlig nevner ----
EP <- unique(readRDS(file.path(OUT,"R1_01_longterm.rds"))$PRIM[
  !is.na(Genes) & Genes != "", .(Genes, q_ep = q)], by = "Genes")
SE <- fread(file.path(OLINK_DIR,"Exercise/MyoGlu_Olink_A1B1_All_Annotated.csv"), sep = ";")
SE <- unique(SE[!is.na(GeneName) & GeneName != "", .(Genes = GeneName, q_se = q)], by = "Genes")
b <- merge(EP, SE, by = "Genes")
e <- b[q_ep < .05, Genes]; s <- b[q_se < .05, Genes]
cat(sprintf("\npanel D: målbare på begge %d | EP-sig %d | serum-sig %d | overlapp %d\n",
            nrow(b), length(e), length(s), length(intersect(e, s))))
# Venn tegnes manuelt: ggvenn viser prosent av UNION, som er samme
# nevner-feil R3 påpeker for akutt/kronisk-overlappet. Her vises i stedet
# overlappet som andel av EP-settet — altså hvor mye av EP-responsen som
# lar seg gjenfinne i serum.
nE <- length(e); nS <- length(s); nO <- length(intersect(e, s))
circ <- function(x0, r = 1.05, n = 200) { t <- seq(0, 2*pi, length.out = n)
  data.table(x = x0 + r*cos(t), y = r*sin(t)) }
cd <- rbind(circ(-0.62)[, set := "EP"], circ(0.62)[, set := "Serum"])
lb <- data.table(
  x = c(-1.15, 0, 1.15, -1.15, 1.15),
  y = c(0, 0, 0, 1.35, 1.35),
  lab = c(sprintf("%d\n(%.1f%%)", nE - nO, 100*(nE-nO)/nE),
          sprintf("%d\n(%.1f%%)", nO, 100*nO/nE),
          sprintf("%d\n(%.1f%%)", nS - nO, 100*(nS-nO)/nS),
          "EP", "Serum"),
  sz = c(4, 4, 4, 4.6, 4.6))
pD <- ggplot() +
  geom_polygon(data = cd, aes(x, y, fill = set, group = set),
               alpha = 0.55, colour = "black", linewidth = 0.6) +
  scale_fill_manual(values = setNames(PAL$venn, c("EP","Serum"))) +
  geom_text(data = lb, aes(x, y, label = lab, size = sz), family = fnt) +
  scale_size_identity() + coord_equal() +
  theme_void() + theme(legend.position = "none",
                       text = element_text(family = fnt))
cat(sprintf("Venn-etiketter: EP-only %d (%.1f%% av EP) | overlapp %d (%.1f%% av EP) | serum-only %d (%.1f%% av serum)\n",
            nE-nO, 100*(nE-nO)/nE, nO, 100*nO/nE, nS-nO, 100*(nS-nO)/nS))
exp_ov <- nE * nS / nrow(b)
cat(sprintf("anrikning: %.1f observert / %.1f forventet = %.1fx (Fisher p = %.3g)\n",
            nO, exp_ov, nO/exp_ov,
            fisher.test(matrix(c(nO, nE-nO, nS-nO, nrow(b)-nE-nS+nO), 2))$p.value))
save_ecv("R1_25_fig5C_venn.pdf", pD, width = 3.5, height = 3.5, unit = "in")
fwrite(d, file.path(OUT, "R1_25_fig5B_stats.csv"))
