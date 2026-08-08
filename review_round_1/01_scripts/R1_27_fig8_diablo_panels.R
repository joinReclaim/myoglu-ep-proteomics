# =====================================================================
# R1_27_fig8_diablo_panels.R — manuskriptets FIGUR 8 (Figures/Figure6.tif)
# (blir figur 7 etter omnummerering)
# Erstatter 202_MixOmics_DIABLO.R:560-700
#
# Uendret stil i panel A: comp1 mot comp2 per blokk, grå linjer som knytter
# samme deltaker før/etter, farge = tidspunkt (#1F77B4 / #FF7F0E), form =
# gruppe, facet_wrap ncol=2, theme_classic.
# Endret:
#   - «Serum ECV proteins» -> «Plasma EP proteins» (Q- to be fixed)
#   - loess-normalisert EP-matrise
#   - blokkene bygges USENTRERT (ingen withinVariation) — se kommentar i koden
#   - 4 fettvevsprøver som manglet HELT er droppet i stedet for å bli
#     erstattet med gruppe x tid-gjennomsnitt (202:259-282 fabrikkerte dem)
#   - n = 48 prøver fra 24 deltakere
#   - circos-cutoff oppgis eksplisitt (R5 minor: var den a priori?)
# Krever: R1_13_diablo.R
# =====================================================================
suppressPackageStartupMessages({library(data.table); library(mixOmics); library(stringr)})
source(file.path(Sys.getenv("ECV_ROOT", "/Users/sindrle/Research/Prosjekter/ECV"),
                 "review_round_1/01_scripts/R1_config.R"))
ROOT <- ECV_ROOT
OUT  <- file.path(ROOT,"review_round_1/02_output"); setwd(ROOT); set.seed(1)
source(file.path(ROOT,"review_round_1/01_scripts/R1_theme.R"))
NCOMP <- 2; KEEPX <- 10; CUTOFF <- 0.6

# ---- blokker bygges USENTRERT ----
# R5 punkt 8: separasjonen er «expected after supervised multi-block feature
# selection with within-participant centring». Uten sentrering faller den
# halvdelen av innvendingen bort, og separasjonen som vises er reell (om enn
# ufullstendig) framfor garantert av dekomposisjonen.
S <- readRDS(file.path(OUT,"R1_00_normalised.rds"))
i <- which(S$si$SampleCode %in% c("A1","B1")); M <- S$Xloe[, i]; s <- S$si[i]
rownames(M) <- make.unique(S$genes); M <- M[rowSums(!is.na(M)) == ncol(M), ]
load(RNASEQ_RDA)
mk <- function(mat){ d <- as.data.frame(mat); sid <- rownames(d)
  tm <- str_extract(sid,"A1|B1"); raw <- as.numeric(str_extract(sid,"\\d+"))
  id <- ifelse(tm=="A1"&raw>=100, raw-100, ifelse(tm=="B1"&raw>=200, raw-200, raw))
  ok <- !is.na(tm)&!is.na(id); d <- d[ok,]; rownames(d) <- paste0(id[ok],"_",tm[ok])
  log2(as.matrix(d)+1) }
SM <- mk(fpkm.sm.fix); AT <- mk(fpkm.at.fix)
ol <- fread(file.path(OLINK_DIR,"MyoGlu_Olink_A1B1_NPX.csv"), sep=";")
olm <- as.matrix(ol[, !c("ID","Time","group"), with=FALSE]); storage.mode(olm) <- "double"
rownames(olm) <- paste0(ol$ID,"_",ol$Time); olm <- olm[, colSums(!is.na(olm))==nrow(olm), drop=FALSE]
key <- paste0(s$ID2,"_",s$SampleCode)
atmiss <- rownames(AT)[rowSums(!is.na(AT))==0]
cat(sprintf("fettvevsprøver som mangler helt og droppes: %d (%s)\n",
            length(atmiss), paste(atmiss, collapse=", ")))
keep <- key %in% rownames(SM) & key %in% rownames(AT) & key %in% rownames(olm) & !(key %in% atmiss)
key <- key[keep]; M <- M[,keep]; s <- s[keep]
topv <- function(Z,n=500){ v <- apply(Z,2,var,na.rm=TRUE); Z[, order(v,decreasing=TRUE)[1:min(n,ncol(Z))], drop=FALSE] }
blkf <- function(Z){ Z <- Z[key,,drop=FALSE]; rownames(Z) <- key; topv(Z) }
X <- list(EP = topv(`rownames<-`(t(M), key)), Serum = blkf(olm),
          Muskel = blkf(SM), Fett = blkf(AT))
X <- lapply(X, function(z) z[, apply(z,2,function(x) all(is.finite(x)) && var(x)>0), drop=FALSE])
# Nivåene navngis som i manuskriptet — circosPlot og cimDiablo viser dem
# rått i «Expression»-legenden, og «A1»/«B1» er interne koder.
# «Baseline»/«Post-training» og ikke «Before»/«After»: cimDiablo sorterer
# legendetiketter alfabetisk men fargeruter i faktorrekkefølge, så etikettene
# må være alfabetisk ordnet kronologisk for at legenden skal stemme.
Y <- factor(s$SampleCode, levels = c("A1","B1"),
            labels = c("Baseline","Post-training"))

names(X) <- c("Plasma EP proteins","Serum proteins","Skeletal muscle mRNA","Adipose tissue mRNA")

# Serum-blokken har OlinkID som kolonnenavn (OID31501 osv). Originalfiguren
# viste disse rått i circos-plottet. Mappes til gennavn via Annotation.csv.
an <- fread(file.path(OLINK_DIR,"Annotation.csv"), sep = ";", quote = "")
gmap <- setNames(an$GeneName, an$OlinkID)
cn <- colnames(X[["Serum proteins"]])
hit <- cn %in% names(gmap)
colnames(X[["Serum proteins"]])[hit] <- make.unique(unname(gmap[cn[hit]]))
cat(sprintf("serum-blokk: %d av %d OlinkID mappet til gennavn\n", sum(hit), length(cn)))
des <- matrix(0.1, length(X), length(X)); diag(des) <- 0
kx <- lapply(X, function(z) rep(min(KEEPX, ncol(z)), NCOMP))
fit <- block.splsda(X, Y, ncomp = NCOMP, keepX = kx, design = des)
cat(sprintf("blokker: %s | prøver: %d (%d deltakere)\n",
            paste(names(X), collapse = ", "), nrow(X[[1]]), uniqueN(s$ID2)))

# ---- panel A ----
d <- rbindlist(lapply(names(X), function(b) data.table(
  comp1 = fit$variates[[b]][,1], comp2 = fit$variates[[b]][,2],
  block = b, ID2 = s$ID2,
  Time = fifelse(s$SampleCode == "A1", "Baseline", "Post-training"),
  Gr   = fifelse(s$Gr == "dysglyc", "Overweight", "Control"))))
d[, block := factor(block, levels = names(X))]
d[, Time  := factor(Time, levels = c("Baseline","Post-training"))]
d[, Gr    := factor(Gr, levels = c("Overweight","Control"))]

# kildedata for panel A, til kodedeponeringen. To komponentskårer per deltaker
# per blokk, subjekt-sentrert og komprimert fra 500 features — ingenting lar seg
# rekonstruere. Deltakerkoden er den samme som i ProteomeXchange-mappingen.
# Panel B deponeres IKKE: der er radene enkeltfeatures fra Olink og RNA-seq.
fwrite(d[, .(subject_id = sprintf("MyoGlu_%03d", ID2), timepoint = Time,
             group = Gr, block, comp1, comp2)],
       file.path(OUT, "R1_27_fig8A_scores.csv"))

lim <- ceiling(max(abs(c(d$comp1, d$comp2)), na.rm = TRUE) * 1.1)

pA <- ggplot(d, aes(comp1, comp2)) +
  geom_hline(yintercept = 0, linewidth = 0.25, colour = "grey80") +
  geom_vline(xintercept = 0, linewidth = 0.25, colour = "grey80") +
  geom_line(aes(group = interaction(block, ID2)), colour = "grey70",
            linewidth = 0.35, alpha = 0.45) +
  geom_point(aes(colour = Time, shape = Gr), size = 2.8, alpha = 0.9) +
  facet_wrap(~ block, ncol = 2) +
  scale_colour_manual(values = PAL$time) +
  coord_cartesian(xlim = c(-lim, lim), ylim = c(-lim, lim)) +
  labs(x = "DIABLO component 1", y = "DIABLO component 2",
       colour = NULL, shape = NULL) +
  theme_classic(base_size = 12) +
  theme(text = element_text(family = fnt), strip.background = element_blank(),
        legend.position = "bottom", legend.box = "vertical",
        axis.text = element_text(colour = "black"))
save_ecv("R1_27_fig8A_individuals.pdf", pA, width = 15, height = 16, unit = "cm")

# ---- panel B: CIM over komponent 1 ----
f <- file.path(FIGDIR, "R1_27_fig8B_cim.pdf")
cairo_pdf(f, width = 9, height = 7)
# color.Y må følge faktornivåene til Y (Before, After) så fargene matcher
# panel A. row.names = FALSE skjuler de rå prøvekodene (1_B1, 21_A1 ...),
# som i originalfiguren.
try(cimDiablo(fit, comp = 1, margins = c(10, 16), size.legend = 0.7,
              color.Y = unname(PAL$time[levels(Y)]),
              row.names = FALSE, legend.position = "topright"), silent = TRUE)
dev.off(); message("figur: ", f)

# ---- panel C: circos ----
f2 <- file.path(FIGDIR, "R1_27_fig8C_circos.pdf")
# større lerret + mindre etiketter: circos ble beskåret i sidene ved 8x8
cairo_pdf(f2, width = 11, height = 11)
par(mar = c(2, 2, 2, 2))
try(circosPlot(fit, cutoff = CUTOFF, line = TRUE, size.labels = 0.9,
               size.variables = 0.55, size.legend = 0.8,
               color.blocks = c("#D95F02","#1B9E77","#7570B3","#E7298A")),
    silent = TRUE)
dev.off(); message("figur: ", f2)
# ---- cutoff-følsomhet (R5 minor: var r >= 0.6 valgt a priori?) ----
# NB: circosPlot bruker IKKE rå feature-feature-korrelasjoner, men en
# similaritetsmatrise utledet av variablenes projeksjon på komponentene.
# Følsomheten må derfor beregnes på DEN matrisen, ikke på rådata.
cat(sprintf("\ncircos-cutoff = %.1f — IKKE valgt a priori (R5 minor).\n", CUTOFF))
pdf(NULL); cp <- circosPlot(fit, cutoff = 0, line = FALSE); invisible(dev.off())
sv <- lapply(names(X), function(b)
  unique(unlist(lapply(seq_len(NCOMP), function(k) selectVar(fit, comp = k)[[b]]$name))))
names(sv) <- names(X)
# NB: circosPlot limer blokknavnet på radetikettene ("GEN EP proteins"), så
# oppslag på featurenavn feiler. Radene ligger i blokkrekkefølge — tilordnes
# posisjonelt, med kontroll på at featurenavnet faktisk hører til blokken.
b2 <- rep(names(sv), sapply(sv, length))
stopifnot(length(b2) == nrow(cp))
ok_map <- mapply(function(rn, b) any(startsWith(rn, sv[[b]])), rownames(cp), b2)
cat(sprintf("blokktilordning verifisert for %d av %d features\n", sum(ok_map), length(ok_map)))
cross <- outer(b2, b2, "!=") & upper.tri(cp)
cat(sprintf("similaritetsmatrise %d x %d | kryssblokk-par: %d\n",
            nrow(cp), ncol(cp), sum(cross)))
for (cc in c(0.4, 0.5, 0.6, 0.7, 0.8))
  cat(sprintf("   r >= %.1f : %4d (%.0f%%)\n", cc,
              sum(abs(cp)[cross] >= cc), 100*mean(abs(cp)[cross] >= cc)))
fwrite(data.table(cutoff = c(.4,.5,.6,.7,.8),
                  n = sapply(c(.4,.5,.6,.7,.8), function(cc) sum(abs(cp)[cross] >= cc)),
                  pct = sapply(c(.4,.5,.6,.7,.8), function(cc) round(100*mean(abs(cp)[cross] >= cc)))),
       file.path(OUT, "R1_27_fig8C_cutoff_sensitivity.csv"))

sel <- rbindlist(lapply(names(X), function(b)
  data.table(block = b, comp = 1, feature = selectVar(fit, comp = 1)[[b]]$name)))
fwrite(sel, file.path(OUT, "R1_27_fig8_selected_features.csv"))
cat(sprintf("utvalgte features komponent 1: %s\n",
            paste(sprintf("%s %d", names(X), sel[, .N, by = block]$N), collapse = " | ")))
