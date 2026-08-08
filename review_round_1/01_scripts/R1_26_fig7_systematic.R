# =====================================================================
# R1_26_fig7_systematic.R — NY FIGUR 7 (blir figur 6 etter omnummerering)
# Erstatter de tolv bokspanelene i Figures/Figure5_new.tif
#
# Begrunnelse: EP-responsene for kandidatene vises allerede i figur 3
# (volcano) og figur 5 (lollipop). Det unike for denne figuren er vevs-
# relasjonen — og den globale analysen viser at det ikke finnes noen
# systematisk sammenheng. Fire eksempler valgt fra en fordeling sentrert
# på null er nettopp det R5 punkt 8 kritiserer. Figuren handler derfor nå
# om det systematiske spørsmålet; enkeltpanelene flyttes til supplementary.
#
# Panel A: ΔEP-protein mot ΔmRNA for ALLE gener med målt transkript, mot
#          permutert null (R5 punkt 8: «compare against a permuted null»)
# Panel B: vevsopprinnelse testet framfor påstått — er responderne anriket
#          for vevsuttrykte gener? (R5 punkt 8)
# Krever: R1_00_setup.R, R1_01_longterm.R
# =====================================================================
suppressPackageStartupMessages({library(data.table); library(stringr)})
source(file.path(Sys.getenv("ECV_ROOT", "/Users/sindrle/Research/Prosjekter/ECV"),
                 "review_round_1/01_scripts/R1_config.R"))
ROOT <- ECV_ROOT
OUT  <- file.path(ROOT,"review_round_1/02_output"); setwd(ROOT); set.seed(1)
source(file.path(ROOT,"review_round_1/01_scripts/R1_theme.R"))
S <- readRDS(file.path(OUT,"R1_00_normalised.rds"))
P <- readRDS(file.path(OUT,"R1_01_longterm.rds"))$PRIM

# ---- ΔEP per deltaker ----
i <- which(S$si$SampleCode %in% c("A1","B1")); M <- S$Xloe[, i]; s <- S$si[i]
nA <- rowSums(!is.na(M[, s$SampleCode=="A1"])); nB <- rowSums(!is.na(M[, s$SampleCode=="B1"]))
M <- M[nA>=.7*26 & nB>=.7*26, ]; g <- S$genes[nA>=.7*26 & nB>=.7*26]
ids <- sort(unique(s$ID2))
dEP <- sapply(ids, function(x){
  a <- M[, s$SampleCode=="A1" & s$ID2==x, drop=FALSE]; b <- M[, s$SampleCode=="B1" & s$ID2==x, drop=FALSE]
  if(!ncol(a)||!ncol(b)) rep(NA_real_, nrow(M)) else b[,1]-a[,1] })
rownames(dEP) <- make.unique(g)

load(RNASEQ_RDA)
dmrna <- function(mat){
  d <- as.data.frame(mat); sid <- rownames(d)
  tm <- str_extract(sid,"A1|B1"); raw <- as.numeric(str_extract(sid,"\\d+"))
  id <- ifelse(tm=="A1"&raw>=100, raw-100, ifelse(tm=="B1"&raw>=200, raw-200, raw))
  k <- !is.na(tm)&!is.na(id); d <- d[k,]; tm <- tm[k]; id <- id[k]
  out <- sapply(ids, function(x){ a <- which(id==x&tm=="A1"); b <- which(id==x&tm=="B1")
    if(!length(a)||!length(b)) rep(NA_real_, ncol(d))
    else log2(as.numeric(d[b[1],])+1) - log2(as.numeric(d[a[1],])+1) })
  rownames(out) <- colnames(d); out }
dSkM <- dmrna(fpkm.sm.fix); dAT <- dmrna(fpkm.at.fix)

rho_set <- function(dm, perm = FALSE){
  sh <- intersect(rownames(dEP), rownames(dm))
  pi_ <- if (perm) sample(seq_along(ids)) else seq_along(ids)
  r <- sapply(sh, function(x){ a <- dEP[x,]; b <- dm[x, pi_]
    ok <- is.finite(a)&is.finite(b)
    if(sum(ok)<10) NA_real_ else suppressWarnings(cor(a[ok],b[ok],method="spearman")) })
  r[!is.na(r)] }
obs <- list(SkM = rho_set(dSkM), AT = rho_set(dAT))
nul <- list(SkM = unlist(lapply(1:5, function(z) rho_set(dSkM, TRUE))),
            AT  = unlist(lapply(1:5, function(z) rho_set(dAT,  TRUE))))
for (t in names(obs)) cat(sprintf("%-4s observert: n=%d, gj.snitt rho %+.4f | null: %+.4f\n",
                                  t, length(obs[[t]]), mean(obs[[t]]), mean(nul[[t]])))

lab <- c(SkM = "Skeletal muscle", AT = "Adipose tissue")
dA <- rbind(
  rbindlist(lapply(names(obs), function(t) data.table(rho = obs[[t]], tissue = lab[t], set = "Observed"))),
  rbindlist(lapply(names(nul), function(t) data.table(rho = nul[[t]], tissue = lab[t], set = "Permuted null"))))
dA[, tissue := factor(tissue, levels = lab)]
pA <- ggplot(dA, aes(rho, colour = set, fill = set)) +
  geom_density(alpha = 0.25, linewidth = 0.7) +
  geom_vline(xintercept = 0, linetype = "dotted") +
  facet_wrap(~ tissue) +
  scale_x_continuous(breaks = c(-0.5, 0, 0.5), limits = c(-0.85, 0.85)) +
  scale_colour_manual(values = c(Observed = cols[1], `Permuted null` = "grey45")) +
  scale_fill_manual(values   = c(Observed = cols[1], `Permuted null` = "grey75")) +
  labs(x = expression(paste("Spearman's rho: ", Delta, "EP protein vs ", Delta, "mRNA")),
       y = "Density") +
  theme_ecv() + theme(legend.position = "top")
save_ecv("R1_26_fig7A_concordance.pdf", pA, width = 17, height = 8, unit = "cm")

# ---- panel B: vevsopprinnelse ----
sig <- unique(P[q < .05, Genes])
dB <- rbindlist(lapply(names(lab), function(t){
  mat <- if (t == "SkM") fpkm.sm.fix else fpkm.at.fix
  e <- colMeans(mat, na.rm = TRUE)
  sh <- intersect(rownames(dEP), names(e))
  data.table(expr = log2(e[sh] + 1), tissue = lab[t],
             grp = fifelse(sh %in% sig, "Responders", "Other measured")) }))
dB[, grp := factor(grp, levels = c("Responders","Other measured"))]
dB[, tissue := factor(tissue, levels = lab)]
st <- dB[, { w <- wilcox.test(expr ~ grp)
  .(p = w$p.value, n_resp = sum(grp=="Responders"), n_oth = sum(grp=="Other measured"),
    med_resp = median(expr[grp=="Responders"]), med_oth = median(expr[grp=="Other measured"])) },
  by = tissue]
cat("\npanel B:\n"); print(st)
pB <- ggplot(dB, aes(grp, expr, fill = grp)) +
  geom_boxplot(outlier.shape = NA, width = 0.55, colour = "black", alpha = 0.75) +
  facet_wrap(~ tissue) +
  scale_fill_manual(values = c(Responders = cols[1], `Other measured` = "grey80")) +
  labs(x = NULL, y = expression(paste("Tissue expression, log"[2], "(FPKM + 1)"))) +
  theme_ecv() + theme(legend.position = "none",
                      axis.text.x = element_text(angle = 20, hjust = 1))
save_ecv("R1_26_fig7B_tissue_origin.pdf", pB, width = 11, height = 9, unit = "cm")
fwrite(st, file.path(OUT, "R1_26_fig7B_stats.csv"))
