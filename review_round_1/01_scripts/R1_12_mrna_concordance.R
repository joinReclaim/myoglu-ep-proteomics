# =====================================================================
# R1_12_mrna_concordance.R — systematisk ΔEP-protein mot ΔmRNA
# Adresserer R5 punkt 8 («correlate ΔEP protein with ΔmRNA for every
# protein with a measured transcript and compare against a permuted
# null»), R2 concern 3 (LTBP4-diskordansen)
#
# Manuskriptet undersøker fire kandidater. Her gjøres den systematiske
# vurderingen diskusjonen selv sier er nødvendig, mot permutert null.
# Krever: R1_00_setup.R, R1_01_longterm.R
# =====================================================================
suppressPackageStartupMessages({library(data.table); library(stringr)})
source(file.path(Sys.getenv("ECV_ROOT", "/Users/sindrle/Research/Prosjekter/ECV"),
                 "review_round_1/01_scripts/R1_config.R"))
ROOT <- ECV_ROOT
OUT  <- file.path(ROOT, "review_round_1/02_output"); setwd(ROOT); set.seed(1)
source(file.path(ROOT, "review_round_1/01_scripts/R1_theme.R"))
S <- readRDS(file.path(OUT, "R1_00_normalised.rds"))
P <- readRDS(file.path(OUT, "R1_01_longterm.rds"))$PRIM

# ---- ΔEP per deltaker (loess, >=70% completeness) ----
i <- which(S$si$SampleCode %in% c("A1","B1")); M <- S$Xloe[, i]; s <- S$si[i]
nA <- rowSums(!is.na(M[, s$SampleCode=="A1"])); nB <- rowSums(!is.na(M[, s$SampleCode=="B1"]))
M <- M[nA>=.7*26 & nB>=.7*26, ]; g <- S$genes[nA>=.7*26 & nB>=.7*26]
ids <- sort(unique(s$ID2))
dEP <- sapply(ids, function(x){
  a <- M[, s$SampleCode=="A1" & s$ID2==x, drop=FALSE]; b <- M[, s$SampleCode=="B1" & s$ID2==x, drop=FALSE]
  if(!ncol(a)||!ncol(b)) rep(NA_real_, nrow(M)) else b[,1]-a[,1] })
rownames(dEP) <- make.unique(g)

# ---- ΔmRNA per deltaker (samme ID-konvensjon som 202_MixOmics_DIABLO.R) ----
load(RNASEQ_RDA)
dmrna <- function(mat){
  d <- as.data.frame(mat); sid <- rownames(d)
  tm <- str_extract(sid, "A1|B1"); raw <- as.numeric(str_extract(sid, "\\d+"))
  id <- ifelse(tm=="A1" & raw>=100, raw-100, ifelse(tm=="B1" & raw>=200, raw-200, raw))
  keep <- !is.na(tm) & !is.na(id); d <- d[keep,]; tm <- tm[keep]; id <- id[keep]
  out <- sapply(ids, function(x){
    a <- which(id==x & tm=="A1"); b <- which(id==x & tm=="B1")
    if(!length(a)||!length(b)) rep(NA_real_, ncol(d))
    else log2(as.numeric(d[b[1],])+1) - log2(as.numeric(d[a[1],])+1) })
  rownames(out) <- colnames(d); out }
dSkM <- dmrna(fpkm.sm.fix); dAT <- dmrna(fpkm.at.fix)
cat(sprintf("ΔEP %d proteiner | ΔmRNA muskel %d | ΔmRNA fett %d | deltakere %d\n\n",
            nrow(dEP), nrow(dSkM), nrow(dAT), length(ids)))

concord <- function(dm, lab){
  sh <- intersect(rownames(dEP), rownames(dm))
  r <- sapply(sh, function(x){
    a <- dEP[x,]; b <- dm[x,]; ok <- is.finite(a)&is.finite(b)
    if(sum(ok)<10) NA_real_ else suppressWarnings(cor(a[ok],b[ok],method="spearman")) })
  r <- r[!is.na(r)]
  null <- replicate(200, { pi_ <- sample(seq_along(ids))
    q <- sapply(names(r), function(x){ a <- dEP[x,]; b <- dm[x,pi_]
      ok <- is.finite(a)&is.finite(b)
      if(sum(ok)<10) NA_real_ else suppressWarnings(cor(a[ok],b[ok],method="spearman")) })
    mean(q, na.rm=TRUE) })
  pv <- (1+sum(abs(null) >= abs(mean(r))))/(length(null)+1)
  cat(sprintf("=== %s ===\n", lab))
  cat(sprintf("  %d gener med både protein og transkript\n", length(r)))
  cat(sprintf("  gjennomsnittlig rho: %+.4f | permutert null: %+.4f (SD %.4f) | p = %.3f\n",
              mean(r), mean(null), sd(null), pv))
  cat(sprintf("  andel positiv: %.0f%% | |rho|>0.5: %.1f%% \n",
              100*mean(r>0), 100*mean(abs(r)>0.5)))
  list(r=r, null=null, lab=lab) }
cm <- concord(dSkM, "SKJELETTMUSKEL"); ca <- concord(dAT, "FETTVEV")

cat("\n=== KANDIDATENE ===\n")
kk <- rbindlist(lapply(c("LTBP4","VCL","MEP1B","HSP90B1","TPI1"), function(x)
  data.table(Genes=x, rho_muskel=round(unname(cm$r[x]),3), rho_fett=round(unname(ca$r[x]),3))))
print(kk)

cat("\n=== ER RESPONDERNE ANRIKET FOR VEVSUTTRYKTE GENER? ===\n")
sig <- unique(P[q<.05, Genes])
for (nm in c("muskel","fett")) {
  mat <- if(nm=="muskel") fpkm.sm.fix else fpkm.at.fix
  expr <- colMeans(mat, na.rm=TRUE)
  sh <- intersect(rownames(dEP), names(expr))
  a <- expr[intersect(sh, sig)]; b <- expr[setdiff(sh, sig)]
  w <- wilcox.test(a, b)
  cat(sprintf("  %-7s median FPKM respondere %.1f vs øvrige målte %.1f (n=%d/%d, p=%.3g)\n",
              nm, median(a,na.rm=TRUE), median(b,na.rm=TRUE), length(a), length(b), w$p.value)) }

fwrite(data.table(Genes=names(cm$r), rho_SkM=cm$r,
                  rho_AT=ca$r[match(names(cm$r), names(ca$r))]),
       file.path(OUT,"R1_12_mrna_protein_concordance.csv"))

# ---- figur ----
d <- rbind(data.table(rho=cm$r, vev="Skjelettmuskel"), data.table(rho=ca$r, vev="Fettvev"))
p <- ggplot(d, aes(rho, fill=vev)) +
  geom_histogram(bins=60, alpha=.75, position="identity", colour=NA) +
  geom_vline(xintercept=0, linetype="dotted") +
  scale_fill_manual(values=cols) +
  labs(x=expression(paste("Spearman rho: ",Delta,"EP-protein mot ",Delta,"mRNA")), y="Antall gener") +
  theme_ecv() + theme(legend.position="top")
save_ecv("R1_12_mRNA_protein_concordance.pdf", p, width=12, height=8)
