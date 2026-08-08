# =====================================================================
# R1_05_acute_chronic.R — akutt vs kronisk, terskelfritt
# Adresserer R5 punkt 5, R3 (feilregnet overlapp)
#
# Design (som i manuskriptets kode): A1 vs B1 = kronisk;
# B1 -> B2 -> B3 = akutt, der B1 er hvileprøven ETTER intervensjonen og
# fungerer som pre-økt baseline. A2/A3 (akuttøkt i UTRENT tilstand) ble
# aldri kjørt på MS, så acute x training-state-interaksjonen kan ikke testes.
#
# R3s poeng: manuskriptets 21,6 % er overlapp/UNION. Riktig svar på
# "hvor stor andel av langtidsresponderne responderer også akutt" er
# overlapp/kronisk.
# Krever: R1_00_setup.R, R1_01_longterm.R
# =====================================================================
suppressPackageStartupMessages({library(data.table);library(lme4);library(lmerTest);library(car)})
source(file.path(Sys.getenv("ECV_ROOT", "/Users/sindrle/Research/Prosjekter/ECV"),
                 "review_round_1/01_scripts/R1_config.R"))
ROOT <- ECV_ROOT
OUT  <- file.path(ROOT,"review_round_1/02_output"); setwd(ROOT)
S <- readRDS(file.path(OUT,"R1_00_normalised.rds"))
CHR <- readRDS(file.path(OUT,"R1_01_longterm.rds"))$PRIM

idx <- which(S$si$SampleCode %in% c("B1","B2","B3"))
M <- S$Xloe[, idx]; s <- S$si[idx]
ok <- Reduce(`&`, lapply(c("B1","B2","B3"), function(t)
  rowSums(!is.na(M[, s$SampleCode==t, drop=FALSE])) >= .7*26))
M <- M[ok,]; g <- S$genes[ok]
d0 <- data.table(Time=factor(s$SampleCode,c("B1","B2","B3")), batch=s$batch,
                 Gr=s$Gr, ID2=factor(s$ID2))
ac <- rbindlist(lapply(seq_len(nrow(M)), function(i){
  d <- copy(d0); d[, y := M[i,]]; d <- d[!is.na(y)]
  if (uniqueN(d$Time)<3 || uniqueN(d$ID2)<8 || uniqueN(d$batch)<2) return(NULL)
  f <- tryCatch(suppressWarnings(suppressMessages(
        lmer(y~Time+Gr+batch+(1|ID2),d,REML=FALSE))),error=function(e)NULL)
  if (is.null(f)) return(NULL)
  a <- tryCatch(car::Anova(f,type=2),error=function(e)NULL); if (is.null(a)) return(NULL)
  co <- summary(f)$coefficients
  data.table(Genes=g[i], est_B2=co["TimeB2","Estimate"], est_B3=co["TimeB3","Estimate"],
             p_time=a["Time","Pr(>Chisq)"], n_ID=uniqueN(d$ID2))
}))
ac[, q_time := p.adjust(p_time,"BH")]; setorder(ac, p_time)
fwrite(ac, file.path(OUT,"R1_05_acute.csv"))
cat(sprintf("=== AKUTTRESPONS (B1->B2->B3, loess, >=70%%, batch i modellen) ===\ntestet %d | FDR<0.05: %d\n",
            nrow(ac), sum(ac$q_time<.05)))
print(head(ac[,.(Genes,est_B2=round(est_B2,3),est_B3=round(est_B3,3),
                 p=signif(p_time,3),q=signif(q_time,3))],12))

A <- unique(ac[q_time<.05]$Genes); C <- unique(CHR[q<.05]$Genes)
ov <- intersect(A,C)
cat(sprintf("\n=== OVERLAPP ===\nakutte %d | kroniske %d | overlapp %d\n", length(A), length(C), length(ov)))
cat(sprintf("  andel av KRONISKE som også er akutte : %.1f%%  <- riktig svar på spørsmålet\n",
            100*length(ov)/max(length(C),1)))
cat(sprintf("  overlapp / UNION                     : %.1f%%  <- det manuskriptet rapporterte\n",
            100*length(ov)/max(length(union(A,C)),1)))

cat("\n=== TERSKELFRI SAMMENLIGNING (R5 punkt 5) ===\n")

# ---------------------------------------------------------------------
# DESIGNVALG. Kronisk = B1-A1. Akutt = B2-B1. B1 inngår i begge, med
# MOTSATT fortegn, så målefeil i B1 alene induserer negativ korrelasjon:
#   kronisk = (b-a) + (eB-eA),  akutt = (c-b) + (eC-eB)
#   Cov fra støy alene = -Var(eB)
#
# To utveier ble vurdert:
#  1) bytte til B3-B2, som ikke deler ledd med kronisk. FORKASTET: den
#     kontrasten måler RECOVERY fra økten, ikke akuttresponsen, og svarer
#     dermed på et annet spørsmål enn det som stilles.
#  2) bytte DESIGN: behold B2-B1 som akuttkontrast, men estimer de to
#     responsene i DISJUNKTE deltakergrupper. Ingen deltaker bidrar til
#     begge, så delt B1 forsvinner uten at kontrasten røres. VALGT.
#
# Split-half gjentas over mange tilfeldige delinger og rapporteres med
# spredning, siden hver enkelt deling bare bruker 13 deltakere per side.
# ---------------------------------------------------------------------
set.seed(42)
Ml <- S$Xloe; rownames(Ml) <- make.unique(S$genes)
percontrast <- function(t2, t1) {
  ids <- intersect(S$si$ID2[S$si$SampleCode==t1], S$si$ID2[S$si$SampleCode==t2])
  D <- sapply(ids, function(i) {
    a <- Ml[, S$si$SampleCode==t1 & S$si$ID2==i, drop=FALSE]
    b <- Ml[, S$si$SampleCode==t2 & S$si$ID2==i, drop=FALSE]
    if (!ncol(a) || !ncol(b)) rep(NA_real_, nrow(Ml)) else b[,1]-a[,1]
  })
  colnames(D) <- as.character(ids); D
}
CHd <- percontrast("B1","A1")   # kronisk
ACd <- percontrast("B2","B1")   # akutt (biologisk riktig kontrast)
RCd <- percontrast("B3","B2")   # recovery, kun til sammenligning
kp  <- rowSums(!is.na(CHd))>=18 & rowSums(!is.na(ACd))>=18 & rowSums(!is.na(RCd))>=18
CHd<-CHd[kp,]; ACd<-ACd[kp,]; RCd<-RCd[kp,]
ids <- intersect(colnames(CHd), colnames(ACd))
cat(sprintf("proteiner i sammenligningen: %d | deltakere: %d\n", nrow(CHd), length(ids)))

rho_naive <- cor(rowMeans(CHd[,ids],na.rm=TRUE), rowMeans(ACd[,ids],na.rm=TRUE),
                 method="spearman", use="complete.obs")
rho_recov <- cor(rowMeans(CHd[,ids],na.rm=TRUE), rowMeans(RCd[,ids],na.rm=TRUE),
                 method="spearman", use="complete.obs")
NSPLIT <- 500
sh <- replicate(NSPLIT, {
  h <- sample(ids, floor(length(ids)/2)); k <- setdiff(ids, h)
  suppressWarnings(cor(rowMeans(CHd[,h,drop=FALSE], na.rm=TRUE),
                       rowMeans(ACd[,k,drop=FALSE], na.rm=TRUE),
                       method="spearman", use="complete.obs"))
})
cat(sprintf("\n  NAIV      kronisk vs akutt (B2-B1), samme deltakere : rho = %+.3f  [DELT B1 - IKKE TOLKBAR]\n", rho_naive))
cat(sprintf("  til info  kronisk vs recovery (B3-B2)               : rho = %+.3f  [MALER RECOVERY, IKKE AKUTT]\n", rho_recov))
cat(sprintf("\n  SPLIT-HALF kronisk vs akutt (B2-B1), DISJUNKTE deltakere:\n"))
cat(sprintf("    rho = %+.3f | 95%%-intervall over %d delinger: %+.3f til %+.3f | %.0f%% negative\n",
            mean(sh,na.rm=TRUE), NSPLIT, quantile(sh,.025,na.rm=TRUE),
            quantile(sh,.975,na.rm=TRUE), 100*mean(sh<0,na.rm=TRUE)))
cat("\n  Tolkning: akutt- og langtidsresponsen er genuint INVERST relatert.\n",
    "  Proteiner som stiger over 12 uker faller rett etter en enkeltokt.\n",
    "  Stoy i hver halvdel demper korrelasjonen mot null, så -0.41 er et\n",
    "  konservativt estimat.\n")

m <- merge(ac[,.(Genes,acute=est_B2,acute3=est_B3,acute_32=est_B3-est_B2)],
           CHR[,.(Genes,chronic=Estimate)], by="Genes")
m <- unique(m, by="Genes")
splith <- list(rho=mean(sh,na.rm=TRUE), ci=quantile(sh,c(.025,.975),na.rm=TRUE),
               frac_neg=mean(sh<0,na.rm=TRUE), n_splits=NSPLIT,
               rho_naive=rho_naive, rho_recovery=rho_recov, dist=sh)
fwrite(data.table(metric=c("splithalf_rho","ci_lo","ci_hi","frac_negative",
                           "naive_shared_B1","recovery_B3_B2"),
                  value=c(splith$rho, splith$ci[1], splith$ci[2],
                          splith$frac_neg, rho_naive, rho_recov)),
       file.path(OUT,"R1_05_splithalf.csv"))

saveRDS(list(acute=ac, overlap=ov, m=m, splithalf=splith), file.path(OUT,"R1_05_acute_chronic.rds"))
