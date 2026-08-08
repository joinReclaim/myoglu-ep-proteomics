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
m <- merge(ac[,.(Genes,acute=est_B2,acute3=est_B3,acute_32=est_B3-est_B2)],
           CHR[,.(Genes,chronic=Estimate)], by="Genes")
m <- unique(m, by="Genes")
cat("ADVARSEL: B1 inngår i BEGGE estimatene for kontrastene B2-B1 og B3-B1\n",
    "(kronisk = B1-A1, akutt = B2-B1). Målefeil i B1 går inn med motsatt\n",
    "fortegn og induserer negativ korrelasjon uavhengig av biologi. Kun\n",
    "kontrasten B3-B2 deler ingen ledd med kronisk og kan tolkes.\n\n")
for (nm in c("acute","acute3","acute_32")) {
  pe <- cor(m[[nm]], m$chronic); sp <- cor(m[[nm]], m$chronic, method="spearman")
  ct <- cor.test(m[[nm]], m$chronic)
  lab <- switch(nm, acute="B2-B1 (0 min)   [DELT LEDD]",
                    acute3="B3-B1 (2 t)     [DELT LEDD]",
                    acute_32="B3-B2 (2 t vs 0) [RENT]")
  cat(sprintf("  kronisk vs %-28s n=%d: Pearson r=%+.3f (p=%.2g), Spearman rho=%+.3f, samme fortegn %.0f%%\n",
              lab, nrow(m), pe, ct$p.value, sp, 100*mean(m[[nm]]*m$chronic>0)))
}
cat("\nFortolkning: en positiv korrelasjon betyr at akutt- og langtidsprogrammet\ndeler retning; nær null betyr at de er uavhengige. Dette erstatter\nlisteoverlapp, som er et styrkedrevet likhetsmål.\n")
saveRDS(list(acute=ac, overlap=ov, m=m), file.path(OUT,"R1_05_acute_chronic.rds"))
