# =====================================================================
# R1_04_group_time.R — group x time-interaksjon + lipidjustering
# Adresserer R5 punkt 6, R3 punkt 6, R1 punkt 6
#
# NB: gruppene ble dannet på BMI OG glykemi i fellesskap, så en signifikant
# interaksjon kan IKKE tilskrives adipositet alene.
# Krever: R1_00_setup.R
# =====================================================================
suppressPackageStartupMessages({library(data.table); library(lme4); library(lmerTest)})
source(file.path(Sys.getenv("ECV_ROOT", "/Users/sindrle/Research/Prosjekter/ECV"),
                 "review_round_1/01_scripts/R1_config.R"))
ROOT <- ECV_ROOT
OUT  <- file.path(ROOT, "review_round_1/02_output"); setwd(ROOT)
S <- readRDS(file.path(OUT, "R1_00_normalised.rds"))

idx <- which(S$si$SampleCode %in% c("A1","B1"))
M <- S$Xloe[, idx]; s <- S$si[idx]
nA <- rowSums(!is.na(M[, s$SampleCode=="A1"])); nB <- rowSums(!is.na(M[, s$SampleCode=="B1"]))
ok <- nA >= .7*26 & nB >= .7*26; M <- M[ok,]; g <- S$genes[ok]
d0 <- data.table(SampleCode=factor(s$SampleCode,c("A1","B1")), Gr=s$Gr,
                 batch=s$batch, ID2=factor(s$ID2))

res <- rbindlist(lapply(seq_len(nrow(M)), function(i){
  d <- copy(d0); d[, y := M[i,]]; d <- d[!is.na(y)]
  if (uniqueN(d$ID2) < 8 || uniqueN(d$batch) < 2) return(NULL)
  f0 <- tryCatch(suppressWarnings(suppressMessages(lmer(y~SampleCode+Gr+batch+(1|ID2),d,REML=FALSE))),error=function(e)NULL)
  f1 <- tryCatch(suppressWarnings(suppressMessages(lmer(y~SampleCode*Gr+batch+(1|ID2),d,REML=FALSE))),error=function(e)NULL)
  if (is.null(f0)||is.null(f1)) return(NULL)
  a <- anova(f0,f1)
  co <- summary(f1)$coefficients
  nm <- grep("SampleCodeB1:", rownames(co), value=TRUE)
  data.table(Genes=g[i], int_est=if(length(nm)) co[nm[1],"Estimate"] else NA_real_,
             p_int=a$`Pr(>Chisq)`[2], n_ID=uniqueN(d$ID2))
}))
res[, q_int := p.adjust(p_int,"BH")]
setorder(res, p_int)
fwrite(res, file.path(OUT,"R1_04_group_time.csv"))
cat(sprintf("=== GROUP x TIME ===\ntestet %d | nominell p<0.05: %d (forventet ~%.0f) | FDR<0.05: %d\n",
            nrow(res), sum(res$p_int<.05), .05*nrow(res), sum(res$q_int<.05)))
print(head(res[,.(Genes,int_est=round(int_est,3),p_int=signif(p_int,3),q_int=signif(q_int,3))],10))

# power: minste påvisbare interaksjon, 13 vs 13, paret
sdd <- median(apply(M,1,function(x){
  a<-x[s$SampleCode=="A1"];b<-x[s$SampleCode=="B1"];sd(b-a,na.rm=TRUE)}),na.rm=TRUE)
d80 <- (qnorm(.975)+qnorm(.80))*sdd*sqrt(1/13+1/13)
cat(sprintf("\nPOWER: median SD(delta) = %.3f log2. Minste interaksjon påvisbar\nved 80%% power, alfa 0.05 (uten multiplisitet): %.3f log2.\nMed FDR over %d proteiner er reell påvisbar effekt vesentlig større.\n",
            sdd, d80, nrow(res)))

# lipidjustering
cat("\n=== OVERLEVER FENOTYPEASSOSIASJONENE JUSTERING FOR dHDL? ===\n")
PH <- readRDS(file.path(OUT,"R1_02_phenotypes.rds"))
ph <- as.data.table(S$pheno); setnames(ph,"timepoint","Time",skip_absent=TRUE)
pa <- ph[Time=="A1"][order(ID2)]; pb <- ph[Time=="B1"][order(ID2)]
ids <- as.numeric(colnames(PH$dp)); if(all(is.na(ids))) ids <- PH$dph$ID2
pa <- pa[match(PH$dph$ID2,ID2)]; pb <- pb[match(PH$dph$ID2,ID2)]
dHDL <- pb$pl_HDL - pa$pl_HDL; dChol <- pb$pl_TOT.Chol - pa$pl_TOT.Chol
out <- rbindlist(lapply(c("LTBP4","VCL","MEP1B","HSP90B1","TPI1"), function(gn){
  if(!gn %in% rownames(PH$dp)) return(NULL)
  x <- PH$dp[gn,]; y <- PH$dph$VO2max.kg; okk <- !is.na(x)&!is.na(y)&!is.na(dHDL)
  m1 <- summary(lm(rank(y[okk])~rank(x[okk])))
  m2 <- summary(lm(rank(y[okk])~rank(x[okk])+dHDL[okk]+dChol[okk]))
  data.table(Genes=gn, p_ujustert=coef(m1)[2,4], p_lipidjustert=coef(m2)[2,4], n=sum(okk))
}))
print(out[,.(Genes,p_ujustert=signif(p_ujustert,3),p_lipidjustert=signif(p_lipidjustert,3),n)])
