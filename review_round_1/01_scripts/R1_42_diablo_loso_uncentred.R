# =====================================================================
# R1_42_diablo_loso_uncentred.R — leave-one-subject-out for den
# USENTRERTE firblokkmodellen, altså den figur 8 faktisk viser.
#
# Bakgrunn: tallet 91,7 % har sirkulert i dokumentasjonen uten kodemessig
# opphav. R1_13 kryssvaliderer den SENTRERTE modellen (100 %); R1_27, som
# lager figuren, kryssvaliderer ikke i det hele tatt. Med koden offentlig
# deponert kan et tall uten skript ikke stå i tilsvaret.
#
# Blokkene bygges nøyaktig som i R1_27 (usentrert, topp-500 varians per
# blokk, samme filtrering og samme 48 prøver fra 24 deltakere).
# LOSO: hold én DELTAKER ute (begge prøvene), tren på resten, klassifiser
# de to holdt ute. Deltakernivå, ikke prøvenivå — prøvenivå ville lekket
# den andre prøven fra samme person inn i treningssettet.
#
# Krever: R1_00_setup.R
# =====================================================================
suppressPackageStartupMessages({library(data.table); library(mixOmics); library(stringr)})
source(file.path(Sys.getenv("ECV_ROOT", "/Users/sindrle/Research/Prosjekter/ECV"),
                 "review_round_1/01_scripts/R1_config.R"))
ROOT <- ECV_ROOT
OUT  <- file.path(ROOT,"review_round_1/02_output"); setwd(ROOT); set.seed(1)
NCOMP <- 2; KEEPX <- 10

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
keep <- key %in% rownames(SM) & key %in% rownames(AT) & key %in% rownames(olm) & !(key %in% atmiss)
key <- key[keep]; M <- M[,keep]; s <- s[keep]
topv <- function(Z,n=500){ v <- apply(Z,2,var,na.rm=TRUE); Z[, order(v,decreasing=TRUE)[1:min(n,ncol(Z))], drop=FALSE] }
blkf <- function(Z){ Z <- Z[key,,drop=FALSE]; rownames(Z) <- key; topv(Z) }
X <- list(EP = topv(`rownames<-`(t(M), key)), Serum = blkf(olm),
          Muskel = blkf(SM), Fett = blkf(AT))
X <- lapply(X, function(z) z[, apply(z,2,function(x) all(is.finite(x)) && var(x)>0), drop=FALSE])
Y <- factor(s$SampleCode, levels=c("A1","B1"), labels=c("Baseline","Post-training"))
des <- matrix(0.1, length(X), length(X)); diag(des) <- 0
kx <- lapply(X, function(z) rep(min(KEEPX, ncol(z)), NCOMP))
cat(sprintf("prover %d | deltakere %d | blokker %s\n", nrow(X[[1]]), uniqueN(s$ID2),
            paste(sapply(X, ncol), collapse="/")))

ids <- unique(s$ID2)
pred <- rep(NA_character_, length(Y)); predm <- rep(NA_character_, length(Y))
for (u in ids) {
  te <- which(s$ID2 == u); tr <- setdiff(seq_along(Y), te)
  Xtr <- lapply(X, function(z) z[tr,,drop=FALSE])
  Xte <- lapply(X, function(z) z[te,,drop=FALSE])
  kxt <- lapply(Xtr, function(z) rep(min(KEEPX, ncol(z)), NCOMP))
  f <- tryCatch(suppressWarnings(suppressMessages(
         block.splsda(Xtr, Y[tr], ncomp=NCOMP, keepX=kxt, design=des))), error=function(e) NULL)
  if (is.null(f)) next
  pr <- tryCatch(predict(f, newdata=Xte), error=function(e) NULL)
  if (is.null(pr)) next
  # Majoritetsstemme gir NA ved 2-2 mellom fire blokker — det skjedde for 11 av
  # 48 prover. A rapportere treffsikkerhet kun blant de loste ville vaert en
  # selektiv nevner. Vi bruker derfor VEKTET stemme, som mixOmics loser med
  # komponentenes forklaringskraft, og beholder majoritetsstemmen som kontroll.
  wv <- tryCatch(pr$WeightedVote[["max.dist"]][, NCOMP], error=function(e) NULL)
  mv <- pr$MajorityVote[["max.dist"]][, NCOMP]
  pred[te]  <- as.character(if (!is.null(wv)) wv else mv)
  predm[te] <- as.character(mv)
}
ok <- !is.na(pred)
acc <- mean(pred[ok] == as.character(Y)[ok])
okm <- !is.na(predm)
cat(sprintf("\n--- majoritetsstemme (kontroll) ---\n  lost: %d av %d | korrekt blant loste: %.1f%% | korrekt av ALLE 48: %.1f%%\n",
            sum(okm), length(Y), 100*mean(predm[okm]==as.character(Y)[okm]),
            100*sum(predm[okm]==as.character(Y)[okm])/length(Y)))
bt  <- binom.test(sum(pred[ok]==as.character(Y)[ok]), sum(ok), 0.5)
cat(sprintf("\n=== LOSO, USENTRERT, fire blokker ===\n"))
cat(sprintf("  klassifisert: %d av %d prover (%d deltakere)\n", sum(ok), length(Y), uniqueN(s$ID2)))
cat(sprintf("  korrekt     : %.1f%% (%d av %d)\n", 100*acc, sum(pred[ok]==as.character(Y)[ok]), sum(ok)))
cat(sprintf("  binomialtest mot 50%%: p = %.3g\n", bt$p.value))
print(table(sant=as.character(Y)[ok], predikert=pred[ok]))

fwrite(data.table(model="four blocks, uncentred, LOSO by participant",
                  n_samples=sum(ok), n_participants=uniqueN(s$ID2),
                  n_correct=sum(pred[ok]==as.character(Y)[ok]),
                  accuracy=round(acc,4), binom_p=bt$p.value),
       file.path(OUT,"R1_42_loso_uncentred.csv"))
fwrite(data.table(sample=key[ok], ID2=s$ID2[ok], true=as.character(Y)[ok], predicted=pred[ok]),
       file.path(OUT,"R1_42_loso_predictions.csv"))
cat("\nskrevet: R1_42_loso_uncentred.csv\n")
