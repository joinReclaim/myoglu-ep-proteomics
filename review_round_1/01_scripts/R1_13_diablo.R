# =====================================================================
# R1_13_diablo.R — DIABLO med kryssvalidering og stabilitetsvurdering
# Adresserer R5 punkt 8, R3 punkt 7, R2 concern 6
#
# Innvendingene mot 202_MixOmics_DIABLO.R:
#   - separasjon av pre/post etter SUPERVISED feature-seleksjon på
#     subjekt-sentrerte data er forventet, ikke et funn
#   - top-500 per blokk forhåndsvalgt på varians, keepX=10 fast,
#     ingen perf()/tune(), cutoff 0.6 satt manuelt
#   - hele manglende fettvevsprøver ble imputert med gruppe x tid-
#     gjennomsnitt, altså FABRIKERTE prøver. Gjøres ikke her.
#
# Her: kryssvalidert klassifikasjonsfeil mot permutert nullfordeling,
# og stabilitet av feature-seleksjon på tvers av folds.
# Krever: R1_00_setup.R
# =====================================================================
suppressPackageStartupMessages({library(data.table); library(mixOmics); library(stringr)})
source(file.path(Sys.getenv("ECV_ROOT", "/Users/sindrle/Research/Prosjekter/ECV"),
                 "review_round_1/01_scripts/R1_config.R"))
ROOT <- ECV_ROOT
OUT  <- file.path(ROOT,"review_round_1/02_output"); setwd(ROOT); set.seed(1)
source(file.path(ROOT,"review_round_1/01_scripts/R1_theme.R"))
S <- readRDS(file.path(OUT,"R1_00_normalised.rds"))

TOPN <- 500; NCOMP <- 2; KEEPX <- 10; FOLDS <- 5; NREP <- 10; NPERM <- 30

# ---- blokker: EP (loess), muskel, fett. INGEN imputering av hele prøver ----
i <- which(S$si$SampleCode %in% c("A1","B1")); M <- S$Xloe[,i]; s <- S$si[i]
rownames(M) <- make.unique(S$genes)
M <- M[rowSums(!is.na(M)) == ncol(M), ]
load(RNASEQ_RDA)
mk <- function(mat){
  d <- as.data.frame(mat); sid <- rownames(d)
  tm <- str_extract(sid,"A1|B1"); raw <- as.numeric(str_extract(sid,"\\d+"))
  id <- ifelse(tm=="A1"&raw>=100, raw-100, ifelse(tm=="B1"&raw>=200, raw-200, raw))
  ok <- !is.na(tm)&!is.na(id); d <- d[ok,]; rownames(d) <- paste0(id[ok],"_",tm[ok])
  log2(as.matrix(d)+1) }
SM <- mk(fpkm.sm.fix); AT <- mk(fpkm.at.fix)

# ---- serum/Olink-blokk (fjerde blokk, som i originalen) ----
# NB: DIABLO kobler blokker på PRØVER, ikke på features. Serumblokken
# trenger derfor ingen proteinoverlapp med EP-blokken.
ol <- fread(file.path(OLINK_DIR,"MyoGlu_Olink_A1B1_NPX.csv"), sep=";")
olm <- as.matrix(ol[, !c("ID","Time","group"), with=FALSE])
storage.mode(olm) <- "double"
rownames(olm) <- paste0(ol$ID, "_", ol$Time)
olm <- olm[, colSums(!is.na(olm)) == nrow(olm), drop=FALSE]
key <- paste0(s$ID2,"_",s$SampleCode)
keep <- key %in% rownames(SM) & key %in% rownames(AT) & key %in% rownames(olm)
# 4 fettvevsprøver mangler HELT. Originalskriptet fylte dem med gruppe x tid-
# gjennomsnitt, altså fabrikerte prøver. Her droppes de i stedet.
at_missing <- rownames(AT)[rowSums(!is.na(AT)) == 0]
cat(sprintf("fettvevsprøver som mangler helt og droppes: %d (%s)\n",
            length(at_missing), paste(at_missing, collapse=", ")))
keep <- keep & !(key %in% at_missing)
key <- key[keep]; M <- M[,keep]; s <- s[keep]
cat(sprintf("komplette prøver i alle tre blokker: %d (%s)\n", length(key),
            paste(names(table(s$SampleCode)), table(s$SampleCode), collapse=", ")))

# subjekt-sentrering (som originalen) — merk at dette i seg selv gjør pre/post
# lettere å skille, og er grunnen til at separasjon ikke er et funn
ctr <- function(X, id){ for(u in unique(id)){ j <- id==u
  X[j,] <- sweep(X[j,,drop=FALSE],2,colMeans(X[j,,drop=FALSE],na.rm=TRUE)) }; X }
topv <- function(X,n){ v <- apply(X,2,var,na.rm=TRUE); X[, order(v,decreasing=TRUE)[1:min(n,ncol(X))]] }
EPm <- t(M); rownames(EPm) <- key          # blokkene MÅ ha identiske radnavn
SMm <- SM[key, , drop=FALSE]; rownames(SMm) <- key
ATm <- AT[key, , drop=FALSE]; rownames(ATm) <- key
OLm <- olm[key, , drop=FALSE]; rownames(OLm) <- key
X <- list(EP = ctr(topv(EPm, TOPN), s$ID2),
          Serum = ctr(topv(OLm, TOPN), s$ID2),
          Muskel = ctr(topv(SMm, TOPN), s$ID2),
          Fett = ctr(topv(ATm, TOPN), s$ID2))
X <- lapply(X, function(z) z[, apply(z,2,function(x) all(is.finite(x)) && var(x) > 0), drop=FALSE])
cat("blokkstørrelser etter filtrering: ",
    paste(sprintf("%s %dx%d", names(X), sapply(X,nrow), sapply(X,ncol)), collapse=" | "), "\n")
X <- lapply(X, function(z){ rownames(z) <- key; z })
Y <- factor(s$SampleCode, c("A1","B1")); names(Y) <- key
stopifnot(all(sapply(X, function(z) identical(rownames(z), key))))
des <- matrix(0.1, length(X), length(X)); diag(des) <- 0
kx <- lapply(X, function(z) rep(min(KEEPX, ncol(z)), NCOMP))

fitm <- block.splsda(X, Y, ncomp=NCOMP, keepX=kx, design=des)
pf <- perf(fitm, validation="Mfold", folds=FOLDS, nrepeat=NREP, progressBar=FALSE)
ber <- pf$WeightedVote.error.rate$max.dist["Overall.BER", ]
cat(sprintf("\n=== KRYSSVALIDERT KLASSIFIKASJONSFEIL (BER, %d-fold x %d) ===\n", FOLDS, NREP))
for (k in seq_along(ber)) cat(sprintf("  komponent %d: BER = %.3f\n", k, ber[k]))

cat(sprintf("\n=== PERMUTERT NULL (%d omstokkinger av tidspunkt-etiketten) ===\n", NPERM))
nullber <- sapply(seq_len(NPERM), function(z){
  Yp <- factor(sample(as.character(Y)), c("A1","B1"))
  f <- tryCatch(block.splsda(X, Yp, ncomp=NCOMP, keepX=kx, design=des), error=function(e)NULL)
  if(is.null(f)) return(NA_real_)
  p <- tryCatch(perf(f, validation="Mfold", folds=FOLDS, nrepeat=2, progressBar=FALSE), error=function(e)NULL)
  if(is.null(p)) return(NA_real_)
  min(p$WeightedVote.error.rate$max.dist["Overall.BER", ], na.rm=TRUE) })
nullber <- nullber[!is.na(nullber)]
obs <- min(ber, na.rm=TRUE)
cat(sprintf("  observert beste BER: %.3f | null: median %.3f [%.3f, %.3f]\n",
            obs, median(nullber), quantile(nullber,.05), quantile(nullber,.95)))
cat(sprintf("  permutasjons-p = %.3f\n", (1+sum(nullber<=obs))/(length(nullber)+1)))

saveRDS(list(perf=pf, ber=ber, nullber=nullber, X=X, Y=Y, key=key, s=s),
        file.path(OUT,"R1_13_diablo.rds"))

# ---- KORREKT KRYSSVALIDERING: hold ute HELE DELTAKERE, ikke enkeltprøver ----
# Subjekt-sentrering med to tidspunkt gjør prøvene til eksakte speilbilder
# (x_A1 = -x_B1). Splittes de på prøvenivå havner en deltakers to prøver i
# hver sin fold, og testprøven er da minus treningsprøven — perfekt
# klassifisering er garantert av transformasjonen, ikke av biologi.
subj <- s$ID2
cat("\n=== LEAVE-ONE-SUBJECT-OUT (begge prøver holdes ute samlet) ===\n")
pred <- sapply(unique(subj), function(u){
  tr <- subj != u; te <- subj == u
  Xtr <- lapply(X, function(z) z[tr,,drop=FALSE]); Xte <- lapply(X, function(z) z[te,,drop=FALSE])
  f <- tryCatch(block.splsda(Xtr, Y[tr], ncomp=NCOMP, keepX=kx, design=des), error=function(e) NULL)
  if (is.null(f)) return(NA_real_)
  pr <- tryCatch(predict(f, newdata=Xte), error=function(e) NULL)
  if (is.null(pr)) return(NA_real_)
  cl <- pr$WeightedVote$max.dist[, NCOMP]
  mean(as.character(cl) == as.character(Y[te]), na.rm=TRUE) })
pred <- pred[!is.na(pred)]
cat(sprintf("  korrekt klassifisert: %.1f%% (%d deltakere)\n", 100*mean(pred), length(pred)))
cat(sprintf("  binomialtest mot 50%%: p = %.3g\n",
            binom.test(round(sum(pred)*2), 2*length(pred), 0.5)$p.value))

cat("\n=== STABILITET AV FEATURE-SELEKSJON ===\n")
st_all <- tryCatch(pf$features$stable, error=function(e) NULL)
if (is.null(st_all)) cat("  ikke tilgjengelig i denne mixOmics-versjonen\n") else {
  flat <- function(z) if (is.list(z)) unlist(lapply(z, flat), use.names=TRUE) else z
  for (b in names(X)) {
    v <- tryCatch(flat(st_all[[b]]), error=function(e) NULL)
    if (is.null(v)) v <- tryCatch(flat(lapply(st_all, `[[`, b)), error=function(e) NULL)
    if (is.null(v) || !length(v)) { cat(sprintf("  %s: ikke tilgjengelig\n", b)); next }
    v <- sort(v[is.finite(v)], decreasing=TRUE)
    cat(sprintf("  %-7s %d features valgt >=50%% av folds, %d >=80%%\n",
                b, sum(v>=.5), sum(v>=.8))) } }

d <- data.table(BER=nullber)
p <- ggplot(d, aes(BER)) +
  geom_histogram(bins=20, fill="grey70", colour=NA) +
  geom_vline(xintercept=min(ber, na.rm=TRUE), colour=cols[1], linewidth=1) +
  labs(x="Balansert klassifikasjonsfeil (kryssvalidert)", y="Antall permutasjoner") +
  theme_ecv()
save_ecv("R1_13_DIABLO_permutation.pdf", p, width=10, height=8)
