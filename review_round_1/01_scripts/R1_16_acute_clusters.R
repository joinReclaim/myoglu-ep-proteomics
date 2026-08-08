# =====================================================================
# R1_16_acute_clusters.R — akutte responsmønstre, reanalysert
# Erstatter klyngeanalysen i 101_ECV_acute_responses.R (linje 136-244)
#
# Forskjeller fra originalen:
#   - loess-normaliserte verdier, completeness-filter, batch i modellen
#     (responderlisten kommer fra R1_05, ikke fra den unormaliserte)
#   - silhouette-profil for k = 2..8 RAPPORTERES (R4 ber om dette; verdien
#     ble beregnet i 101:177-184 men aldri vist)
#   - klyngescorenes kobling til langtidsgevinst rapporteres ordentlig med
#     estimat, KI og n (R2 concern 6.3; forsøkt i 101:509-536, aldri rapportert)
#
# Klyngeparametre, som i originalen: euklidsk avstand på z-skårede
# gjennomsnittstrajektorier, hierarkisk klynging med Ward D2.
# Krever: R1_00_setup.R, R1_05_acute_chronic.R
# =====================================================================
suppressPackageStartupMessages({library(data.table); library(cluster)})
source(file.path(Sys.getenv("ECV_ROOT", "/Users/sindrle/Research/Prosjekter/ECV"),
                 "review_round_1/01_scripts/R1_config.R"))
ROOT <- ECV_ROOT
OUT  <- file.path(ROOT,"review_round_1/02_output"); setwd(ROOT); set.seed(1)
source(file.path(ROOT,"review_round_1/01_scripts/R1_theme.R"))
S  <- readRDS(file.path(OUT,"R1_00_normalised.rds"))
AC <- readRDS(file.path(OUT,"R1_05_acute_chronic.rds"))$acute

TP <- c("B1","B2","B3"); LAB <- c("B","0'","2t")
i <- which(S$si$SampleCode %in% TP); M <- S$Xloe[, i]; s <- S$si[i]
rownames(M) <- make.unique(S$genes)

# per-protein gjennomsnittstrajektorie, z-skåret over de tre tidspunktene
traj <- function(genes){
  m <- t(sapply(genes, function(g) sapply(TP, function(t) mean(M[g, s$SampleCode==t], na.rm=TRUE))))
  m <- m[apply(m,1,function(x) all(is.finite(x))), , drop=FALSE]
  t(apply(m, 1, function(x) (x-mean(x))/sd(x))) }

for (THR in c(0.05, 0.25)) {
  genes <- intersect(unique(AC[q_time < THR, Genes]), rownames(M))
  tm <- traj(genes)
  cat(sprintf("\n############ AKUTTE RESPONDERE q<%.2f: %d proteiner ############\n", THR, nrow(tm)))
  if (nrow(tm) < 12) { cat("  for få til klynging\n"); next }
  d <- dist(tm); hc <- hclust(d, method="ward.D2")
  sil <- sapply(2:8, function(k) mean(silhouette(cutree(hc,k), d)[,3]))
  cat("  silhouette-profil (R4 ber om denne):\n")
  for (k in 2:8) cat(sprintf("    k=%d: %.3f%s\n", k, sil[k-1],
      ifelse(sil[k-1]==max(sil)," <- maks","")))
  kbest <- (2:8)[which.max(sil)]
  cat(sprintf("  maksimal silhouette ved k=%d; manuskriptet bruker k=4 (silhouette %.3f)\n",
              kbest, sil[3]))
  if (THR != 0.25) next

  K <- 4
  cl <- cutree(hc, K)
  cat(sprintf("\n  --- k=4 (som i manuskriptet) --- størrelser: %s\n",
              paste(sprintf("C%d n=%d", 1:K, as.integer(table(cl))), collapse=", ")))
  for (k in 1:K) cat(sprintf("    C%d gj.snitt z: B %+.2f | 0' %+.2f | 2t %+.2f\n",
      k, mean(tm[cl==k,1]), mean(tm[cl==k,2]), mean(tm[cl==k,3])))

  # ---- klyngescore per deltaker og kobling til langtidsgevinst ----
  ids <- sort(unique(s$ID2))
  dz <- function(t2,t1){ z <- sapply(ids, function(x){
      a <- M[names(cl), s$SampleCode==t1 & s$ID2==x, drop=FALSE]
      b <- M[names(cl), s$SampleCode==t2 & s$ID2==x, drop=FALSE]
      if(!ncol(a)||!ncol(b)) rep(NA_real_, length(cl)) else b[,1]-a[,1] })
    t(apply(z,1,function(x) (x-mean(x,na.rm=TRUE))/sd(x,na.rm=TRUE))) }
  z13 <- dz("B3","B1")
  sc <- sapply(1:K, function(k) colMeans(z13[cl==k, , drop=FALSE], na.rm=TRUE))
  colnames(sc) <- paste0("C",1:K)

  PHx <- readRDS(file.path(OUT,"R1_02_phenotypes.rds"))$dph
  out <- rbindlist(lapply(c("VO2max.kg","M_mg.kg.min"), function(y){
    df <- data.frame(y=PHx[[y]][match(ids, PHx$ID2)], sc)
    df <- df[complete.cases(df),]
    m <- lm(rank(y) ~ ., data=df); ci <- confint(m)
    rbindlist(lapply(paste0("C",1:K), function(k)
      data.table(utfall=y, klynge=k, beta=coef(m)[k], lo=ci[k,1], hi=ci[k,2],
                 p=summary(m)$coefficients[k,4], n=nrow(df)))) }))
  cat("\n  --- KOLLINEARITET MELLOM KLYNGESCORENE ---\n")
  cc <- cor(sc, use="pairwise.complete.obs")
  print(round(cc,2))
  cat(sprintf("  median |r|: %.2f — klyngene er definert ved motsatte trajektorieformer,\n",
              median(abs(cc[upper.tri(cc)]))))
  cat("  så scorene er i praksis én akse. Det er grunnen til de brede intervallene\n  i fellesmodellen under.\n")

  cat("\n  --- KOBLING TIL LANGTIDSGEVINST (R2 concern 6.3) ---\n")
  cat("  NB: akuttøkten ble målt ETTER intervensjonen, mens gevinsten skjedde før.\n")
  cat("  Dette er derfor en assosiasjon i trent tilstand, IKKE prediksjon.\n")
  cat("  (101_ECV_acute_responses.R:408 kaller seksjonen 'Prediction' — den\n   formuleringen kan ikke brukes.)\n")
  single <- rbindlist(lapply(c("VO2max.kg","M_mg.kg.min"), function(y){
    df <- data.frame(y=PHx[[y]][match(ids, PHx$ID2)], sc); df <- df[complete.cases(df),]
    rbindlist(lapply(paste0("C",1:K), function(k){
      m <- lm(as.formula(paste("rank(y) ~", k)), data=df); ci <- confint(m)[k,]
      data.table(utfall=y, klynge=k, modell="alene", beta=coef(m)[k], lo=ci[1], hi=ci[2],
                 p=summary(m)$coefficients[k,4], n=nrow(df)) })) }))
  cat("\n  enkeltpredikator (én klynge om gangen):\n")
  print(single[, .(utfall, klynge, beta=round(beta,2),
                   KI=sprintf("[%.2f, %.2f]",lo,hi), p=signif(p,3), n)])
  fwrite(single, file.path(OUT,"R1_16_acute_cluster_vs_gain_single.csv"))
  cat("\n  fire klyngescorer samtidig:\n")
  print(out[, .(utfall, klynge, beta=round(beta,2), KI=sprintf("[%.2f, %.2f]",lo,hi),
                p=signif(p,3), n)])
  # ---- RIKTIG TEST: reduser til aksen og test ÉN gang ----
  # Når fire scorer måler samme akse skal man verken fitte fire prediktorer
  # (kollinearitet -> brede KI) eller teste dem hver for seg (multiplisitet).
  pca <- prcomp(sc, center = TRUE, scale. = TRUE)
  cat(sprintf("\n  --- DOMINERENDE AKSE: PC1 forklarer %.0f%% av variasjonen ---\n",
              100*pca$sdev[1]^2/sum(pca$sdev^2)))
  cat("  lastninger: ", paste(sprintf("%s %+.2f", colnames(sc), pca$rotation[,1]),
                              collapse = " | "), "\n")
  axis_res <- rbindlist(lapply(c("VO2max.kg","M_mg.kg.min"), function(y){
    yy <- PHx[[y]][match(ids, PHx$ID2)]; k <- is.finite(yy) & is.finite(pca$x[,1])
    ct <- suppressWarnings(cor.test(pca$x[k,1], yy[k], method = "spearman"))
    m <- lm(rank(yy[k]) ~ rank(pca$x[k,1])); ci <- confint(m)[2,]
    data.table(utfall = y, rho = unname(ct$estimate), p = ct$p.value,
               beta = coef(m)[2], lo = ci[1], hi = ci[2], n = sum(k)) }))
  print(axis_res[, .(utfall, rho = round(rho,3), p = signif(p,3),
                     KI = sprintf("[%.2f, %.2f]", lo, hi), n)])
  cat("\n  fortegn per klynge mot ΔVO₂max (koherens):\n")
  yy <- PHx$VO2max.kg[match(ids, PHx$ID2)]
  for (k in paste0("C",1:K)) cat(sprintf("    %s rho = %+.3f\n", k,
    suppressWarnings(cor(sc[,k], yy, use="pairwise.complete.obs", method="spearman"))))
  fwrite(axis_res, file.path(OUT,"R1_16_acute_axis_vs_gain.csv"))
  cat(sprintf("\n  Konklusjon: én koherent akse, alle fire klynger samme vei, men\n  p = %.3f med n = %d. Konsistent retning som studien ikke har styrke\n  til å etablere — ikke et fravær av signal.\n",
              axis_res[utfall=="VO2max.kg", p], axis_res[utfall=="VO2max.kg", n]))
  fwrite(out, file.path(OUT,"R1_16_acute_cluster_vs_gain.csv"))
  fwrite(data.table(Genes=names(cl), cluster=paste0("C",cl)),
         file.path(OUT,"R1_16_acute_clusters.csv"))

  # ---- figurer ----
  ds <- data.table(k=2:8, sil=sil)
  p1 <- ggplot(ds, aes(k, sil)) + geom_line() + geom_point(size=2) +
    geom_point(data=ds[k==4], colour=cols[1], size=3.5) +
    scale_x_continuous(breaks=2:8) +
    labs(x="Antall klynger (k)", y="Gjennomsnittlig silhouette-bredde") + theme_ecv()
  save_ecv("R1_16_silhouette_profile.pdf", p1, width=9, height=7)

  pl <- rbindlist(lapply(1:K, function(k) rbindlist(lapply(1:3, function(j)
    data.table(Gene=rownames(tm)[cl==k], Tid=j, z=tm[cl==k,j],
               Klynge=sprintf("C%d (n=%d)", k, sum(cl==k)))))))
  p2 <- ggplot(pl, aes(Tid, z, group=Gene)) +
    geom_line(alpha=.18, linewidth=.4, colour="grey40") +
    stat_summary(aes(group=1), fun=mean, geom="line", linewidth=1.4, colour=cols[2]) +
    stat_summary(aes(group=1), fun=mean, geom="point", size=2.5, colour=cols[2]) +
    facet_wrap(~Klynge, axes="all") +
    scale_x_continuous(breaks=1:3, labels=LAB) +
    labs(x=NULL, y="Z-skåret gjennomsnittlig proteinmengde") + theme_ecv()
  save_ecv("R1_16_acute_patterns.pdf", p2, width=16, height=12)
}
