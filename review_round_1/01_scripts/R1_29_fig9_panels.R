# =====================================================================
# R1_29_fig9_panels.R — manuskriptets FIGUR 9 (Figures/Figure7.tif)
# (blir figur 8 etter omnummerering)
# Erstatter 101_ECV_acute_responses.R:194-244 og 363-394
#
# Panel A: fire akutte forløpsmønstre. Uendret stil: tynne linjer per protein,
#   tykk klyngegjennomsnitt, facet_wrap, x-akse B / 0' / 2t, Tableau-farger.
#   ENDRET: k=4 beholdes for TOLKBARHET, ikke fordi silhouette peker dit —
#   profilen maksimeres ved k=2 (0,705) eller k=3. Må stå i legenden (R4).
# Panel B: Venn. ENDRET: begrenset til proteiner testet i BEGGE analysene, og
#   prosenter relative til hvert sett (ikke union — union er samme nevnerfeil
#   R3 påpeker for det opprinnelige 21,6 %).
# Panel C: NYTT. Akutt mot kronisk effekt for de overlappende (R5 punkt 5).
#   Kun de overlappende vises: hele skyen ville overdrevet anti-korrelasjonen,
#   siden B2−B1 og B1−A1 deler B1. Statistikken i legenden skal være
#   splitt-halv-estimatet (rho = −0,42), ikke noe avlest fra dette panelet.
# Krever: R1_00_setup.R, R1_01_longterm.R, R1_05_acute_chronic.R
# =====================================================================
suppressPackageStartupMessages({library(data.table); library(cluster); library(ggrepel)})
source(file.path(Sys.getenv("ECV_ROOT", "/Users/sindrle/Research/Prosjekter/ECV"),
                 "review_round_1/01_scripts/R1_config.R"))
ROOT <- ECV_ROOT
OUT  <- file.path(ROOT,"review_round_1/02_output"); setwd(ROOT); set.seed(1)
source(file.path(ROOT,"review_round_1/01_scripts/R1_theme.R"))
S  <- readRDS(file.path(OUT,"R1_00_normalised.rds"))
AC <- readRDS(file.path(OUT,"R1_05_acute_chronic.rds"))$acute
CH <- readRDS(file.path(OUT,"R1_01_longterm.rds"))$PRIM
TP <- c("B1","B2","B3"); LAB <- c("B", "0'", "2h"); K <- 4

# ---- panel A ----
i <- which(S$si$SampleCode %in% TP); M <- S$Xloe[, i]; s <- S$si[i]
rownames(M) <- make.unique(S$genes)
g <- intersect(unique(AC[q_time < 0.25, Genes]), rownames(M))
tm <- t(sapply(g, function(x) sapply(TP, function(t) mean(M[x, s$SampleCode==t], na.rm=TRUE))))
tm <- tm[apply(tm,1,function(x) all(is.finite(x))), , drop=FALSE]
tm <- t(apply(tm, 1, function(x) (x-mean(x))/sd(x)))
hc <- hclust(dist(tm), method="ward.D2"); cl <- cutree(hc, K)
sil <- sapply(2:8, function(k) mean(silhouette(cutree(hc,k), dist(tm))[,3]))
cat(sprintf("panel A: %d proteiner | klynger: %s\n", nrow(tm),
            paste(sprintf("C%d n=%d",1:K,as.integer(table(cl))), collapse=", ")))
cat(sprintf("  silhouette k=2..8: %s (k=4 gir %.3f)\n",
            paste(sprintf("%.2f",sil),collapse=" "), sil[3]))
pl <- rbindlist(lapply(1:K, function(k) rbindlist(lapply(1:3, function(j)
  data.table(Gene=rownames(tm)[cl==k], Tid=j, z=tm[cl==k,j],
             Cluster=sprintf("C%d (n=%d)", k, sum(cl==k)), ck=paste0("C",k))))))
pA <- ggplot(pl, aes(Tid, z, group=Gene, colour=ck)) +
  geom_line(alpha=.18, linewidth=.4) +
  stat_summary(aes(group=1), fun=mean, geom="line", linewidth=1.5) +
  stat_summary(aes(group=1), fun=mean, geom="point", size=2.8) +
  facet_wrap(~Cluster, axes="all") +
  scale_colour_manual(values = PAL$cluster) +
  scale_x_continuous(breaks=1:3, labels=LAB) +
  labs(x=NULL, y="Z-scored mean protein abundance") +
  theme_ecv() + theme(legend.position="none")
save_ecv("R1_29_fig9A_patterns.pdf", pA, width=15, height=12, unit="cm")

# ---- panel B ----
J <- merge(unique(AC[,.(Genes,akutt=est_B2,q_a=q_time)],by="Genes"),
           unique(CH[,.(Genes,kronisk=Estimate,q_c=q)],by="Genes"), by="Genes")
# q<0.25 gjennomgående i figuren: panel A kan ikke klynge 36 proteiner, og
# «akutte respondere» kan ikke bety to ulike tall i samme figur.
# Primærtallet (q<0.05: 125/35/12, 9.6 % av kronisk) oppgis i teksten.
THR_FIG <- 0.25
a <- J[q_a<THR_FIG, Genes]; c_ <- J[q_c<THR_FIG, Genes]; ov <- intersect(a,c_)
strict <- intersect(J[q_a<.05,Genes], J[q_c<.05,Genes])
nA_<-length(a); nC<-length(c_); nO<-length(ov)
ex <- nA_*nC/nrow(J)
ft <- fisher.test(matrix(c(nO,nC-nO,nA_-nO,nrow(J)-nA_-nC+nO),2))
cat(sprintf("\npanel B: testet i begge %d | kronisk %d | akutt %d | overlapp %d\n",
            nrow(J), nC, nA_, nO))
cat(sprintf("  %.1f%% av kroniske | %.1f%% av akutte | forventet %.1f -> %.1fx (p=%.3g)\n",
            100*nO/nC, 100*nO/nA_, ex, nO/ex, ft$p.value))
circ <- function(x0,r=1.05,n=200){ t<-seq(0,2*pi,length.out=n); data.table(x=x0+r*cos(t), y=r*sin(t)) }
cd <- rbind(circ(-0.62)[,set:="Chronic"], circ(0.62)[,set:="Acute"])
lb <- data.table(x=c(-1.15,0,1.15,-1.15,1.15), y=c(0,0,0,1.35,1.35),
  # overlappet har TO gyldige nevnere — begge vises, rekkefølgen defineres
  # i legenden (kronisk / akutt). Ellers går ikke høyre sirkel opp til 100 %.
  lab=c(sprintf("%d\n(%.1f%%)",nC-nO,100*(nC-nO)/nC),
        sprintf("%d\n(%.1f%% / %.1f%%)",nO,100*nO/nC,100*nO/nA_),
        sprintf("%d\n(%.1f%%)",nA_-nO,100*(nA_-nO)/nA_), "Chronic","Acute"), sz=c(4,3.5,4,4.6,4.6))
pB <- ggplot() +
  geom_polygon(data=cd, aes(x,y,fill=set,group=set), alpha=.55, colour="black", linewidth=.6) +
  scale_fill_manual(values=setNames(PAL$venn, c("Chronic","Acute"))) +
  geom_text(data=lb, aes(x,y,label=lab,size=sz), family=fnt) +
  scale_size_identity() + coord_equal() +
  theme_void() + theme(legend.position="none", text=element_text(family=fnt))
save_ecv("R1_29_fig9B_venn.pdf", pB, width=9, height=8, unit="cm")

# ---- panel C ----
d <- J[Genes %in% ov]
cat(sprintf("panel C: %d overlappende, %d med motsatt fortegn\n",
            nrow(d), sum(d$akutt*d$kronisk<0)))
d[, strikt := Genes %in% strict]
cat(sprintf("  av disse %d også signifikante ved q<0.05\n", sum(d$strikt)))
pC <- ggplot(d, aes(kronisk, akutt)) +
  geom_hline(yintercept=0, linetype="dotted") + geom_vline(xintercept=0, linetype="dotted") +
  geom_point(data=d[!(strikt)], colour="grey72", size=1.6) +
  geom_point(data=d[(strikt)], colour=cols[1], size=2.4) +
  geom_text_repel(data=d[(strikt)], aes(label=Genes), size=2.9, family=fnt,
                  min.segment.length=0, segment.size=.2, seed=1, max.overlaps=Inf) +
  scale_x_continuous(expand = expansion(mult = 0.18)) +
  labs(x=expression(paste("Chronic training effect, ",Delta,"log"[2],"(int)")),
       y=expression(paste("Acute exercise effect, ",Delta,"log"[2],"(int)"))) +
  theme_ecv()
save_ecv("R1_29_fig9C_acute_vs_chronic.pdf", pC, width=14, height=9, unit="cm")
fwrite(d[order(-abs(kronisk))], file.path(OUT,"R1_29_fig9_overlap.csv"))
