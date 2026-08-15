# =====================================================================
# R1_43_candidate_coclustering.R — holder firergruppen sammen etter
# reanalysen?
#
# BAKGRUNN OG HVORFOR DETTE ER RIKTIG TEST.
# De fire kandidatene (LTBP4, VCL, MEP1B, HSP90B1) ble IKKE valgt som de
# fire sterkeste mot dVO2max. De ble valgt fordi de klumpet seg sammen
# visuelt i protein x fenotype-heatmappet (manuskriptets figur 2, nå
# figur 4). Det er en annen seleksjonsregel, og den har konsekvenser:
#
#   - Permutasjonstesten i R1_02 velger de fire hoyeste |rho| mot dVO2max
#     på nytt i hver permutasjon. Den tester dermed en prosedyre som
#     aldri ble brukt, og gir p = 0.23 for noe ingen har påstått.
#   - Rangeringen mot dVO2max (LTBP4 nr 4, HSP90B1 15, MEP1B 16, VCL 93
#     av 1539) er heller ikke det relevante kriteriet.
#
# Riktig test av den FAKTISKE seleksjonen: er de fire mer like hverandre
# i sitt korrelasjonsmonster over alle 14 fenotyper enn en tilfeldig
# firergruppe? Det er «klumper de seg sammen i heatmappet», tallfestet.
#
# MERK om sirkularitet: de ble valgt på heatmappet fra den OPPRINNELIGE
# analysen og testes her på det REANALYSERTE. Samme deltakere og samme
# proteiner, men annen preprosessering — så dette er robusthet mot
# normaliseringsendringen, ikke uavhengig replikasjon. Det skal sies.
#
# Krever: R1_23_fig4_heatmap.R (R1_23_fig4_matrix.csv)
# =====================================================================
suppressPackageStartupMessages({library(data.table)})
source(file.path(Sys.getenv("ECV_ROOT", "/Users/sindrle/Research/Prosjekter/ECV"),
                 "review_round_1/01_scripts/R1_config.R"))
ROOT <- ECV_ROOT
OUT  <- file.path(ROOT, "review_round_1/02_output"); setwd(ROOT); set.seed(1)
NPERM <- 10000
CAND  <- c("LTBP4","VCL","MEP1B","HSP90B1")

m  <- fread(file.path(OUT, "R1_23_fig4_matrix.csv"))
M  <- as.matrix(m[, -1]); rownames(M) <- m$phenotype
cand <- intersect(CAND, colnames(M))
stopifnot(length(cand) == 4)
cat(sprintf("heatmap-matrise: %d fenotyper x %d proteiner\n", nrow(M), ncol(M)))

C   <- cor(M, use = "pairwise.complete.obs")
pr  <- combn(cand, 2)
obs <- mean(apply(pr, 2, function(p) C[p[1], p[2]]))

nullv <- replicate(NPERM, {
  g <- sample(colnames(M), length(cand))
  mean(apply(combn(g, 2), 2, function(p) C[p[1], p[2]]))
})
pval <- (1 + sum(nullv >= obs)) / (NPERM + 1)

cat(sprintf("\n=== KOKLYNGING OVER 14 FENOTYPER ===\n"))
cat(sprintf("  observert profil-likhet, de fire : %+.3f\n", obs))
cat(sprintf("  tilfeldige firergrupper          : median %+.3f | 95%%-persentil %+.3f\n",
            median(nullv), quantile(nullv, .95)))
cat(sprintf("  permutasjons-p (%d trekninger)  : %.4f\n", NPERM, pval))
cat("\nparvise likheter:\n")
for (i in seq_len(ncol(pr)))
  cat(sprintf("  %-8s - %-8s %+.3f\n", pr[1,i], pr[2,i], C[pr[1,i], pr[2,i]]))

fwrite(data.table(metric = c("observed_mean_similarity","null_median","null_p95",
                             "perm_p","n_perm","n_phenotypes","n_proteins"),
                  value = c(round(obs,4), round(median(nullv),4),
                            round(unname(quantile(nullv,.95)),4), pval,
                            NPERM, nrow(M), ncol(M))),
       file.path(OUT, "R1_43_candidate_coclustering.csv"))
fwrite(data.table(protein_1 = pr[1,], protein_2 = pr[2,],
                  similarity = round(apply(pr, 2, function(p) C[p[1],p[2]]), 4)),
       file.path(OUT, "R1_43_candidate_pairs.csv"))
cat("\nskrevet: R1_43_candidate_coclustering.csv\n")
