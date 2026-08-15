# =====================================================================
# R1_40_vesiclepedia_plasma.R — Vesiclepedia-bakgrunnsraten, PLASMA-begrenset
#
# Erstatter R1_38. Verifiseringen viste at R1_38 annoterte mot ALLE humane
# Vesiclepedia-oppføringer uansett kilde. Manuskriptets påstand er derimot
# spesifikt om PLASMA: «278 of the exercise-responsive EP proteins have
# previously been identified in human plasma extracellular vesicles.»
# Mot alle kilder blir dekningen ~93 %, drevet av store cellelinje-
# eksperimenter (de største bidrar 5 191 proteiner hver). Det er ikke
# tallet manuskriptet påstår.
#
# Datakilder (lastet ned 8. august 2026, Vesiclepedia versjon 5.1):
#   VESICLEPEDIA_PROTEIN_MRNA_DETAILS_5.1.txt.gz
#   vesiclepedia_human_plasma_experiments.txt — 326 eksperiment-IDer med
#     SPECIES = "Homo sapiens" og SAMPLE = "Plasma", hentet fra
#     VESICLEPEDIA_EXPERIMENT_DETAILS_5.1.txt via nettleser (Cloudflare
#     blokkerer curl). Vesikkeltyper: 162 extracellular vesicles,
#     85 exosomes, 36 microparticles, 18 small EV, 17 microvesicles, øvrige få.
#
# R5 minor: «Please report the same proportion among the non-significant
# proteins; if comparable, the annotation is uninformative.»
#
# Krever: R1_01_longterm.R
# =====================================================================
suppressPackageStartupMessages({library(data.table)})
source(file.path(Sys.getenv("ECV_ROOT", "/Users/sindrle/Research/Prosjekter/ECV"),
                 "review_round_1/01_scripts/R1_config.R"))
ROOT <- ECV_ROOT
OUT  <- file.path(ROOT, "review_round_1/02_output"); setwd(ROOT)

VP  <- fread(cmd = paste("gzcat", shQuote(file.path(OUT, "VESICLEPEDIA_PROTEIN_MRNA_DETAILS_5.1.txt.gz"))))
setnames(VP, make.names(names(VP)))
plasma_ids <- as.integer(strsplit(readLines(
  file.path(OUT, "vesiclepedia_human_plasma_experiments.txt")), ",")[[1]])
cat(sprintf("Vesiclepedia 5.1: %d rader | humane plasma-eksperimenter: %d\n",
            nrow(VP), length(plasma_ids)))

# NB: CONTENT TYPE finnes som "protein", "Protein" og "protein " (mellomrom).
# Et likhetsfilter paa "protein" dropper 26 % av radene.
hum <- VP[SPECIES == "Homo sapiens" & tolower(trimws(CONTENT.TYPE)) == "protein"]
any_src <- unique(toupper(hum$GENE.SYMBOL))
plasma  <- unique(toupper(hum[EXPERIMENT.ID %in% plasma_ids]$GENE.SYMBOL))
cat(sprintf("unike gensymboler: alle kilder %d | kun plasma %d\n",
            length(any_src), length(plasma)))

PRIM <- fread(file.path(OUT, "R1_01_longterm_primary.csv"))
PRIM <- PRIM[!is.na(Genes) & Genes != ""]
setorder(PRIM, q); PRIM <- unique(PRIM, by = "Genes")
PRIM[, sym := toupper(Genes)]

tst <- function(lab, universe, hits) {
  sig <- universe[q <  0.05]; ns <- universe[q >= 0.05]
  a <- sum(sig$sym %in% hits); b <- nrow(sig) - a
  c_ <- sum(ns$sym  %in% hits); d <- nrow(ns)  - c_
  ft <- fisher.test(matrix(c(a, b, c_, d), 2))
  cat(sprintf("\n%s\n  signifikante  %4d/%4d = %5.1f%%\n  ikke-signif.  %4d/%4d = %5.1f%%\n  OR %.2f (95%% KI %.2f-%.2f), p = %.3g\n",
              lab, a, a+b, 100*a/(a+b), c_, c_+d, 100*c_/(c_+d),
              ft$estimate, ft$conf.int[1], ft$conf.int[2], ft$p.value))
  data.table(annotation = lab, sig_hit = a, sig_n = a+b, sig_pct = round(100*a/(a+b),1),
             ns_hit = c_, ns_n = c_+d, ns_pct = round(100*c_/(c_+d),1),
             OR = round(unname(ft$estimate),3), p = ft$p.value)
}

cat("\n================ PRIMÆRANALYSEN (reanalyse, 127 respondere) ================")
r1 <- tst("Vesiclepedia, ALLE kilder", PRIM, any_src)
r2 <- tst("Vesiclepedia, KUN humant plasma", PRIM, plasma)

res <- rbind(r1, r2)
fwrite(res, file.path(OUT, "R1_40_vesiclepedia_plasma.csv"))
cat("\nskrevet: R1_40_vesiclepedia_plasma.csv\n")
