# =====================================================================
# R1_37_go_acute_chronic.R — GO/KEGG: akutt vs kronisk respons
# Adresserer R2 punkt 4.2
#
# R1_03_pathways.R kjørte ORA for den KRONISKE responsen mot riktig
# bakgrunn (det målte proteomet, ikke helgenomet) og separat for opp- og
# nedregulerte. Her gjentas NØYAKTIG samme oppsett for den AKUTTE
# responsen (R1_05_acute.csv, q_time < 0.05), og de to stilles opp mot
# hverandre.
#
# UNIVERSPROBLEMET (håndteres eksplisitt):
#   Kronisk (R1_01, A1 vs B1) og akutt (R1_05, B1->B2->B3) tester IKKE
#   samme proteinsett — 70 %-dekningskravet slår ulikt ut, og lmer-
#   modellene faller fra på litt ulike proteiner. Derfor gjøres to
#   analyser:
#     (a) "own"    — hver respons mot SITT EGET testede sett. Kronisk-
#                    delen reproduserer da R1_03 eksakt (verifiseres
#                    under), men de to søylene er ikke direkte
#                    sammenlignbare fordi nevneren er ulik.
#     (b) "common" — SNITTET av de to testede settene brukes som universe
#                    for BEGGE, og begge genlister trimmes til snittet.
#                    Dette er den gyldige hode-mot-hode-sammenligningen
#                    og den som brukes i figur og hovedtabell.
#
# Krever: R1_01_longterm.R, R1_03_pathways.R, R1_05_acute_chronic.R
# =====================================================================
suppressPackageStartupMessages({
  library(data.table); library(clusterProfiler); library(org.Hs.eg.db)
  library(ggplot2)
})
source(file.path(Sys.getenv("ECV_ROOT", "/Users/sindrle/Research/Prosjekter/ECV"),
                 "review_round_1/01_scripts/R1_config.R"))
ROOT <- ECV_ROOT
OUT  <- file.path(ROOT, "review_round_1/02_output"); setwd(ROOT)
S <- readRDS(file.path(OUT, "R1_00_normalised.rds"))
source(file.path(ROOT, "review_round_1/01_scripts/R1_theme.R"))

set.seed(1)
cat(sprintf("design: %d injeksjoner, %d deltakere (%s)\n",
            ncol(S$Xloe), uniqueN(S$si$ID2),
            paste(sprintf("%s=%d", names(table(S$si$SampleCode)),
                          as.integer(table(S$si$SampleCode))), collapse = " ")))

# --- inn: de to responsene -------------------------------------------
CHR <- readRDS(file.path(OUT, "R1_01_longterm.rds"))$PRIM       # A1 vs B1
AC  <- fread(file.path(OUT, "R1_05_acute.csv"))                 # B1->B2->B3

univ_chr <- unique(CHR$Genes)
univ_ac  <- unique(AC$Genes)
univ_cmn <- intersect(univ_chr, univ_ac)
cat(sprintf("\n=== UNIVERSER (symboler) ===\nkronisk testet %d | akutt testet %d | snitt %d | kun kronisk %d | kun akutt %d\n",
            length(univ_chr), length(univ_ac), length(univ_cmn),
            length(setdiff(univ_chr, univ_ac)), length(setdiff(univ_ac, univ_chr))))

# Signaturer. Kronisk: fortegn på Estimate (som i R1_03).
# Akutt: Time-leddet er en omnibus-test med 2 frihetsgrader, så "retning"
# defineres av fortegnet rett etter økt (est_B2 = B2-B1), som i R1_05.
chr_sig <- CHR[q < 0.05]
ac_sig  <- AC[q_time < 0.05]
cat(sprintf("\nkronisk FDR<0.05: %d (opp %d / ned %d)\nakutt   FDR<0.05: %d (opp %d / ned %d)  [retning = est_B2]\n",
            uniqueN(chr_sig$Genes), uniqueN(chr_sig[Estimate > 0]$Genes),
            uniqueN(chr_sig[Estimate < 0]$Genes), uniqueN(ac_sig$Genes),
            uniqueN(ac_sig[est_B2 > 0]$Genes), uniqueN(ac_sig[est_B2 < 0]$Genes)))
cat(sprintf("akutt: fortegn B2 og B3 er enige for %d av %d proteiner\n",
            sum(sign(ac_sig$est_B2) == sign(ac_sig$est_B3)), nrow(ac_sig)))
ov <- intersect(unique(ac_sig$Genes), unique(chr_sig$Genes))
cat(sprintf("proteinoverlapp akutt/kronisk: %d (%s)\n", length(ov),
            paste(sort(ov), collapse = ", ")))

# --- mapping til Entrez ----------------------------------------------
to_entrez <- function(x) unique(na.omit(suppressMessages(
  bitr(x, "SYMBOL", "ENTREZID", org.Hs.eg.db))$ENTREZID))
U <- list(chr_own = to_entrez(univ_chr),
          ac_own  = to_entrez(univ_ac),
          common  = to_entrez(univ_cmn))
cat(sprintf("\nEntrez-mappede universer: kronisk %d | akutt %d | felles %d\n",
            length(U$chr_own), length(U$ac_own), length(U$common)))

# --- ORA, samme kall som R1_03 ---------------------------------------
# pvalueCutoff/qvalueCutoff = 1 slik at HELE termtabellen kommer ut; BH
# beregnes uansett over alle testede termer, så p.adjust er identisk med
# et kall med cutoff 0.05. Signifikans defineres etterpå med
# clusterProfilers standard (p.adjust < 0.05 & qvalue < 0.2), som gir
# nøyaktig samme termer som R1_03 — verifiseres nedenfor.
run_ora <- function(genes, Uni, resp, lab, utag) {
  G <- to_entrez(genes); if (!length(G)) return(NULL)
  rbindlist(lapply(c("BP", "CC", "KEGG"), function(ont) {
    e <- tryCatch({
      if (ont == "KEGG")
        enrichKEGG(G, universe = Uni, organism = "hsa", pAdjustMethod = "BH",
                   pvalueCutoff = 1, qvalueCutoff = 1)
      else
        enrichGO(G, OrgDb = org.Hs.eg.db, keyType = "ENTREZID", ont = ont,
                 universe = Uni, pAdjustMethod = "BH", pvalueCutoff = 1,
                 qvalueCutoff = 1, readable = TRUE)
    }, error = function(e) NULL)
    if (is.null(e)) return(NULL)
    d <- as.data.table(as.data.frame(e))
    if (!nrow(d)) return(NULL)
    d[, `:=`(response = resp, set = lab, ontology = ont, universe = utag)][]
  }), fill = TRUE)
}

# genlister per universe-variant (trimmet til det aktuelle universet)
sets_for <- function(keep) list(
  list(resp = "chronic", set = "all",  g = intersect(unique(chr_sig$Genes), keep)),
  list(resp = "chronic", set = "up",   g = intersect(unique(chr_sig[Estimate > 0]$Genes), keep)),
  list(resp = "chronic", set = "down", g = intersect(unique(chr_sig[Estimate < 0]$Genes), keep)),
  list(resp = "acute",   set = "all",  g = intersect(unique(ac_sig$Genes), keep)),
  list(resp = "acute",   set = "up",   g = intersect(unique(ac_sig[est_B2 > 0]$Genes), keep)),
  list(resp = "acute",   set = "down", g = intersect(unique(ac_sig[est_B2 < 0]$Genes), keep)))

jobs <- c(
  lapply(sets_for(univ_chr)[1:3], function(x) c(x, list(U = U$chr_own, utag = "own"))),
  lapply(sets_for(univ_ac)[4:6],  function(x) c(x, list(U = U$ac_own,  utag = "own"))),
  lapply(sets_for(univ_cmn),      function(x) c(x, list(U = U$common,  utag = "common"))))
cat("\ngenlister som sendes til ORA (etter trimming til universet):\n")
for (j in jobs) cat(sprintf("  %-8s %-5s %-6s n=%d\n", j$resp, j$set, j$utag, length(j$g)))

FULL <- rbindlist(lapply(jobs, function(j)
  run_ora(j$g, j$U, j$resp, j$set, j$utag)), fill = TRUE)
FULL[, sig := p.adjust < 0.05 & qvalue < 0.2]
cat(sprintf("\ntermer testet totalt: %d | signifikante: %d\n", nrow(FULL), sum(FULL$sig)))
if (!"KEGG" %in% FULL$ontology)
  cat("MERK: ingen KEGG-resultater — sjekk nettverkstilgang til rest.kegg.jp\n")

# --- verifikasjon mot R1_03 ------------------------------------------
old <- fread(file.path(OUT, "R1_03_pathways_measured_universe.csv"))
old[, set := fcase(set == "alle", "all", set == "opp", "up", set == "ned", "down")]
new <- FULL[response == "chronic" & universe == "own" & sig == TRUE]
k <- function(d) paste(d$ID, d$set, d$ontology)
cat(sprintf("\n=== VERIFIKASJON mot R1_03 (kronisk, eget universe) ===\nR1_03 %d termer | her %d | felles %d | kun R1_03 %d | kun her %d\n",
            nrow(old), nrow(new), length(intersect(k(old), k(new))),
            length(setdiff(k(old), k(new))), length(setdiff(k(new), k(old)))))
if (length(setdiff(k(old), k(new))) || length(setdiff(k(new), k(old))))
  cat("  ADVARSEL: avvik fra R1_03 — sjekk clusterProfiler/GO-versjon\n")

# --- effekten av å bytte universe ------------------------------------
cat("\n=== ANTALL SIGNIFIKANTE TERMER: eget vs felles universe ===\n")
print(dcast(FULL[sig == TRUE, .N, by = .(response, set, universe, ontology)],
            response + set + ontology ~ universe, value.var = "N", fill = 0))

# --- hode mot hode på FELLES universe ---------------------------------
CM <- FULL[universe == "common"]
W <- dcast(CM, ID + Description + ontology + set ~ response,
           value.var = c("p.adjust", "GeneRatio", "Count", "FoldEnrichment"))
setnames(W, gsub("^p\\.adjust_", "q_", names(W)))
for (cc in c("q_acute", "q_chronic")) if (!cc %in% names(W)) W[, (cc) := NA_real_]
W[, sig_acute   := !is.na(q_acute)   & q_acute   < 0.05]
W[, sig_chronic := !is.na(q_chronic) & q_chronic < 0.05]
# clusterProfilers qvalue-filter tas med for å matche R1_03-definisjonen
sigkey <- FULL[sig == TRUE & universe == "common", paste(ID, set, response)]
W[, sig_acute   := paste(ID, set, "acute")   %in% sigkey]
W[, sig_chronic := paste(ID, set, "chronic") %in% sigkey]
W[, status := fcase(sig_acute & sig_chronic, "shared",
                    sig_acute,               "acute only",
                    sig_chronic,             "chronic only",
                    default = "neither")]

cat("\n=== FELLES / UNIKE TERMER (felles universe) ===\n")
print(dcast(W[status != "neither", .N, by = .(set, ontology, status)],
            set + ontology ~ status, value.var = "N", fill = 0))
for (s in c("all", "up", "down")) {
  d <- W[set == s & status != "neither"]
  if (!nrow(d)) { cat(sprintf("\n--- %s: ingen signifikante termer i noen av responsene ---\n", s)); next }
  cat(sprintf("\n--- %s (%d termer signifikante i minst én respons) ---\n", s, nrow(d)))
  print(head(d[order(pmin(q_acute, q_chronic, na.rm = TRUE))][
    , .(ontology, Description, status, qA = signif(q_acute, 2),
        qC = signif(q_chronic, 2), nA = Count_acute, nC = Count_chronic)], 25))
}

# --- er forskjellen systematisk? (terskelfritt) -----------------------
cat("\n=== SYSTEMATISK FORSKJELL? — terskelfri rangkorrelasjon ===\n")
cat("Spearman mellom -log10(p) for de SAMME termene i de to responsene.\n",
    "Høy rho = de to responsene handler om det samme; nær null = de er\n",
    "uavhengige programmer. Kun termer testet i begge inngår.\n", sep = "")
# ---------------------------------------------------------------------
# ⚠ RETNINGSSPLITTET KAN IKKE SAMMENLIGNES PÅ TVERS AV DE TO RESPONSENE.
# Kronisk retning = fortegn(B1-A1). Akutt retning = fortegn(B2-B1).
# B1 inngår i begge med MOTSATT fortegn, så de to retningsetikettene er
# antijustert ved konstruksjon, ikke ved biologi:
#   alle 12 proteiner signifikante i begge har motsatt fortegn ved B2
#   akutt-opp ∩ kronisk-opp = 0, akutt-ned ∩ kronisk-ned = 0
#   Spearman(kronisk estimat, est_B2) over hele universet = -0.74
#
# ORA innen HVER respons for seg er gyldig — settene er definert innenfor
# sin egen analyse. Det som IKKE kan tolkes er "up" mot "up" og "down" mot
# "down" PÅ TVERS. Radene for set = "up" og "down" under er derfor
# deskriptive for hver respons, og rho-tallene for dem skal ikke siteres
# som mål på om de to programmene deler retning.
#
# Spørsmålet "deler akutt og kronisk respons retning" besvares i stedet
# terskelfritt i R1_05_acute_chronic.R med split-half over disjunkte
# deltakere (rho = -0.41), som ikke deler noen måling.
# ---------------------------------------------------------------------
cat("\n  ⚠ set = 'up'/'down' kan IKKE sammenlignes paa tvers av responsene:\n",
    "    B1 inngaar i begge kontrastene med motsatt fortegn, saa retnings-\n",
    "    etikettene er antijustert ved konstruksjon. Se R1_05 for det\n",
    "    gyldige svaret (split-half, rho = -0.41).\n", sep = "")
rho <- rbindlist(lapply(c("all", "up", "down"), function(s)
  rbindlist(lapply(unique(CM$ontology), function(o) {
    d <- W[set == s & ontology == o & !is.na(q_acute) & !is.na(q_chronic)]
    dd <- merge(CM[response == "acute"   & set == s & ontology == o, .(ID, pA = pvalue)],
                CM[response == "chronic" & set == s & ontology == o, .(ID, pC = pvalue)], by = "ID")
    if (nrow(dd) < 20) return(NULL)
    ct <- suppressWarnings(cor.test(-log10(dd$pA), -log10(dd$pC), method = "spearman"))
    data.table(set = s, ontology = o, n_terms = nrow(dd),
               rho = round(unname(ct$estimate), 3), p = signif(ct$p.value, 3))
  }), fill = TRUE)), fill = TRUE)
if (nrow(rho)) print(rho) else cat("  for få felles testede termer\n")

# hvor godt gjør de kroniske termene det i akuttanalysen, og motsatt
cat("\nkryssoppslag — FDR for den ene responsens termer i den andre:\n")
for (s in c("all", "up", "down"))
  for (dir in list(c("chronic", "acute"), c("acute", "chronic"))) {
    a <- dir[1]; b <- dir[2]
    d <- W[set == s & get(paste0("sig_", a)) == TRUE]
    if (!nrow(d)) { cat(sprintf("  [%s] %s: ingen signifikante termer\n", s, a)); next }
    qb <- d[[paste0("q_", b)]]
    cat(sprintf("  [%-4s] %d termer sign. i %-7s -> i %-7s median FDR %.2g, min %.2g, %d/%d med FDR<0.05 (%d ikke testbar)\n",
                s, nrow(d), a, b, median(qb, na.rm = TRUE),
                suppressWarnings(min(qb, na.rm = TRUE)),
                sum(!is.na(qb) & qb < 0.05), nrow(d), sum(is.na(qb))))
  }

# --- hva bærer de akutte termene? --------------------------------------
# Med bare 35 signifikante proteiner er hver term båret av en håndfull
# gener; skriv dem ut så leseren kan se hvor tynt grunnlaget er.
acsig <- FULL[universe == "common" & response == "acute" & sig == TRUE]
if (nrow(acsig)) {
  cat("\n=== GENER BAK DE AKUTT-SIGNIFIKANTE TERMENE (felles universe) ===\n")
  for (i in seq_len(nrow(acsig))) {
    ids <- strsplit(acsig$geneID[i], "/")[[1]]
    sym <- if (acsig$ontology[i] == "KEGG")
      suppressMessages(bitr(ids, "ENTREZID", "SYMBOL", org.Hs.eg.db))$SYMBOL else ids
    cat(sprintf("  %-5s %-38s %s  FDR=%.3g  genes: %s\n", acsig$ontology[i],
                substr(acsig$Description[i], 1, 38), acsig$GeneRatio[i],
                acsig$p.adjust[i], paste(sym, collapse = ", ")))
  }
} else cat("\nIngen akutt-signifikante termer å bryte ned.\n")

# --- CSV ---------------------------------------------------------------
outtab <- rbind(
  W[status != "neither", .(universe = "common", set, ontology, ID, Description,
                           status, GeneRatio_acute, GeneRatio_chronic,
                           Count_acute, Count_chronic,
                           FE_acute = round(FoldEnrichment_acute, 2),
                           FE_chronic = round(FoldEnrichment_chronic, 2),
                           q_acute = signif(q_acute, 3), q_chronic = signif(q_chronic, 3))],
  FULL[universe == "own" & sig == TRUE,
       .(universe = "own", set, ontology, ID, Description,
         status = paste(response, "only universe"),
         GeneRatio_acute = ifelse(response == "acute", GeneRatio, NA_character_),
         GeneRatio_chronic = ifelse(response == "chronic", GeneRatio, NA_character_),
         Count_acute = ifelse(response == "acute", Count, NA_integer_),
         Count_chronic = ifelse(response == "chronic", Count, NA_integer_),
         FE_acute = ifelse(response == "acute", round(FoldEnrichment, 2), NA_real_),
         FE_chronic = ifelse(response == "chronic", round(FoldEnrichment, 2), NA_real_),
         q_acute = ifelse(response == "acute", signif(p.adjust, 3), NA_real_),
         q_chronic = ifelse(response == "chronic", signif(p.adjust, 3), NA_real_))],
  fill = TRUE)
setorder(outtab, universe, set, ontology, q_chronic, q_acute)
fwrite(outtab, file.path(OUT, "R1_37_go_acute_vs_chronic.csv"))
cat(sprintf("\nskrev %s (%d rader)\n",
            file.path(OUT, "R1_37_go_acute_vs_chronic.csv"), nrow(outtab)))

# --- FIGUR: de to responsene side om side ------------------------------
# Alle termer som er signifikante i MINST ÉN respons (felles universe),
# med FDR fra BEGGE responsene på samme rad. Maks 8 termer per retnings-
# sett for lesbarhet, men alle akutt-signifikante termer beholdes alltid.
fd <- rbindlist(lapply(c("all", "up", "down"), function(s) {
  d <- W[set == s & status != "neither"][order(pmin(q_acute, q_chronic, na.rm = TRUE))]
  if (!nrow(d)) return(NULL)
  rbind(head(d, 8), d[sig_acute == TRUE])[!duplicated(ID)]
}), fill = TRUE)
if (nrow(fd)) {
  L <- melt(fd[, .(ID, set, Description, ontology, status, q_acute, q_chronic)],
            id.vars = c("ID", "set", "Description", "ontology", "status"),
            variable.name = "response", value.name = "q")
  L[, response := factor(fifelse(response == "q_acute", "Acute", "Chronic"),
                         levels = c("Acute", "Chronic"))]
  L[, panel := factor(fcase(set == "all", "All responders", set == "up", "Increased",
                            default = "Decreased"),
                      levels = c("All responders", "Increased", "Decreased"))]
  L[, lab := sprintf("%s  [%s]", Description,
                     fcase(ontology == "BP", "BP", ontology == "CC", "CC", default = "KEGG"))]
  L[, tested := !is.na(q)]
  L[, neglog10q := -log10(q)]
  L[tested == FALSE, neglog10q := 0]               # termen er ikke testbar i den responsen
  ordr <- L[, .(o = max(neglog10q)), by = .(panel, key = paste(set, ID))][order(panel, o)]
  L[, ykey := factor(paste(set, ID), levels = ordr$key)]
  # NB: navngitt vektor — med facet + scales="free_y" matcher ggplot en
  # UNAVN'T labels-vektor posisjonelt innenfor hvert panel, som gir feil
  # etiketter. Navngitt vektor matches på nivånavn og er trygg.
  ul    <- L[!duplicated(ykey)]
  ylabs <- setNames(ul$lab, as.character(ul$ykey))
  cols2 <- unname(PAL$venn)

  p <- ggplot(L, aes(x = neglog10q, y = ykey, colour = response)) +
    geom_line(aes(group = ykey), colour = "grey70", linewidth = 0.4) +
    geom_vline(xintercept = -log10(0.05), linetype = "dashed",
               colour = "grey40", linewidth = 0.4) +
    geom_point(aes(shape = tested), size = 2.6) +
    facet_grid(panel ~ ., scales = "free_y", space = "free_y") +
    scale_y_discrete(labels = ylabs) +
    scale_colour_manual(values = c(Acute = cols2[2], Chronic = cols2[1])) +
    scale_shape_manual(values = c(`TRUE` = 16, `FALSE` = 1), guide = "none") +
    theme_ecv(10) +
    theme(legend.position = "top", strip.text.y = element_text(angle = 0)) +
    labs(x = expression("-log"[10]*"(FDR)"), y = NULL,
         caption = paste0("Over-representation of the acute (B1→B2→B3) and chronic (pre→post) responses\n",
                          "against the same background: the ", length(U$common),
                          " proteins tested in both analyses. Dashed line: FDR = 0.05.\n",
                          "Open symbols: term not testable in that response (no genes from that list)."))
  save_ecv("R1_37_go_acute_vs_chronic.pdf", p, width = 24, height = 18, unit = "cm")
} else {
  cat("\nINGEN termer signifikante i noen av responsene på felles universe — figur droppes.\n")
}

cat("\n=== KONKLUSJON (tall å sitere) ===\n")
cat(sprintf("akutt: %d/%d proteiner FDR<0.05 mot kronisk %d/%d\n",
            uniqueN(ac_sig$Genes), length(univ_ac), uniqueN(chr_sig$Genes), length(univ_chr)))
cat(sprintf("felles universe %d proteiner; signifikante termer: akutt %d, kronisk %d, felles %d\n",
            length(U$common), sum(W$sig_acute), sum(W$sig_chronic),
            sum(W$status == "shared")))
