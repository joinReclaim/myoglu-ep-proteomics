# =====================================================================
# R1_20_fig3B_pathways.R — manuskriptets FIGUR 3, PANEL B (Figures/Figure1.tif)
# Erstatter 006_ECV_pathway_analysis.R:301-340
#
# Uendret fra originalen: facet_grid(source ~ analysis, scales/space free_y),
# punktstørrelse = Count, form = ontologi, x = -log10(FDR), theme_classic,
# fem termer per blokk, 19 x 10 cm.
# Endret: universe = MÅLT proteom (1522) i stedet for helgenom (18986), og
# separate kolonner for alle / opp / ned (R2 concern 4).
# Krever: R1_03_pathways.R
# =====================================================================
suppressPackageStartupMessages({library(data.table)})
source(file.path(Sys.getenv("ECV_ROOT", "/Users/sindrle/Research/Prosjekter/ECV"),
                 "review_round_1/01_scripts/R1_config.R"))
ROOT <- ECV_ROOT
OUT  <- file.path(ROOT,"review_round_1/02_output"); setwd(ROOT)
source(file.path(ROOT,"review_round_1/01_scripts/R1_theme.R"))

d <- fread(file.path(OUT,"R1_03_pathways_measured_universe.csv"))
d[, analysis := factor(fcase(set == "alle", "All responders",
                             set == "opp",  "Increased",
                             set == "ned",  "Decreased"),
                       levels = c("All responders","Increased","Decreased"))]
d[, source := factor(fcase(ontology == "BP",   "GO BP",
                           ontology == "CC",   "GO CC",
                           ontology == "KEGG", "KEGG"),
                     levels = c("GO BP","GO CC","KEGG"))]
d[, neglog10q := -log10(p.adjust)]
top <- d[order(p.adjust), head(.SD, 5), by = .(analysis, source)]
cat("termer per panel:\n"); print(dcast(top[, .N, by=.(source, analysis)],
                                        source ~ analysis, value.var = "N", fill = 0))

p <- ggplot(top, aes(x = neglog10q, y = reorder(Description, neglog10q))) +
  geom_point(aes(size = Count, shape = source), alpha = 0.85) +
  facet_grid(source ~ analysis, scales = "free_y", space = "free_y") +
  scale_x_continuous(breaks = c(0, 3, 6, 9), limits = c(0, 10.5)) +
  theme_classic() +
  theme(text = element_text(size = 12, family = fnt),
        plot.title = element_text(hjust = 0.5),
        axis.title = element_text(size = 12),
        axis.text = element_text(colour = "black", size = 12, family = fnt),
        axis.text.x = element_text(colour = "black", size = 12, family = fnt),
        axis.line = element_line(colour = "black", linewidth = 0.75)) +
  labs(x = expression("-log"[10]*"(FDR)"), y = NULL, size = "Proteins", shape = NULL)
save_ecv("R1_20_fig3B_pathways.pdf", p, width = 26, height = 11, unit = "cm")
fwrite(top[order(analysis, source, p.adjust),
           .(analysis, source, Description, GeneRatio, BgRatio,
             FoldEnrichment = round(FoldEnrichment,1), q = signif(p.adjust,3), Count)],
       file.path(OUT,"R1_20_fig3B_terms.csv"))
