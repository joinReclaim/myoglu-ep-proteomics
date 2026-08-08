# =====================================================================
# R1_22_fig3D_module_pathways.R — manuskriptets FIGUR 3, PANEL D
# Erstatter 006_ECV_pathway_analysis.R:342-425 («Red cluster»)
#
# Uendret: samme dotplot-layout som panel B (facet_grid source ~ analysis,
# størrelse = Count, form = ontologi, x = -log10(FDR)), theme_classic.
# Endret: modulene kommer fra USENTRERT k=2-klynging i R1_21 (ikke k=7,
#   ikke sentrert — sentreringens begrunnelse falt bort da det globale
#   skiftet viste seg å være en lasteffekt som loess allerede fjerner), og
# universe = MÅLT proteom.
# Krever: R1_01_longterm.R, R1_21_fig3C_heatmap.R
# =====================================================================
suppressPackageStartupMessages({library(data.table); library(clusterProfiler); library(org.Hs.eg.db)})
source(file.path(Sys.getenv("ECV_ROOT", "/Users/sindrle/Research/Prosjekter/ECV"),
                 "review_round_1/01_scripts/R1_config.R"))
ROOT <- ECV_ROOT
OUT  <- file.path(ROOT,"review_round_1/02_output"); setwd(ROOT)
source(file.path(ROOT,"review_round_1/01_scripts/R1_theme.R"))
P  <- readRDS(file.path(OUT,"R1_01_longterm.rds"))$PRIM
md <- fread(file.path(OUT,"R1_21_fig3C_modules.csv"))

U <- unique(na.omit(suppressMessages(bitr(unique(P$Genes), "SYMBOL","ENTREZID", org.Hs.eg.db))$ENTREZID))
cat(sprintf("universe (målt proteom): %d | moduler: %s\n", length(U),
            paste(sprintf("%s n=%d", names(table(md$module_uncentred)),
                          as.integer(table(md$module_uncentred))), collapse = ", ")))

run <- function(genes, ont) {
  G <- unique(na.omit(suppressMessages(bitr(genes,"SYMBOL","ENTREZID",org.Hs.eg.db))$ENTREZID))
  e <- tryCatch(if (ont == "KEGG")
        enrichKEGG(G, universe = U, organism = "hsa", pAdjustMethod = "BH", pvalueCutoff = .05)
      else enrichGO(G, OrgDb = org.Hs.eg.db, keyType = "ENTREZID", ont = ont,
        universe = U, pAdjustMethod = "BH", pvalueCutoff = .05, readable = TRUE),
      error = function(e) NULL)
  if (is.null(e) || !nrow(as.data.frame(e))) return(NULL)
  as.data.table(as.data.frame(e))[, ontology := ont]
}
res <- rbindlist(lapply(sort(unique(md$module_uncentred)), function(m)
  rbindlist(lapply(c("BP","CC","KEGG"), function(o) {
    d <- run(md[module_uncentred == m, Genes], o)
    if (!is.null(d)) d[, module := m]; d }), fill = TRUE)), fill = TRUE)

if (!nrow(res)) { cat("ingen signifikante termer\n"); quit(save="no") }
res[, source := factor(fcase(ontology=="BP","GO BP", ontology=="CC","GO CC",
                             ontology=="KEGG","KEGG"), levels=c("GO BP","GO CC","KEGG"))]
res[, analysis := factor(module, levels = sort(unique(md$module_uncentred)))]
res[, neglog10q := -log10(p.adjust)]
top <- res[order(p.adjust), head(.SD, 5), by = .(analysis, source)]
cat("\ntermer per modul:\n"); print(top[, .N, by = .(analysis, source)])
cat("\ntopp termer:\n")
print(top[order(analysis, p.adjust)][, .(analysis, source, Description,
        GeneRatio, BgRatio, q = signif(p.adjust,3), Count)][1:min(12,.N)])

p <- ggplot(top, aes(x = neglog10q, y = reorder(Description, neglog10q))) +
  geom_point(aes(size = Count, shape = source), alpha = 0.85) +
  facet_grid(source ~ analysis, scales = "free_y", space = "free_y") +
  # headroom så de største punktene ikke kuttes av panelkanten
  scale_x_continuous(breaks = seq(0, 15, 5), limits = c(0, 19),
                     expand = expansion(mult = 0.06)) +
  theme_classic() +
  theme(text = element_text(size = 12, family = fnt),
        plot.title = element_text(hjust = 0.5),
        axis.title = element_text(size = 12),
        axis.text = element_text(colour = "black", size = 12, family = fnt),
        axis.text.x = element_text(colour = "black", size = 12, family = fnt),
        axis.line = element_line(colour = "black", linewidth = 0.75)) +
  labs(x = expression("-log"[10]*"(FDR)"), y = NULL, size = "Proteins", shape = NULL)
save_ecv("R1_22_fig3D_module_pathways.pdf", p, width = 26, height = 13, unit = "cm")
fwrite(top, file.path(OUT,"R1_22_fig3D_terms.csv"))
