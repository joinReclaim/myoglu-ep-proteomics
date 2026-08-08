# =====================================================================
# R1_19_fig3A_volcano.R — manuskriptets FIGUR 3, PANEL A (Figures/Figure1.tif)
# Erstatter 002_ECV_vulcano.R:237-284
#
# Uendret fra originalen: farger (red/black), punktstørrelse 0.3, theme_minimal-
# varianten, aksetitler, stiplet vertikal null-linje, etikettregel (rank_q +
# rank_eff, tre opp og tre ned), 4 x 3 TOMMER via cairo_pdf.
# Endret: nye data (127 sig, 76 opp / 51 ned av 1545 testet), x-akse strammet
# til det nye effektområdet.
# N=-annotasjonene lå ikke i originalkoden — de settes i PowerPoint.
# Krever: R1_01_longterm.R
# =====================================================================
suppressPackageStartupMessages({library(data.table); library(ggrepel)})
source(file.path(Sys.getenv("ECV_ROOT", "/Users/sindrle/Research/Prosjekter/ECV"),
                 "review_round_1/01_scripts/R1_config.R"))
ROOT <- ECV_ROOT
OUT  <- file.path(ROOT,"review_round_1/02_output"); setwd(ROOT)
source(file.path(ROOT,"review_round_1/01_scripts/R1_theme.R"))

A <- readRDS(file.path(OUT,"R1_01_longterm.rds"))$PRIM
A[, q_plot := pmax(q, .Machine$double.xmin)][, neglog10q := -log10(q_plot)]
A[, rank_combined := rank(q_plot) + rank(-abs(Estimate))]
ups <- head(A[Estimate > 0][order(rank_combined)]$Genes, 3)
dns <- head(A[Estimate < 0][order(rank_combined)]$Genes, 3)
cat(sprintf("sig %d (opp %d / ned %d) av %d testet\n",
            sum(A$q < .05), sum(A$q < .05 & A$Estimate > 0),
            sum(A$q < .05 & A$Estimate < 0), nrow(A)))
cat("etiketter:", paste(c(ups, dns), collapse = ", "), "\n")

p <- ggplot(A, aes(Estimate, neglog10q)) +
  geom_point(data = A[q <= 0.05], colour = PAL$volcano["sig"], size = 0.3) +
  geom_point(data = A[q >  0.05], colour = PAL$volcano["ns"],  size = 0.3) +
  scale_x_continuous(expression(paste(Delta, "log"[2], "(int)")),
                     # aksen går til +/-1.3 så alle punkter så vidt får plass
                     # (dataområde -0.75 til 1.18), men merkes kun -1 .. 1
                     expand = xp, limits = c(-1.3, 1.3),
                     breaks = seq(-1, 1, 0.5)) +
  scale_y_continuous(expression("-log"[10] * "(Q)"), expand = xp) +
  geom_vline(xintercept = 0, linetype = "dotted", colour = "black", linewidth = 0.5) +
  coord_cartesian(clip = "off") +
  geom_text_repel(data = A[Genes %in% c(ups, dns)], aes(label = Genes),
                  size = 3, nudge_y = 0.15, segment.size = 0.2, na.rm = TRUE,
                  min.segment.length = 0, seed = 1, family = fnt) +
  theme_ecv_min()
save_ecv("R1_19_fig3A_volcano.pdf", p, width = 4, height = 3, unit = "in")
