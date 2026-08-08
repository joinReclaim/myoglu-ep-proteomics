# =====================================================================
# R1_theme.R — figurstil hentet DIREKTE fra originalskriptene
#
# ALL FIGURTEKST SKAL VÆRE PÅ ENGELSK. Manuskriptet er engelsk.
#
# Palettene varierer per figurtype i originalkoden — bruk riktig én:
#   PAL$volcano   002_ECV_vulcano.R:238-239
#   PAL$heat      005:231, 006:216 (pheatmap)
#   PAL$group_xy  005:265, 436, 461 (XY-plott med gruppefarge)
#   PAL$group_box 101:32-33, 003:26-27 (boksplott)
#   PAL$venn      101:389, 201:33
#   PAL$time      202:622-625 (DIABLO)
#   PAL$cluster   utledet fra Figure7.tif (Tableau) — BEKREFT
# =====================================================================
suppressPackageStartupMessages({library(ggplot2); library(cowplot)})
source(file.path(Sys.getenv("ECV_ROOT", "/Users/sindrle/Research/Prosjekter/ECV"),
                 "review_round_1/01_scripts/R1_config.R"))

fnt <- "Helvetica"
pd  <- position_dodge(0.5)
xp  <- c(0.05, 0.05)
cols <- c("#C04140", "#4C5688")   # deklarert i originalskriptene

PAL <- list(
  volcano   = c(sig = "red", ns = "black"),
  heat      = colorRampPalette(c("blue", "white", "red"))(100),
  group_xy  = c(dysglyc = "red", control = "black"),
  group_box = c("midnightblue", "darkred"),
  venn      = c("#4C78A8", "#F58518"),
  # NB: etikettene MÅ være alfabetisk sortert i kronologisk rekkefølge.
  # cimDiablo sorterer legendetiketter alfabetisk men tegner fargeruter i
  # faktorrekkefølge — med "Before/After" blir legenden feil. "Baseline"
  # før "Post-training" (B < P) gjør at begge stemmer.
  time      = c("Baseline" = "#1F77B4", "Post-training" = "#FF7F0E"),
  cluster   = c(C1 = "#1F77B4", C2 = "#FF7F0E", C3 = "#2CA02C", C4 = "#9467BD")
)
FIGDIR <- "/Users/sindrle/Research/Prosjekter/ECV/review_round_1/03_figures"

# theme_minimal-varianten (volcano, 002:251-265)
theme_ecv_min <- function(base = 12) {
  theme_minimal() +
    theme(text = element_text(size = base, family = fnt),
          panel.grid = element_blank(), panel.border = element_blank(),
          legend.position = "none", legend.title = element_blank(),
          strip.background = element_blank(), axis.ticks = element_line(),
          axis.title = element_text(size = base),
          axis.text = element_text(colour = "black", size = base, family = fnt),
          axis.line = element_line(colour = "black", linewidth = 0.75))
}
# theme_classic-varianten (005, 101 m.fl.)
theme_ecv <- function(base = 12) {
  theme_classic() +
    theme(text = element_text(size = base, family = fnt),
          panel.grid = element_blank(), panel.border = element_blank(),
          legend.title = element_blank(), strip.background = element_blank(),
          plot.title = element_text(hjust = 0.5), axis.ticks = element_line(),
          axis.title = element_text(size = base),
          axis.text = element_text(colour = "black", size = base, family = fnt),
          axis.line = element_line(colour = "black", linewidth = 0.75))
}
# 101_ECV_acute_responses.R og 003_ECV_FHL1.R la på gjennomsiktig bakgrunn.
# 005-avledede figurer gjorde IKKE det. Legg denne på kun der originalen hadde den.
theme_transparent <- function() {
  theme(plot.background  = element_rect(fill = "transparent", colour = NA),
        panel.background = element_rect(fill = "transparent", colour = NA),
        legend.background = element_rect(fill = "transparent", colour = NA),
        legend.box.background = element_rect(fill = "transparent", colour = NA))
}
# unit = "in" eller "cm" — originalene bruker BEGGE, se per figur
save_ecv <- function(file, plot = ggplot2::last_plot(), width, height, unit = "cm") {
  dir.create(FIGDIR, showWarnings = FALSE, recursive = TRUE)
  cowplot::ggsave2(file.path(FIGDIR, file), plot = plot, width = width,
                   height = height, unit = unit, device = cairo_pdf, family = fnt)
  message("figur: ", file.path(FIGDIR, file))
}
