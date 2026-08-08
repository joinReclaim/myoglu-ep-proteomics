# =====================================================================
# run_all.R — entry point
#
# Run from the root of this repository:
#     Rscript run_all.R
#
# R1_00_setup.R must run first; it builds the cache
# (review_round_1/02_output/R1_00_normalised.rds) that every other
# script reads. The remaining scripts are otherwise independent of
# each other and can be run individually.
# =====================================================================

if (!nzchar(Sys.getenv("ECV_ROOT"))) Sys.setenv(ECV_ROOT = normalizePath("."))
ROOT <- Sys.getenv("ECV_ROOT")
SCR  <- file.path(ROOT, "review_round_1/01_scripts")

if (!file.exists(file.path(SCR, "R1_config.R")))
  stop("ECV_ROOT does not point at the repository root: ", ROOT)

source(file.path(SCR, "R1_config.R"))

# --- check external datasets before starting -------------------------
ext <- c(PHENO = PHENO_FILE, OLINK = OLINK_DIR, RNASEQ = RNASEQ_RDA)
missing <- ext[!file.exists(ext)]
if (length(missing)) {
  message("\nExternal datasets not found:\n",
          paste0("  ", names(missing), "  ", missing, collapse = "\n"),
          "\n\nSee README, section 'Data availability'. Scripts that need a\n",
          "missing dataset will fail; the rest will run.\n")
}
if (!file.exists(file.path(ROOT, "pg.report_Frode_plasma.csv")))
  message("Mass spectrometry matrix not found at ", ROOT,
          "/pg.report_Frode_plasma.csv\n",
          "Download it from ProteomeXchange (see README) and place it here.\n")

# --- run order -------------------------------------------------------
ORDER <- c(
  # preprocessing — must be first
  "R1_00_setup.R",
  # primary analyses
  "R1_01_longterm.R",            # 12-week responder list
  "R1_02_phenotypes.R",          # phenotype matrix, FDR, permutation
  "R1_03_pathways.R",            # ORA against the measured proteome
  "R1_04_group_time.R",          # group x time interaction
  "R1_05_acute_chronic.R",       # acute vs chronic, threshold-free
  "R1_11_serum_ep_overlap.R",    # EP fraction vs serum PEA          [needs OLINK]
  "R1_12_mrna_concordance.R",    # EP protein vs tissue mRNA         [needs RNASEQ]
  "R1_13_diablo.R",              # multi-block integration
  "R1_15_partial_phenotype.R",   # mutually adjusted associations
  "R1_16_acute_clusters.R",      # acute response patterns
  "R1_17_chronic_clusters.R",    # correlation clustering of responders
  # figure panels
  "R1_19_fig3A_volcano.R",
  "R1_20_fig3B_pathways.R",
  "R1_21_fig3C_heatmap.R",
  "R1_22_fig3D_module_pathways.R",
  "R1_23_fig4_heatmap.R",
  "R1_24_fig5AB.R",
  "R1_25_fig5CD.R",
  "R1_26_fig7_systematic.R",
  "R1_27_fig8_diablo_panels.R",
  "R1_29_fig9_panels.R")

for (s in ORDER) {
  cat("\n", strrep("=", 70), "\n", s, "\n", strrep("=", 70), "\n", sep = "")
  ok <- tryCatch({ source(file.path(SCR, s), echo = FALSE); TRUE },
                 error = function(e) { message("FAILED: ", conditionMessage(e)); FALSE })
  if (!ok && s == "R1_00_setup.R")
    stop("R1_00_setup.R failed; nothing downstream can run.")
}
cat("\nDone. Output in review_round_1/02_output, figures in review_round_1/03_figures\n")
