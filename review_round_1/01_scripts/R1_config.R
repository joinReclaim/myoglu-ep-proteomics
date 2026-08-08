# =====================================================================
# R1_config.R — all file paths in one place
#
# Sourced first by every R1_* script. Override with environment
# variables if the datasets live elsewhere:
#   export ECV_ROOT=/path/to/repository        # this repository
#   export ECV_PHENO=/path/to/MyoGlu_pheno.csv
#   export ECV_OLINK=/path/to/olink
#   export ECV_RNASEQ=/path/to/expression.RData
#
# EXTERNAL DATASETS — see README for availability status:
#   PHENO    MyoGlu clinical phenotypes (clamp, VO2max, MRI/MRS)
#   OLINK    Olink Explore 3072 NPX, serum, same participants
#   RNASEQ   RNA-seq from m. vastus lateralis and subcutaneous adipose tissue
#
# NOTE: this file is the only script that differs from the version used to
# produce the published results. Two path definitions pointing to a second,
# unpublished cohort have been removed; they were not used by any script in
# this repository. All other scripts are byte-identical to those that were run.
# =====================================================================
ge <- function(v, default) { x <- Sys.getenv(v); if (nzchar(x)) x else default }

ECV_ROOT <- ge("ECV_ROOT",  getwd())
OUT      <- file.path(ECV_ROOT, "review_round_1/02_output")
FIGDIR   <- file.path(ECV_ROOT, "review_round_1/03_figures")

PHENO_FILE <- ge("ECV_PHENO",
  file.path(ECV_ROOT, "external/MyoGlu_pheno_2019-23-08.csv"))
OLINK_DIR  <- ge("ECV_OLINK",  file.path(ECV_ROOT, "external/Olink"))
RNASEQ_RDA <- ge("ECV_RNASEQ",
  file.path(ECV_ROOT, "external/00_Expression_correlations.RData"))
