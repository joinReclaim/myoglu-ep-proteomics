# Analysis code — plasma extracellular particle-enriched proteome after exercise training

Code and figure source data accompanying:

> Lee-Ødegård S, *et al.* *[title]*. **Cardiovascular Diabetology** (submitted).
> ProteomeXchange accession **PXD______**.

This repository contains the full analysis for the revised manuscript: 22
numbered scripts, the numeric source data behind the figure panels, and the
vector figure files themselves.

---

## Quick start

```bash
git clone https://github.com/joinReclaim/myoglu-ep-proteomics.git && cd myoglu-ep-proteomics

export ECV_ROOT=$(pwd)
export ECV_PHENO=/path/to/MyoGlu_pheno_2019-23-08.csv
export ECV_OLINK=/path/to/olink
export ECV_RNASEQ=/path/to/00_Expression_correlations.RData

# place pg.report_Frode_plasma.csv and Frode_Norheim_sample_info.csv
# (from ProteomeXchange) in the repository root

Rscript run_all.R
```

`run_all.R` reports which external datasets are missing and then runs everything
it can. Individual scripts can also be run on their own, provided
`R1_00_setup.R` has been run once — it builds the cache
(`review_round_1/02_output/R1_00_normalised.rds`) that all the others read.

**R packages:** `data.table`, `lme4`, `lmerTest`, `car`, `ggplot2`,
`clusterProfiler`, `org.Hs.eg.db`, `mixOmics`, `ComplexHeatmap`, `ggrepel`,
`eulerr`.

---

## Layout

```
run_all.R                          entry point, run order
review_round_1/01_scripts/         24 R scripts (22 analysis + config + theme)
review_round_1/02_output/          28 CSVs — numeric source data for the panels
review_round_1/03_figures/         21 PDFs — vector figure panels as produced by R
```

All paths are resolved in `01_scripts/R1_config.R` and can be overridden with
the environment variables above. Figure styling — colour palettes, fonts, panel
theme — is centralised in `01_scripts/R1_theme.R`.

---

## Data availability

| Input | Status |
|---|---|
| Mass spectrometry raw files, DIA-NN output, sample metadata | **Public** — ProteomeXchange **PXD______** |
| Clinical phenotypes (hyperinsulinaemic-euglycaemic clamp, VO₂max, strength, MRI/MRS) | MyoGlu trial, NCT01803568 — on reasonable request |
| Serum Olink Explore 3072 (NPX) | MyoGlu trial — on reasonable request |
| Skeletal muscle and adipose tissue RNA-seq | MyoGlu trial — on reasonable request |

The scripts are therefore not runnable end to end from this repository alone.
The proteomics data are public; the clinical and transcriptomic data are
governed by the consent under which participants enrolled and by the approval
from the Regional Committee for Medical and Health Research Ethics, and are
shared on request within those limits.

`review_round_1/02_output/` exists precisely so this does not block anyone: it
holds the numeric values behind the figure panels, so they can be reproduced or
redrawn without access to individual-level data. **No file in this repository
contains clinical phenotype data.** With one exception, noted below, every file
is aggregate — one row per protein, pathway or cluster.

### Two panels whose source data is not deposited

`R1_27_fig8A_scores.csv` is the one participant-level file: two latent-variable
scores per participant, timepoint and block. They are subject-centred
projections of 500 features onto two components, from which nothing can be
reconstructed, so the panel can be redrawn by anyone. The participant codes match
those in the ProteomeXchange deposit.

Two panels cannot be supported this way:

| Panel | Underlying values | Why not deposited |
|---|---|---|
| Fig. 5A, 5B | Participant-level protein and phenotype pairs | Clinical phenotypes |
| Fig. 8B | Selected features × samples | Individual-level serum Olink and tissue RNA-seq values |

For these the repository provides the summary statistics annotated on the panels
(`R1_24_fig5A_rho.csv`, `R1_25_fig5B_stats.csv`,
`R1_27_fig8_selected_features.csv`) and the code that draws them, but not the
participant-level input. This follows from the consent rather than from any
choice about openness: anyone granted access to the restricted datasets can
regenerate both panels by running the scripts unchanged.

---

## What the code does differently from the originally submitted analysis

The reanalysis prepared during revision differs from the first submission in
four respects, all implemented here:

1. **Normalisation.** Cyclic loess against a row-median reference, fitted on
   observed pairs only (NA-safe), applied across all 104 injections. The first
   submission used no between-sample normalisation.
2. **Missing values.** No imputation. Proteins are required to have ≥70% valid
   values within each timepoint being contrasted; the threshold is varied from
   50% to 100% as a sensitivity analysis in `R1_01_longterm.R`.
3. **Batch.** The nanoLC failure split the injection sequence cleanly in two.
   Batch enters as a **model term** in the mixed model rather than being removed
   from the matrix beforehand.
4. **Multiplicity.** Benjamini-Hochberg FDR across the entire protein ×
   phenotype matrix, and a permutation test for the composite score in which
   the *selection step is repeated inside each permutation*.

Effect estimates for the main contrast are stable across all normalisation ×
completeness combinations; the sensitivity table is printed by
`R1_01_longterm.R`.

---

## Not included

Seven scripts that draw on a second, as yet unpublished cohort are withheld
until that study is published. They support a cross-platform methodological
comparison discussed in the response to reviewers and underpin no result
reported in the manuscript; their absence does not affect reproducibility of
anything here.

`01_scripts/R1_config.R` is the only file that differs from the version that was
run: two path definitions pointing at that cohort were removed. No script in
this repository referenced them. Every other script is byte-identical to the
code that produced the reported numbers.

---

## A note on the comments

Several scripts carry long header comments documenting why a step is as it is —
the normalisation assumption, the derivation of the batch variable, the
identity of the sample-to-file mapping, and a set of analytical pitfalls
encountered during the reanalysis. These are deliberate: they are the fullest
record of the reasoning behind the pipeline. Some are written in Norwegian.

---

## Licence

MIT — see `LICENSE`. The licence covers the code. Data availability is governed
by the terms in the table above.

## Citation

Repository: https://github.com/joinReclaim/myoglu-ep-proteomics

Archived at Zenodo: **DOI to be inserted on release**.
