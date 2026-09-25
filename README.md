<!-- badges: start -->
<!-- For more info: https://usethis.r-lib.org/reference/badges.html -->
[![Lifecycle: experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
<!-- badges: end -->

# Spectacled Eider Population Viability Analysis (PVA)

## Overview

This repository contains the code, data references, and reporting for a Bayesian population viability analysis (PVA) of spectacled eiders (*Somateria fischeri*), built from an integrated population model (IPM). The analysis estimates the probability of quasi-extinction (a breeding population below 250 individuals, per the Spectacled Eider Recovery Plan; USFWS 1996) for the Yukon-Kuskokwim Delta (YKD) and Arctic Coastal Plain (ACP) populations under alternative future scenarios of environmental lead contamination and winter sea-ice conditions. Results support the U.S. Fish and Wildlife Service Species Status Assessment (SSA) for the species.

The IPM combines aerial survey counts, a multistate mark–recapture (survival) model based on banded birds from Kigigak Island, nest-success data, and expert-elicited priors into a single joint model, fit in JAGS. Population projections run to 2100.

**Two versions of the analysis exist:**

- **2021 version** — the original PVA, presented as Appendix B of the 2021 SSA (USFWS 2021). Implemented in R Markdown with an accompanying JAGS IPM.
- **2026 update** — a revision for the new SSA. The most significant change is a revised mark–recapture component that accommodates missing banding and resighting years (2016–2018 and 2020) within the m-array formulation; this required a dedicated simulation-validation effort. Other revisions include updated observed data (counts, nest success, and sea ice through the mid-2020s), revised priors to better match the expert elicitation, a shared ice effect across age classes, and a reparameterized detection-deviation term. See `NEWS.md` for a version-by-version change log and the report (below) for full methods.

> **Note on repository contents.** The file list under [Usage](#usage) describes the files relevant to the 2026 analysis (and their 2021 predecessors where applicable). The working repository may contain additional exploratory or superseded files not yet described here; these are pending cleanup.

## Installation

No installation of this repository is required. To reproduce the analysis you will need:

- **R** (developed under R 4.6.1) and **JAGS** (v4.3.1), called from R via the **jagsUI** package (v1.6.3).
- R packages used across the scripts and report: `jagsUI`, `tidyverse`, `knitr`, `cowplot`, `kableExtra`, `AKaerial`, `readxl`, `zoo`, `mgcv`, and `viridis` (via `scale_*_viridis_d`). 
- The `AKaerial` package (USFWS) supplies the aerial survey count data programmatically; see [https://github.com/USFWS/AKaerial](https://github.com/USFWS/AKaerial). 

## Usage

To reproduce the full analysis: (1) fit the YKD scenarios with `IPM_YKD_26.R`, (2) fit the ACP scenarios with `IPM_ACP_26.R`, then (3) render the report `SPEI_PVA_IPM_2026.qmd`, which reads the saved model fits and produces all figures and tables. Model fitting is computationally intensive (40,000 iterations × 4 chains per scenario, run in parallel).

The files below are grouped by role. Files marked *(referenced, not yet described)* are called by the scripts or report but were not describe in this draft of the README.

### Reports

| File | Description |
|------|-------------|
| `SPEI_PVA_IPM_2026.qmd` | Quarto source for the 2026 report. Contains the full methods, runs the post-processing/figures/tables from saved model fits, and renders to a (watermarked, line-numbered draft) PDF. This is the primary human-readable document for the 2026 update. |
| `SPEI_PVA_IPM_20201120.Rmd` | R Markdown source for the original 2021 report (Appendix B of the 2021 SSA). Retained for reference and for comparison of methods/priors between versions. |
| `references.bib` | BibTeX reference database used by the .qmd to produce the reference list in the 2026 report. |

### Models and analysis scripts

| File | Description |
|------|-------------|
| `IPM_YKD_26.R` | 2026 YKD analysis driver. Defines the JAGS model text and writes four scenario-specific model files (see *Generated JAGS files* below); prepares all input data (counts via `AKaerial`, the m-array from encounter histories, sea-ice covariates and their negative-binomial dispersion, nest-success covariates, and the breeding-propensity prior); and runs the YKD scenarios, saving fits to `MS_Scenarios/`. Includes the `marray()` helper (builds the multistate m-array from encounter histories) and a `plot_results()` helper. |
| `IPM_ACP_26.R` | 2026 ACP analysis driver. Same structure as the YKD script but for the ACP: a single JAGS model file (no lead submodel), ice-only scenarios, ACP-specific count handling (external detection correction, 1988–2006 entered as missing), ACP-specific density-dependence and initial-population priors, and shared Kigigak survival/fecundity. Saves fits to `MS_Scenarios/`. |
| `IPM_YKD.R` | Original (2021) YKD IPM script. Predecessor to `IPM_YKD_26.R`. |
| `IPM_MS.R` | Intermediate manuscript version derived from `IPM_YKD.R` (winter 2021–22 revisions), itself the basis for the 2026 scripts. |
| `QE_Prob.R` | Computes quasi-extinction probabilities by year from the saved posterior fits; `source()`d by the report to build the QE tables and cumulative-probability figures. |

### Generated JAGS model files

These are written by the R driver scripts (via `cat()`); they are outputs of the scripts rather than hand-edited sources, but are listed here because the report and fitting steps depend on them.

| File | Description |
|------|-------------|
| `YKD_IPM_L2008.jags` | YKD IPM, lead-decline-from-2008 scenario. |
| `YKD_IPM_L2026.jags` | YKD IPM, lead-decline-from-2026 scenario. |
| `YKD_IPM_Lconst.jags` | YKD IPM, constant-lead scenario. |
| `YKD_IPM_current.jags` | YKD IPM, "current conditions" scenario (lead declines 2008→2026, constant thereafter). |
| `ACP_IPM.jags` | ACP IPM (no lead submodel; ice scenario set by input data). |

### Prior-development scripts

These scripts derive informed priors used in the IPM; they are run separately (not part of the main fitting pipeline) and their outputs are hard-coded or otherwise carried into the IPM scripts.

| File | Description |
|------|-------------|
| `DDpriorsEO.R` | Derives the density-dependence priors from the expert-elicited carrying capacity (using the deterministic Lefkovitch projection). Source of the `Gamma` priors on the density-dependence coefficients for both populations. |
| `ObsPriors.R` | Derives the observation-model detection-deviation prior (the `Beta` latent variable for `d_t`) from expert-elicitation question 9. |
| `BP_prior_EO_revised.R` | Derives the joint multivariate-normal prior (mean vector and covariance) for the breeding-propensity coefficients (intercept, linear and quadratic ice effects). Source of `muBP`/`SigmaBP`. |
| `nb_regress.R` | Fits the negative-binomial relationship used to set the dispersion/size parameter for the sea-ice projections. |

### Mark–recapture (CJS) component and simulation validation

| File / folder | Description |
|---------------|-------------|
| `CMR_component/` | Standalone multistate mark–recapture (CJS) model and outputs, fit outside the IPM. Used to derive informed survival/detection/maturation priors and to provide the "independent" survival estimates compared against the IPM in the report. |
| `CMR_component/CMR.missing.RDS` | Saved fit of the standalone CJS model (with missing-year handling) read by the report to plot independent survival estimates. |
| `CMR_component/CMR_simulation/` | Simulation study validating the missing-year m-array implementation, including the cross-check against an independently coded individual-encounter-history model. Basis for the limitation noted in the report (weak local identifiability of detection following very-low-detection years). |

### Input data (`input_data/`)

Aerial survey counts for the 2026 analysis are pulled programmatically from the `AKaerial` package rather than from CSVs; several static inputs remain as files. Data is not maintained on GitHub or currently on a public repository. Contact the Maintainer for data. 

| File | Description |
|------|-------------|
| `eh_na_1992-2025_adults_ducklings.csv` | Individual encounter histories (Kigigak Island, 1992–2025), for adults and ducklings; processed into the multistate m-array. |
| `fecundity.csv` | Annual nest-success and clutch-size estimates (Kigigak Island) used to build the nest-success covariate and inform the fecundity prior. |
| `sea_ice_vars_1979_2025.xlsx` | Observed winter sea-ice metrics (1979–2025); source of the updated observed extreme-ice-day series. |
| `extreme.sea.ice.csv` | Sea-ice series including RCP4.5/RCP8.5 model-averaged projections, merged with the observed record for covariate construction and projection. |
| `YKD_counts_2026.csv` | YKD count series written out from the `AKaerial` pull in `IPM_YKD_26.R`, so the report can read it at render time without re-querying the package. Same source as the model inputs (not an independent series). |
| `ACP_counts_2026.csv` | ACP count series written out from the `AKaerial` pull in `IPM_ACP_26.R`, as above. |
| `lead.exposure.csv` | Lead-exposure scenario curves used for the lead-scenario figure. |
| `YKD_SPEI.csv`, `ACP_SPEI.csv` | 2021-era count inputs. Superseded by `AKaerial` pulls in 2026; retained for reference. |
| `extreme.sea.ice.csv` (2021 use), `YKD.QEresults.csv`, `ACP.QEresults.csv`, `YKD.NBresults.csv`, `ACP.NBresults.csv` | 2021-era precomputed results tables read by the original `.Rmd`. Superseded in 2026 by on-the-fly computation from saved fits. |

### Model output (`MS_Scenarios/`)

Saved JAGS fits (`.rds`), one per scenario, produced by the driver scripts and read by the report.

| File | Description |
|------|-------------|
| `YKD.L2008.4.5.rds`, `YKD.L2026.4.5.rds`, `YKD.Lconst.4.5.rds` | YKD fits, RCP4.5, under the three lead scenarios. |
| `YKD.L2008.8.5.rds`, `YKD.L2026.8.5.rds`, `YKD.Lconst.8.5.rds` | YKD fits, RCP8.5, under the three lead scenarios. |
| `YKD.current.rds` | YKD "current conditions" fit. |
| `ACP.4.5.rds`, `ACP.8.5.rds`, `ACP.current.rds` | ACP fits under the three ice scenarios. |
| `jags.data8.5.rds`, `jags.data4.5.rds`, `jags.data.current.rds` | Saved JAGS data bundles for the ACP scenarios (written by `IPM_ACP_26.R`). |

### Figures (`figures/`)

Figures rendered by the report are written here in `.svg`, `.png`, and `.pdf` (e.g., `counts`, `counts2`, `nest`, `lead`, `ice`, `traj`, `extinct`, `medYKD`, `medACP`, `sur`, `beta`, `qeykd`, `qeacp`). The life-history diagram (`figures/Picture1.png`) is a static input, not generated by the report. 

## Getting help

Contact the [project maintainer](mailto:erik_osnas@fws.gov) for help with this repository. 

## Contribute

Contact the project maintainer for information about contributing to this repository. Submit a [GitHub Issue](https://github.com/USFWS/SpectacledEiderPVA/issues) to report a bug or request a feature or enhancement. 

-----

![](https://i.creativecommons.org/l/zero/1.0/88x31.png) This work is licensed under a [Creative Commons Zero Universal v1.0 License](https://creativecommons.org/publicdomain/zero/1.0/).
