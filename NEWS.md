# NEWS

All notable changes to the Spectacled Eider PVA/IPM are documented here, most recent first.

## 2026 update (for the 2026 SSA)

The 2026 revision updates the observed data through the mid-2020s and makes several structural and prior changes to the IPM. The single most significant change is a revised mark–recapture component that accommodates missing banding and resighting years within the m-array formulation.

### Mark–recapture (survival) component

- **Missing-year handling (major change).** The multistate mark–recapture component was revised to accommodate years with no banding and/or resighting effort (2016, 2017, 2018, and 2020). Detection is multiplied by a binary resighting-effort indicator so that observation is suspended in no-effort years while survival continues to operate across them; the affected probability mass routes automatically into the "not resighted" category of the m-array. There was no known precedent for doing this in the m-array format used here, so a dedicated simulation-validation effort was undertaken (`CMR_component/CMR_simulation/`).
- **Indexing bug fixed / code clarified.** A bug in the earlier CJS indexing was identified and corrected; the JAGS code was also revised to make the effort indexing more transparent, with a programmatic check on the index. 
- **Extended encounter series.** The encounter-history series was extended to 1992–2025 (`input_data/eh_na_1992-2025_adults_ducklings.csv`).
- **Shared ice effect on survival.** Duckling and post-hatch age classes now share a single pair of ice coefficients (linear and quadratic), differing only in their intercepts, to induce correlation in survival across age classes from their shared overwintering distribution. The 2021 model estimated separate juvenile and adult ice coefficients.
- **Simplifications from the simulation work.** Maturation probability ($\alpha$) and mean detection ($\mu_p$) were shown to be identifiable from the data and are now given informed priors from the standalone CJS fit rather than expert elicitation; $\alpha$ is treated as time-constant.

### Priors

- **Survival intercepts** updated from the standalone CJS refit on the extended series: $\mu_{\phi_0}$ from Beta(15, 45) → Beta(20, 60); $\mu_{\phi_a}$ from Beta(90, 10) → Beta(84, 16).
- **Maturation** $\alpha$: from a truncated inverse-gamma (expert-elicited) → Beta(13, 21) (from CJS model).
- **Mean detection** $\mu_p$: from Beta(1, 1) → Beta(47, 53) (from CJS model).
- **Breeding-propensity coefficients**: three independent univariate priors → a single joint multivariate-normal prior (`BP_prior_EO_revised.R`) that preserves the elicited correlations among the intercept and the linear and quadratic ice effects.
- **Density dependence**: reparameterized from a normal prior on $\log(\beta_{DD}^{-1})$ (carrying-capacity scale) to a Gamma prior placed directly on the coefficient (`DDpriorsEO.R`); the density term now enters on a per-1000-female scale. YKD uses Gamma(3.1, 251.648); ACP uses Gamma(3.02, 158.30).
- **Observation detection deviation** ($d_t$): from a truncated normal with an estimated variance → a fixed distribution derived from a latent Beta variable mapped to a detection range (`ObsPriors.R`, elicitation question 9). YKD uses Beta(4, 4) on a 0.2–0.8 range; ACP uses Beta(2.8, 2.8) on a 0.5–1.0 range centered on 0.75.

### Data updates

- **Counts** now pulled programmatically from the `AKaerial` package (YKD 1988–2026; ACP 2007–2026) rather than from static CSVs. The ACP count timeline is aligned to the mark–recapture series by entering 1988–2006 as missing.
- **YKD missing count years** are now 2011, 2015, and 2020 (2020 added; the 2015 low estimate is retained as missing, with an updated rationale relating the low count to nesting activity/nest success rather than solely a novice observer).
- **ACP missing count years** include the 1988–2006 pre-survey period plus 2020 and 2021 (COVID-related survey gap).
- **Nest success** extended with estimates from Friendly et al. (2025) for 2019 and 2021; remaining recent years without data use an informed placeholder prior. 
- **Sea ice** updated through the 2024–25 winter (`sea_ice_vars_1979_2025.xlsx`); the negative-binomial dispersion parameter for ice projections is now re-estimated per scenario (RCP4.5 ≈ 3.46, RCP8.5 ≈ 3.40, current ≈ 3.33) via an `mgcv::gam` fit, replacing the fixed 2021 values (3.75 / 3.5). 

### ACP model

- **Survival priors no longer ACP-specific.** The 2021 ACP model used ACP-specific survival intercepts (Beta(2.3, 4), Beta(26, 3)) reflecting an elicited belief that ACP survival exceeds YKD survival. The 2026 ACP shares the Kigigak-derived survival intercepts with the YKD, on the grounds that (1) it is a simplification, (2) the mark–recapture likelihood overwhelms the prior difference, and (3) the shared Kigigak baseline represents low—but non-zero—lead exposure, appropriate for the ACP given evidence of low-level ACP lead exposure. The YKD/ACP survival difference is now expressed through the YKD's lead term rather than through separate priors.
- **Observer levels**: ACP now uses eight observer-team levels; YKD uses six.

### Reporting and outputs

- Report reimplemented in **Quarto** (`SPEI_PVA_IPM_2026.qmd`), replacing the 2021 R Markdown. QE tables and figures are now computed on the fly from saved fits (`QE_Prob.R`) rather than from precomputed CSVs.
- **Quasi-extinction horizons** reported as cumulative probability of QE before 2066 (40 years) and before 2100 (74 years); posterior sample size for QE calculation increased to 80,000. Median (rather than mean) trajectories are emphasized in the projection figures.
- Software versions: JAGS 4.3.1, jagsUI 1.6.3, R 4.6.1. MCMC settings 40,000 iterations × 4 chains. 

## 2021 version (Appendix B, 2021 SSA)

- Initial Bayesian PVA derived from an integrated population model for the YKD and ACP spectacled eider populations (USFWS 2021).
- Combined aerial survey counts, a multistate mark–recapture survival model (Kigigak Island), Kigigak nest-success/clutch-size data, and expert-elicited priors (IDEA protocol) in a joint model fit in JAGS.
- YKD scenarios spanned three lead-exposure futures × two climate (RCP4.5/RCP8.5) scenarios plus a current-conditions model; ACP scenarios varied climate only. Projections to 2100; quasi-extinction defined as a breeding population below 250 individuals (USFWS 1996).
- Implemented in R Markdown (`SPEI_PVA_IPM_20201120.Rmd`) with static CSV inputs and precomputed results tables.
