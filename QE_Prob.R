quasi_extinction <- function(posterior, threshold, last_observed_col = 39) {
  # --------------------------------------------------------------------------
  # Calculate quasi-extinction probability from a JAGS posterior projection.
  #
  # ARGUMENTS
  #   posterior         : matrix of posterior population sizes.
  #                       - rows = MCMC samples (draws from the posterior)
  #                       - cols = time points (years), chronological order
  #   threshold         : numeric, the quasi-extinction threshold (X). A sample
  #                       is quasi-extinct the first year it drops below X.
  #   last_observed_col : integer, the column index of the LAST observed
  #                       (data-informed) year. All columns up to and including
  #                       this one are dropped before calculation, so that the
  #                       first retained column is the first FORECAST year.
  #                       Default 0 = no trimming (matrix is already all forecast
  #                       years).
  #
  # WHY TRIMMING MATTERS
  #   Your posterior object contains past years where the population was
  #   observed via counts (data-informed years) followed by projected years.
  #   This function computes the probability of *future* quasi-extinction, so
  #   column 1 of the analyzed matrix must be the first year AFTER the last
  #   observed year. If observed years are left in, any sample already below the
  #   threshold during the observed period would be miscounted as "extinct at
  #   year 1," mixing past dynamics into a forward-looking risk estimate.
  #   (Assumes observed counts were above the threshold; if not, handle that
  #   case separately as a modeling decision.)
  #
  #   Example: if the last observed year is column 20, pass
  #            last_observed_col = 20, and forecasting starts at column 21.
  #
  # RETURNS a list with:
  #   first_below       : per-sample column index (within the FORECAST window)
  #                       of the first year below X; NA if never. Map to calendar
  #                       years if needed, e.g. forecast_years[first_below].
  #   prob_ever         : P(quasi-extinction ever occurs in the horizon)
  #   prob_by_year      : cumulative P(quasi-extinct by year t), per forecast col
  #   median_first_year : median first-crossing year among samples that cross
  #
  # NOTE: uses '<'. Switch to '<=' if hitting exactly X counts as
  #       quasi-extinct.
  # --------------------------------------------------------------------------
  
  # Trim observed (data-informed) columns so col 1 = first forecast year
  if (last_observed_col > 0) {
    if (last_observed_col >= ncol(posterior)) {
      stop("last_observed_col leaves no forecast years to analyze.")
    }
    posterior <- posterior[, (last_observed_col + 1):ncol(posterior), drop = FALSE]
  }
  
  n_years <- ncol(posterior)
  
  # For each MCMC sample (row), find the first forecast year below threshold
  first_below <- apply(posterior, 1, function(x) {
    hit <- which(x < threshold)
    if (length(hit) == 0) NA_integer_ else hit[1]
  })
  
  # Cumulative quasi-extinction probability through each forecast year
  prob_by_year <- sapply(seq_len(n_years), function(t) {
    mean(!is.na(first_below) & first_below <= t)
  })
  
  list(
    first_below       = first_below,
    prob_ever         = mean(!is.na(first_below)),
    prob_by_year      = prob_by_year,
    median_first_year = median(first_below, na.rm = TRUE)
  )
}
##------------------------------------------------------------------------------
# load data and call function
# Example:
# posteriors <- list(
#   scenario_A = sims.list$N_A,
#   scenario_B = sims.list$N_B
# )
# 
library(jagsUI)
# results <- lapply(posteriors, quasi_extinction, threshold = 50)
# 
# # e.g. probability of ever hitting quasi-extinction, per scenario
# sapply(results, `[[`, "prob_ever")
# A single test run:
# YKD.L2008.4.5 <- readRDS("MS_Scenarios/YKD.L2008.4.5.rds")$sims.list$Nb
# qe <- quasi_extinction(YKD.L2008.4.5, 250)
# qe
# plot(qe$prob_by_year)

#On the real posteriors:
posteriors <- list(
  YKD.current = readRDS("MS_Scenarios/YKD.current.rds")$sims.list$Nb,
  YKD.L2008.4.5 = readRDS("MS_Scenarios/YKD.L2008.4.5.rds")$sims.list$Nb,
  YKD.L2026.4.5 = readRDS("MS_Scenarios/YKD.L2026.4.5.rds")$sims.list$Nb,
  YKD.constant.4.5 = readRDS("MS_Scenarios/YKD.Lconst.4.5.rds")$sims.list$Nb,
  YKD.L2008.8.5 = readRDS("MS_Scenarios/YKD.L2008.8.5.rds")$sims.list$Nb,
  YKD.L2026.8.5 = readRDS("MS_Scenarios/YKD.L2026.8.5.rds")$sims.list$Nb,
  YKD.constant.8.5 = readRDS("MS_Scenarios/YKD.Lconst.8.5.rds")$sims.list$Nb,
  ACP.current = readRDS("MS_Scenarios/ACP.current.rds")$sims.list$Nb,
  ACP.4.5 = readRDS("MS_Scenarios/ACP.4.5.rds")$sims.list$Nb,
  ACP.8.5 = readRDS("MS_Scenarios/ACP.8.5.rds")$sims.list$Nb
)

results <- lapply(posteriors, quasi_extinction, threshold = 250)
