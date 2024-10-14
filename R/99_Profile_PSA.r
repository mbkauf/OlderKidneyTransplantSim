# Profile PSA
### Script to run PSA
## Install and load Packages
if (!require(doParallel)) install.packages("doParallel")
if (!require(doRNG)) install.packages("doRNG")
if (!require(foreach)) install.packages("foreach")
if (!require(dtplyr)) install.packages("dtplyr")
if (!require(dplyr)) install.packages("dplyr")
if (!require(data.table)) install.packages("data.table")

library(doParallel)
library(doRNG)
library(foreach)
library(dtplyr)
library(dplyr)
library(data.table)

## Load model inputs and functions
source("R/01_model_inputs.R", echo = FALSE)
source("R/02_model_functions.R", echo = FALSE)

comb <- function(...) {
  mapply("rbind", ..., SIMPLIFY = FALSE)
}


my_cluster <- parallel::makeCluster(
  4,
  type = "PSOCK"
)

doParallel::registerDoParallel(cl = my_cluster)
# doRNG::registerDoRNG(seed = 12345)

n <- 4

out <- foreach(i = 1:n, .combine = "comb", .multicombine = TRUE,
               .packages = c("MASS", "copula", "psych",
                             "dplyr", "survival", "dtplyr")) %dopar% {
  ## Run Baseline
  # Waitlist Simulation
  list_sim_df <- list.simulation(
    seed        = round(runif(1, min = 1, max = 99999)),
    data        = sample_df,
    ddtx_coef   = ddtx_coef_trim[i, ],
    ddtx_p      = ddtx_shape_trim[i],
    ldtx_coef   = ldtx_coef_trim[i, ],
    ldtx_gamma  = ldtx_shape_trim[i],
    mort_coef   = mort_coef_trim[i, ],
    mort_shape  = mort_shape_trim[i],
    remove_coef = remove_coef_trim[i, ],
    remove_p    = remove_shape_trim[i]
  )

  # Post-transplant Simulation
  sim_df <- run.simulation(
    data            = list_sim_df,
    seed            = round(runif(1, min = 1, max = 99999)),
    gl30_coef       = mlogit_coef_trim[i, 1:18],
    dgf_coef        = mlogit_coef_trim[i, 19:36],
    mort30_coef     = mlogit_coef_trim[i, 37:54],
    gl_coef         = gl_coef_trim[i, ],
    gl_shape        = gl_shape_trim[i],
    gs_mort_coef    = gs_mort_coef_trim[i, ],
    gs_mort_shape   = gs_mort_shape_trim[i],
    dial_mort_coef  = dial_mort_coef_trim[i, ],
    dial_mort_shape = dial_mort_shape_trim[i],
    gl_mort_coef    = gl_mort_coef_trim[i, ]
  )

  tmp_baseline <- sim_df

  # Make list of costs and QALYs for this run

  m_base_costs_qalys <- f_costs_qalys_psa(data = tmp_baseline, r = disc,
                                          params = l_psa_params, x = 1)

  ## Run Policies
  l_policy_outcomes <- list()
  z <- 0
  for (x in v_tx_bump) {
    for (y in v_kdpi_pct_shift) {
      z <- z + 1
      # Waitlist Simulation
      list_sim_df <- list.simulation(
        seed        = round(runif(1, min = 1, max = 99999)),
        data        = sample_df,
        ddtx_coef   = ddtx_coef_trim[i, ],
        ddtx_p      = ddtx_shape_trim[i],
        ldtx_coef   = ldtx_coef_trim[i, ],
        ldtx_gamma  = ldtx_shape_trim[i],
        mort_coef   = mort_coef_trim[i, ],
        mort_shape  = mort_shape_trim[i],
        remove_coef = remove_coef_trim[i, ],
        remove_p    = remove_shape_trim[i],
        ddtx_rate = x
      )

      # Post-transplant Simulation
      sim_df <- run.simulation(
        data            = list_sim_df,
        seed            = round(runif(1, min = 1, max = 99999)),
        gl30_coef       = mlogit_coef_trim[i, 1:18],
        dgf_coef        = mlogit_coef_trim[i, 19:36],
        mort30_coef     = mlogit_coef_trim[i, 37:54],
        gl_coef         = gl_coef_trim[i, ],
        gl_shape        = gl_shape_trim[i],
        gs_mort_coef    = gs_mort_coef_trim[i, ],
        gs_mort_shape   = gs_mort_shape_trim[i],
        dial_mort_coef  = dial_mort_coef_trim[i, ],
        dial_mort_shape = dial_mort_shape_trim[i],
        gl_mort_coef    = gl_mort_coef_trim[i, ],
        france = TRUE,
        kdpi_pct_shift = y
      )

      tmp_policy <- sim_df

      m_policy_costs_qalys <- f_costs_qalys_psa(data = tmp_policy, r = disc,
                                                params = l_psa_params, x = 1)
      l_policy_outcomes[[z]] <- list(m_policy_costs_qalys,
                                     x,
                                     y)
    }
  }

  list(m_base_costs_qalys, l_policy_outcomes)

}
parallel::stopCluster(cl = my_cluster)
