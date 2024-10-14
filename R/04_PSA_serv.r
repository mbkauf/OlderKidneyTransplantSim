### Script to run PSA
## Install and load Packages
library(doParallel)
library(doRNG)
library(foreach)
library(dtplyr)
library(dplyr)
library(data.table)
library(doSNOW)

## Load model inputs and functions
setwd("//tsclient/Documents/R/OlderKidneyTransplantSim")
source("R/01_model_inputs.R", echo = FALSE)
source("R/02_model_functions.R", echo = FALSE)

##### File paths #####
in_path  <- "//tsclient/Documents/Simulation Model/Model Outputs/Recalibration/"
samp_path  <- "//tsclient/Documents/Simulation Model/Model Inputs/"
out_path <- "//tsclient/Documents/Simulation Model/Results/"

samp_path <- "Z:/PSA/inputs/"
out_path <- "Z:/PSA/results/"

## Generate parameter sets
n_rep <- 10 # Number of repititions of each parameter set
n_cq  <- 50 # number of cost and qaly parameter sets
n     <- 64 * n_rep * n_cq

# Equations
rep_coef <- function(num_rep, m_coef) {
  m_rep_coef <- do.call(rbind, replicate(num_rep, m_coef, simplify = FALSE))
  # m_rep_coef <- m_rep_coef[rep(seq_len(nrow(m_rep_coef)), each = num_total), ]
  return(m_rep_coef)
}

ddtx_coef_trim       <- rep_coef(n_rep * n_cq, ddtx_coef_trim)
ldtx_coef_trim       <- rep_coef(n_rep * n_cq, ldtx_coef_trim)
mort_coef_trim       <- rep_coef(n_rep * n_cq, mort_coef_trim)
remove_coef_trim     <- rep_coef(n_rep * n_cq, remove_coef_trim)
mlogit_coef_trim     <- rep_coef(n_rep * n_cq, mlogit_coef_trim)
gl_coef_trim         <- rep_coef(n_rep * n_cq, gl_coef_trim)
gl_mort_coef_trim    <- rep_coef(n_rep * n_cq, gl_mort_coef_trim)
gs_mort_coef_trim    <- rep_coef(n_rep * n_cq, gs_mort_coef_trim)
dial_mort_coef_trim  <- rep_coef(n_rep * n_cq, dial_mort_coef_trim)

ddtx_shape_trim      <- rep(ddtx_shape_trim, n_rep * n_cq)
ldtx_shape_trim      <- rep(ldtx_shape_trim, n_rep * n_cq)
mort_shape_trim      <- rep(mort_shape_trim, n_rep * n_cq)
remove_shape_trim    <- rep(remove_shape_trim, n_rep * n_cq)
gl_shape_trim        <- rep(gl_shape_trim, n_rep * n_cq)
gs_mort_shape_trim   <- rep(gs_mort_shape_trim, n_rep * n_cq)
dial_mort_shape_trim <- rep(dial_mort_shape_trim, n_rep * n_cq)

##### Policy Inputs #####
v_tx_bump        <- seq(from = 0.05, to = 0.25, by = 0.05)
v_kdpi_pct_shift <- seq(from = 0.05, to = 0.20, by = 0.05)

comb <- function(...) {
  mapply("rbind", ..., SIMPLIFY = FALSE)
}

my_cluster <- parallel::makeCluster(
  80,
  type = "PSOCK"
)

l_psa_params <- get_psa_draws(l_qaly_params = l_qaly_params,
                              l_cost_params = l_cost_params,
                              n_draws = n_cq)

m_psa_params <- do.call(cbind, l_psa_params)
m_psa_params <- m_psa_params[rep(seq_len(nrow(m_psa_params)),
                                 each = n_rep * 64), ]

write.csv(m_psa_params, file = paste0(samp_path, "costs_qalys.csv"),
          row.names = FALSE)

doSNOW::registerDoSNOW(my_cluster)

iters <- n
pb <- txtProgressBar(max = iters, style = 3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress = progress)

out <- foreach(i = 1:n, .combine = "comb", .multicombine = TRUE,
               .options.snow = opts,
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

## Prepare data
z <- 0
j <- 0
x <- 0
v_id <- rep(seq.int(1:64), n_rep)
for (i in 1:20) {
  tmp_cq_counts <- matrix(nrow = n, ncol = 90)
  tmp_p         <- matrix(nrow = n, ncol = 2)
  z <- 0
  for (j in 1:n) {
    z <- z + 1
    pos <- x + j
    tmp_cq_counts[z, ]     <- out[[2]][[pos]][[1]]
    tmp_p[z, 1]     <- out[[2]][[pos]][[2]]
    tmp_p[z, 2]     <- out[[2]][[pos]][[3]]
  }

  tmp_cq_counts <- as.data.frame(tmp_cq_counts) %>%
    mutate(id = v_id) %>%
    group_by(id) %>%
    summarise(across(everything(), ~ sum(.x, na.rm = TRUE))) %>%
    as.matrix()

  assign(paste0("policy_", i, "_cq_counts"), tmp_cq_counts)
  assign(paste0("policy_", i, "_p"), tmp_p)
  x <- x + z
}

base_df <- as.data.frame(out[[1]]) %>%
  mutate(id = v_id) %>%
  group_by(id) %>%
  summarise(across(everything(), ~ sum(.x, na.rm = TRUE))) %>%
  as.matrix()

l_policy_df <- list(policy_1_cq_counts, policy_2_cq_counts, policy_3_cq_counts,
                    policy_4_cq_counts, policy_5_cq_counts, policy_6_cq_counts,
                    policy_7_cq_counts, policy_8_cq_counts, policy_9_cq_counts,
                    policy_10_cq_counts, policy_11_cq_counts,
                    policy_12_cq_counts, policy_13_cq_counts,
                    policy_14_cq_counts, policy_15_cq_counts,
                    policy_16_cq_counts, policy_17_cq_counts,
                    policy_18_cq_counts, policy_19_cq_counts,
                    policy_20_cq_counts)
names(l_policy_df) <-  c("policy_1_cq_counts",  "policy_2_cq_counts",
                         "policy_3_cq_counts",  "policy_4_cq_counts",
                         "policy_5_cq_counts",  "policy_6_cq_counts",
                         "policy_7_cq_counts",  "policy_8_cq_counts",
                         "policy_9_cq_counts",  "policy_10_cq_counts",
                         "policy_11_cq_counts", "policy_12_cq_counts",
                         "policy_13_cq_counts", "policy_14_cq_counts",
                         "policy_15_cq_counts", "policy_16_cq_counts",
                         "policy_17_cq_counts", "policy_18_cq_counts",
                         "policy_19_cq_counts", "policy_20_cq_counts")


out_path <- paste0("results/PSA/")

dir.create(out_path, showWarnings = FALSE)

lapply(names(l_policy_df), function(dname) {
  write.csv(l_policy_df[[dname]],
            file = paste0(out_path, dname, ".csv"), row.names = FALSE)
})

write.csv(base_df, file = paste0(out_path, "base_df.csv"), row.names = FALSE)

print(task_id)


