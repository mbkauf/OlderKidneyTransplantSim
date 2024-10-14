### Script to run PSA
## Install and load Packages
library(doParallel)
library(doRNG)
library(foreach)
library(dtplyr)
library(dplyr)
library(data.table)

## Load model inputs and functions
source("R/01_model_inputs.R", echo = FALSE)
source("R/02_model_functions.R", echo = FALSE)

## Set sherlock dependent variables
if (Sys.getenv("LMOD_SYSHOST") == "sherlock") {
  task_id <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))
  in_path <- "/home/users/mbkauf/OlderKidneyTransplantSim/data/Model Inputs/"
} else {
  # task_id <- 1
  in_path <- ""
}

## Generate parameter sets
n_rep <- 10 # Number of repititions of each parameter set

# Equations
ddtx_coef_trim       <- do.call(rbind, replicate(n_rep, ddtx_coef_trim,
                                                 simplify = FALSE))
ldtx_coef_trim       <- do.call(rbind, replicate(n_rep, ldtx_coef_trim,
                                                 simplify = FALSE))
mort_coef_trim       <- do.call(rbind, replicate(n_rep, mort_coef_trim,
                                                 simplify = FALSE))
remove_coef_trim     <- do.call(rbind, replicate(n_rep, remove_coef_trim,
                                                 simplify = FALSE))
mlogit_coef_trim     <- do.call(rbind, replicate(n_rep, mlogit_coef_trim,
                                                 simplify = FALSE))
gl_coef_trim         <- do.call(rbind, replicate(n_rep, gl_coef_trim,
                                                 simplify = FALSE))
gl_mort_coef_trim    <- do.call(rbind, replicate(n_rep, gl_mort_coef_trim,
                                                 simplify = FALSE))
gs_mort_coef_trim    <- do.call(rbind, replicate(n_rep, gs_mort_coef_trim,
                                                 simplify = FALSE))
dial_mort_coef_trim  <- do.call(rbind, replicate(n_rep, dial_mort_coef_trim,
                                                 simplify = FALSE))
ddtx_shape_trim      <- rep(ddtx_shape_trim, n_rep)
ldtx_shape_trim      <- rep(ldtx_shape_trim, n_rep)
mort_shape_trim      <- rep(mort_shape_trim, n_rep)
remove_shape_trim    <- rep(remove_shape_trim, n_rep)
gl_shape_trim        <- rep(gl_shape_trim, n_rep)
gs_mort_shape_trim   <- rep(gs_mort_shape_trim, n_rep)
dial_mort_shape_trim <- rep(dial_mort_shape_trim, n_rep)

##### Policy Inputs #####
v_tx_bump        <- seq(from = 0.05, to = 0.25, by = 0.05)
v_kdpi_pct_shift <- seq(from = 0.05, to = 0.20, by = 0.05)

# Costs and QALYs
# Save to model inputs folder
s <- 12345 + task_id
set.seed(s)

l_psa_params <- get_psa_draws(l_qaly_params = l_qaly_params,
                              l_cost_params = l_cost_params,
                              n_draws = 1000)

m_psa_params <- do.call(cbind, l_psa_params)
write.csv(m_psa_params,
          file = paste0(in_path, "PSA/", task_id, "/costs_qalys.csv"),
          row.names = FALSE)


# use the environment variable SLURM_NTASKS_PER_NODE to set
# the number of cores to use
comb <- function(...) {
  mapply("rbind", ..., SIMPLIFY = FALSE)
}


my_cluster <- parallel::makeCluster(
  (Sys.getenv("SLURM_NTASKS_PER_NODE")),
  type = "FORK"
)

doParallel::registerDoParallel(cl = my_cluster)

n     <- 64 * n_rep

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

out_path <- paste0("results/PSA/", task_id, "/")

lapply(names(l_policy_df), function(dname) {
  write.csv(l_policy_df[[dname]],
            file = paste0(out_path, dname, ".csv"), row.names = FALSE)
})
write.csv(base_df, file = paste0(out_path, "base_df.csv"), row.names = FALSE)
