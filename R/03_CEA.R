#-------------------------------------------------------------------------------
# Run CEA
# 03_CEA.R
# Date Created: 07/11/2023
# Last Updated:
# Matthew Kaufmann, Stanford University
# Note:


##### Load libraries, model inputs, and functions #####
library(doParallel)
library(doRNG)
library(foreach)
library(dtplyr)
library(dplyr)
library(doSNOW)
library(data.table)
library(survminer)
library(gridExtra)
library(stringr)
library(tidyr)

source("R/01_model_inputs.R", echo = FALSE)
source("R/02_model_functions.R", echo = FALSE)

list2env(l_cost_params, envir = .GlobalEnv)
list2env(l_qaly_params, envir = .GlobalEnv)

##### File paths #####
in_path  <- "//tsclient/Documents/Simulation Model/Model Outputs/Recalibration/"
samp_path  <- "//tsclient/Documents/Simulation Model/Model Inputs/"
out_path <- "//tsclient/Documents/Simulation Model/Results/"

##### Policy Inputs #####
v_tx_bump        <- seq(from = 0.05, to = 0.25, by = 0.05)
v_kdpi_pct_shift <- seq(from = 0.05, to = 0.20, by = 0.05)

##### Run Policy Simulation 100 Times #####
comb <- function(...) {
  mapply("rbind", ..., SIMPLIFY = FALSE)
}

parallel::detectCores()
n_cores <- 80

my_cluster <- parallel::makeCluster(
  n_cores,
  type = "PSOCK"
)

doParallel::registerDoParallel(cl = my_cluster)
doRNG::registerDoRNG(seed = 12345)

n_rep <- 10  # Number of repetitions of each parameter set
ddtx_coef_trim <- do.call(rbind, replicate(n_rep, ddtx_coef_trim, simplify=FALSE))
ldtx_coef_trim <- do.call(rbind, replicate(n_rep, ldtx_coef_trim, simplify=FALSE))
mort_coef_trim <- do.call(rbind, replicate(n_rep, mort_coef_trim, simplify=FALSE))
remove_coef_trim <- do.call(rbind, replicate(n_rep, remove_coef_trim, simplify=FALSE))
mlogit_coef_trim <- do.call(rbind, replicate(n_rep, mlogit_coef_trim, simplify=FALSE))
gl_coef_trim <- do.call(rbind, replicate(n_rep, gl_coef_trim, simplify=FALSE))
gl_mort_coef_trim <- do.call(rbind, replicate(n_rep, gl_mort_coef_trim, simplify=FALSE))
gs_mort_coef_trim <- do.call(rbind, replicate(n_rep, gs_mort_coef_trim, simplify=FALSE))
dial_mort_coef_trim <- do.call(rbind, replicate(n_rep, dial_mort_coef_trim, simplify=FALSE))

ddtx_shape_trim      <- rep(ddtx_shape_trim, n_rep)
ldtx_shape_trim      <- rep(ldtx_shape_trim, n_rep)
mort_shape_trim      <- rep(mort_shape_trim, n_rep)
remove_shape_trim    <- rep(remove_shape_trim, n_rep)
gl_shape_trim        <- rep(gl_shape_trim, n_rep)
gs_mort_shape_trim   <- rep(gs_mort_shape_trim, n_rep)
dial_mort_shape_trim <- rep(dial_mort_shape_trim, n_rep)

n <- 64 * n_rep

out <- foreach(i = 1:n, .combine = "comb", .multicombine = TRUE,
               .packages = c("MASS", "copula", "psych",
                             "dplyr", "survival", "dtplyr")) %dopar% {
  print(i)
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

  m_base_costs_qalys <- f_costs_qalys(data = tmp_baseline, r = disc)

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

      m_policy_costs_qalys <- f_costs_qalys(data = tmp_policy, r = disc)
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
v_id <- rep(seq.int(1:64), 100)
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
serv_path <- "Z:/Deterministic Results/"
lapply(names(l_policy_df), function(dname) {
  write.csv(l_policy_df[[dname]],
            file = paste0(serv_path, dname, ".csv"), row.names = FALSE)
})
write.csv(base_df, file = paste0(serv_path, "base_df_03_11_24.csv"), row.names = FALSE)

### Start here to use previous data ###
serv_path <- "results/Base Case/2024_03_12/" # run if on laptop
pdf_suffix <- paste(format(Sys.time(), "%Y%m%d"), "pdf", sep = ".")
png_suffix <- paste(format(Sys.time(), "%Y%m%d"), "png", sep = ".")

n_rep <- 10

base_df <- read.csv(paste0(serv_path, "base_df.csv"))

l_policy_df_csv <- c("policy_1_cq_counts",  "policy_2_cq_counts",
                     "policy_3_cq_counts",  "policy_4_cq_counts",
                     "policy_5_cq_counts",  "policy_6_cq_counts",
                     "policy_7_cq_counts",  "policy_8_cq_counts",
                     "policy_9_cq_counts",  "policy_10_cq_counts",
                     "policy_11_cq_counts", "policy_12_cq_counts",
                     "policy_13_cq_counts", "policy_14_cq_counts",
                     "policy_15_cq_counts", "policy_16_cq_counts",
                     "policy_17_cq_counts", "policy_18_cq_counts",
                     "policy_19_cq_counts", "policy_20_cq_counts")
for (i in seq_along(l_policy_df_csv)) {
  assign(l_policy_df_csv[i], read.csv(paste0(serv_path, l_policy_df_csv[i],
                                             ".csv")))
}

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

### Analysis starts here ###
event_counts_names <- c("ddtx_count", "ldtx_count", "wait_mort_count",
                        "removal_count", "dgf_count", "gl_count",
                        "ddtx_count_diab0", "ldtx_count_diab0",
                        "wait_mort_count_diab0", "removal_count_diab0",
                        "dgf_count_diab0", "gl_count_diab0", "ddtx_count_diab1",
                        "ldtx_count_diab1", "wait_mort_count_diab1",
                        "removal_count_diab1", "dgf_count_diab1",
                        "gl_count_diab1", "ddtx_count_white",
                        "ldtx_count_white", "wait_mort_count_white",
                        "removal_count_white", "dgf_count_white",
                        "gl_count_white", "ddtx_count_black",
                        "ldtx_count_black", "wait_mort_count_black",
                        "removal_count_black", "dgf_count_black",
                        "gl_count_black", "ddtx_count_hisp", "ldtx_count_hisp",
                        "wait_mort_count_hisp", "removal_count_hisp",
                        "dgf_count_hisp", "gl_count_hisp", "ddtx_count_other",
                        "ldtx_count_other", "wait_mort_count_other",
                        "removal_count_other", "dgf_count_other",
                        "gl_count_other", "ddtx_count_65", "ldtx_count_65",
                        "wait_mort_count_65", "removal_count_65",
                        "dgf_count_65", "gl_count_65", "ddtx_count_70",
                        "ldtx_count_70", "wait_mort_count_70",
                        "removal_count_70", "dgf_count_70", "gl_count_70",
                        "ddtx_count_75", "ldtx_count_75", "wait_mort_count_75",
                        "removal_count_75", "dgf_count_75", "gl_count_75")

v_colnames <- c("tot_cost_health", "tot_cost_societal", "tot_qalys",
                "tot_cost_health_diab_0", "tot_cost_societal_diab_0",
                "tot_qalys_diab_0", "tot_cost_health_diab_1",
                "tot_cost_societal_diab_1", "tot_qalys_diab_1",
                "tot_cost_health_white", "tot_cost_societal_white",
                "tot_qalys_white", "tot_cost_health_black",
                "tot_cost_societal_black", "tot_qalys_black",
                "tot_cost_health_hispanic", "tot_cost_societal_hispanic",
                "tot_qalys_hispanic", "tot_cost_health_other",
                "tot_cost_societal_other", "tot_qalys_other",
                "tot_cost_health_65", "tot_cost_societal_65", "tot_qalys_65",
                "tot_cost_health_70", "tot_cost_societal_70", "tot_qalys_70",
                "tot_cost_health_75", "tot_cost_societal_75", "tot_qalys_75",
                event_counts_names)
v_colnames_qalys <- v_colnames[str_detect(v_colnames, "qalys")]
v_colnames_counts <- v_colnames[str_detect(v_colnames, "count")]


df_sample <- as.data.frame(sample_df)
n_sample <- 100000
n_white <- n_sample - sum(df_sample$black) - sum(df_sample$hispanic) -
  sum(df_sample$other)
n_diab_0 <- n_sample - sum(df_sample$diab_stat)


scale_all   <- function(x) (x / (n_sample * n_rep))
scale_diab0 <- function(x) (x / (n_diab_0 * n_rep))
scale_diab1 <- function(x) (x / (sum(df_sample$diab_stat) * n_rep))
scale_white <- function(x) (x / (n_white * n_rep))
scale_black <- function(x) (x / (sum(df_sample$black) * n_rep))
scale_hisp  <- function(x) (x / (sum(df_sample$hispanic) * n_rep))
scale_other <- function(x) (x / (sum(df_sample$other) * n_rep))
scale_65    <- function(x) (x / (sum(df_sample$can_age_at_listing_c < 5) * n_rep))
scale_70    <- function(x) {
  x / (sum(df_sample$can_age_at_listing_c >= 5 &
             df_sample$can_age_at_listing_c < 10) * n_rep)
}
scale_75    <- function(x) {
  x / (sum(df_sample$can_age_at_listing_c >= 10) * n_rep)
}
scale_q <- function(x) (x / 12)

get_clean_outcomes <- function(data, v_names) {
  data <- as.data.frame(data)
  colnames(data) <- v_names
  data <- as.data.frame(sapply(data, as.numeric))
  return(data)
}

get_cea_results <- function(df, v_names = v_colnames,
                            v_qaly_names = v_colnames_qalys) {
  df_clean   <- get_clean_outcomes(data = df, v_names = v_names)

  df_pp <- df_clean %>%
    dplyr::mutate(across(c("tot_cost_health", "tot_cost_societal", "tot_qalys"),
                         scale_all)) %>%
    dplyr::mutate(across(c("tot_cost_health_diab_0", "tot_cost_societal_diab_0",
                           "tot_qalys_diab_0"),
                         scale_diab0)) %>%
    dplyr::mutate(across(c("tot_cost_health_diab_1", "tot_cost_societal_diab_1",
                           "tot_qalys_diab_1"),
                         scale_diab1)) %>%
    dplyr::mutate(across(c("tot_cost_health_white", "tot_cost_societal_white",
                           "tot_qalys_white"),
                         scale_white)) %>%
    dplyr::mutate(across(c("tot_cost_health_black", "tot_cost_societal_black",
                           "tot_qalys_black"),
                         scale_black)) %>%
    dplyr::mutate(across(c("tot_cost_health_hispanic",
                           "tot_cost_societal_hispanic", "tot_qalys_hispanic"),
                         scale_hisp)) %>%
    dplyr::mutate(across(c("tot_cost_health_other",
                           "tot_cost_societal_other", "tot_qalys_other"),
                         scale_other)) %>%
    dplyr::mutate(across(c("tot_cost_health_65",
                           "tot_cost_societal_65", "tot_qalys_65"),
                         scale_65)) %>%
    dplyr::mutate(across(c("tot_cost_health_70",
                           "tot_cost_societal_70", "tot_qalys_70"),
                         scale_70)) %>%
    dplyr::mutate(across(c("tot_cost_health_75",
                           "tot_cost_societal_75", "tot_qalys_75"),
                         scale_75)) %>%
    dplyr::mutate(across(v_qaly_names, scale_q)) %>%
    dplyr::mutate(dgf_per_10k_tx = (dgf_count/ddtx_count)*10000) %>%
    dplyr::mutate(gl_per_10k_tx = (gl_count/ddtx_count)*10000) %>%
    as.matrix()

  return(df_pp)
}

v_colnames <- c("id", v_colnames)
base_df_pp <- get_cea_results(df = base_df, v_names = v_colnames,
                              v_qaly_names = v_colnames_qalys)
# base_df_pp <- cbind(seq.int(1, 64), base_df_pp)

l_policy_df_clean <- lapply(l_policy_df, function(x) get_cea_results(x))
# l_policy_df_clean <- lapply(l_policy_df_clean, function(x) {x <- x %>% select(-c(id))})
l_inc <- lapply(l_policy_df_clean, function(x) x - base_df_pp)
wtp <- 100000
for (i in seq(l_inc)) {
  l_inc[[i]] <- as.data.frame(l_inc[[i]]) %>%
    mutate(icer_health_all = tot_cost_health / tot_qalys) %>%
    mutate(icer_societal_all = tot_cost_societal / tot_qalys) %>%
    mutate(icer_health_diab_0 = tot_cost_health_diab_0 / tot_qalys_diab_0) %>%
    mutate(icer_societal_diab_0 = tot_cost_societal_diab_0 /
             tot_qalys_diab_0) %>%
    mutate(icer_health_diab_1 = tot_cost_health_diab_1 / tot_qalys_diab_1) %>%
    mutate(icer_societal_diab_1 = tot_cost_societal_diab_1 /
             tot_qalys_diab_1) %>%
    mutate(icer_health_white = tot_cost_health_white / tot_qalys_white) %>%
    mutate(icer_societal_white = tot_cost_societal_white / tot_qalys_white) %>%
    mutate(icer_health_black = tot_cost_health_black / tot_qalys_black) %>%
    mutate(icer_societal_black = tot_cost_societal_black / tot_qalys_black) %>%
    mutate(icer_health_hispanic = tot_cost_health_hispanic /
             tot_qalys_hispanic) %>%
    mutate(icer_societal_hispanic = tot_cost_societal_hispanic /
             tot_qalys_hispanic) %>%
    mutate(icer_health_other = tot_cost_health_other / tot_qalys_other) %>%
    mutate(icer_societal_other = tot_cost_societal_other / tot_qalys_other) %>%
    mutate(icer_health_65 = tot_cost_health_65 / tot_qalys_65) %>%
    mutate(icer_societal_65 = tot_cost_societal_65 / tot_qalys_65) %>%
    mutate(icer_health_70 = tot_cost_health_70 / tot_qalys_70) %>%
    mutate(icer_societal_70 = tot_cost_societal_70 / tot_qalys_70) %>%
    mutate(icer_health_75 = tot_cost_health_75 / tot_qalys_75) %>%
    mutate(icer_societal_75 = tot_cost_societal_75 / tot_qalys_75) %>%
    mutate(nmb_health_all = (tot_qalys * wtp) - tot_cost_health) %>%
    mutate(nmb_societal_all = (tot_qalys * wtp) - tot_cost_societal) %>%
    as.data.frame()
}

l_inc_means <- lapply(l_inc, function(x) colMeans(x))
v_scenarios <- c("Status Quo", "5% Increase", "10% Increase", "15% Increase",
                 "20% Increase", "25% Increase")

df_inc_shift_05 <- data.frame(rbind(rep(0, 115), l_inc_means[[1]],
                                    l_inc_means[[5]], l_inc_means[[9]],
                                    l_inc_means[[13]], l_inc_means[[17]]))
df_inc_shift_10 <- data.frame(rbind(rep(0, 115), l_inc_means[[2]],
                                    l_inc_means[[6]], l_inc_means[[10]],
                                    l_inc_means[[14]], l_inc_means[[18]]))
df_inc_shift_15 <- data.frame(rbind(rep(0, 115), l_inc_means[[3]],
                                    l_inc_means[[7]], l_inc_means[[11]],
                                    l_inc_means[[15]], l_inc_means[[19]]))
df_inc_shift_20 <- data.frame(rbind(rep(0, 115), l_inc_means[[4]],
                                    l_inc_means[[8]], l_inc_means[[12]],
                                    l_inc_means[[16]], l_inc_means[[20]]))

df_inc_shift_05$v_scenarios <- v_scenarios
df_inc_shift_10$v_scenarios <- v_scenarios
df_inc_shift_15$v_scenarios <- v_scenarios
df_inc_shift_20$v_scenarios <- v_scenarios

# df_policy <- as.data.frame(do.call(rbind, l_policy_means))
# df_all <- cbind(rbind(base_means, df_policy), seq(0,20),
#                 c("0%", rep("05%", 4), rep("10%", 4), rep("15%", 4),
#                   rep("20%", 4), rep("25%", 4)),
#                 c("0%", rep(c("05%", "10%", "15%", "20%"), 5)))
# names(df_all)[31:33] <- c("id", "tx_bump", "kdpi_shift")
# df_shift <- split(df_all, "kdpi_shift")

# df_all$tx_bump <- as.factor(df_all$tx_bump)
# df_all$kdpi_shift <- as.factor(df_all$kdpi_shift)

# Quick chart to plot all costs and qalys (not incremental)
library(ggplot2)
library(dampack)
# serv_path <- "C:/Users/mbkauf/Documents/Calibration/Recalibration/"
# ce_plane_health <- ggplot(data = df_all, aes(x = tot_cost_health,
#                                              y = tot_qalys,
#                                              color = tx_bump)) +
#   facet_wrap(.~kdpi_shift) +
#   geom_point()
# ce_plane_health
#
# ice_plane_health <- ggplot(data = df_inc_shift_05,
#                            aes(x = tot_cost_health,
#                                y = tot_qalys,
#                                color = v_scenarios)) +
#   geom_point() +
#   geom_line()
# ice_plane_health

icer_tab <- calculate_icers(cost = df_inc_shift_05$tot_cost_health,
                            effect = df_inc_shift_05$tot_qalys,
                            strategies = v_scenarios)
ice_plane <- plot(icer_tab, label = "all") +
  coord_flip() + theme_bw() + ylab("Inc. Cost ($)") +
  xlab("Inc. Effectiveness (QALYs)") +
  ggtitle("Cost-Effectiveness Frontier - Healthcare Sector Perspective") +
  annotate("text", x = 0.31, y = 1500, label = "$5,500 $/QALY", angle = 39)
set.seed(2)
# Save as pdf and png
# ggsave(filename = paste0(serv_path, "ice_plane_post_" , pdf_suffix),
#        plot = ice_plane)
# ggsave(filename = paste0(serv_path, "ice_plane_post_" , png_suffix),
#        plot = ice_plane)

# Test Boxplot Results


### Base Case Figures for Paper
library(ggplot2)
library(dampack)
library(ggpubr)
library(ggsci)
library(grid)
library(gridExtra)
library(gtable)

fig_path <- "results/Base Case/"

df_05 <- l_inc[[1]] %>%
  mutate(scenario = "5%")
df_10 <- l_inc[[5]] %>%
  mutate(scenario = "10%")
df_15 <- l_inc[[9]] %>%
  mutate(scenario = "15%")
df_20 <- l_inc[[13]] %>%
  mutate(scenario = "20%")
df_25 <- l_inc[[17]] %>%
  mutate(scenario = "25%")

df_shift_05_all <- rbind(df_05, df_10, df_15, df_20, df_25)
df_shift_05_all$scenario <- factor(df_shift_05_all$scenario,
                                   levels = c("5%", "10%",
                                              "15%", "20%",
                                              "25%"))
df_shift_05_all_long <- gather(df_shift_05_all, nmb, value,
                               nmb_health_all:nmb_societal_all)
df_shift_05_all_long$nmb <- factor(df_shift_05_all_long$nmb,
                                   levels = c("nmb_health_all",
                                              "nmb_societal_all"),
                                   labels = c("Healthcare Sector",
                                              "Modified Healthcare Sector"))

## Figure 2: Incremental cost-effectiveness frontier for base-case results from
##           healthcare sector perspective
# icer_tab <- calculate_icers(cost = df_inc_shift_05$tot_cost_health,
#                             effect = df_inc_shift_05$tot_qalys,
#                             strategies = v_scenarios)
#
# fig2a <- plot(icer_tab, label = "all") +
#   coord_flip() + theme_bw() + ylab("Inc. Cost ($)") +
#   xlab("Inc. Effectiveness (QALYs)") +
#   ggtitle("Cost-Effectiveness Frontier - Healthcare Sector Perspective") +
#   annotate("text", x = 0.3, y = 1400, label = "$5,500 $/QALY", angle = 43.5)
# set.seed(2)
# ggsave(filename = paste0(serv_path, "ice_plane_post_" , pdf_suffix),
#        plot = fig2a)
# ggsave(filename = paste0(serv_path, "ice_plane_post_" , png_suffix),
#        plot = fig2a)
#
# icer_tab_soc <- calculate_icers(cost = df_inc_shift_05$tot_cost_societal,
#                             effect = df_inc_shift_05$tot_qalys,
#                             strategies = v_scenarios)
# fig2b <- plot(icer_tab_soc, label = "all") +
#   coord_flip() + theme_bw() + ylab("Inc. Cost ($)") +
#   scale_y_continuous(labels = scales::comma) +
#   xlab("Inc. Effectiveness (QALYs)") +
#   ggtitle("Cost-Effectiveness Frontier - Societal Sector Perspective")
#
# fig_ice_health <- fig2a +
#   labs(title = "Healthcare Sector Perspective", tag = "A")
#
# fig_ice_soc    <- fig2b +
#   labs(title = "Societal Perspective", tag = "B")
#
# ice_combined <- ggarrange(fig_ice_health, fig_ice_soc, nrow = 1,
#                           legend = "bottom")
# set.seed(2)
# ggsave(filename = paste0(serv_path, "ice_combined_", pdf_suffix),
#        plot = ice_combined, width = 10, height = 6)
# ggsave(filename = paste0(serv_path, "ice_combined_", png_suffix),
#        plot = ice_combined, width = 10, height = 6)

## Figure 2 custom
format.qaly  <- function(x, ...) {
  paste0("$", formatC(round(as.numeric(x), digits=-2), format="f", digits=0, big.mark=","),
         "/QALY")
}

icer_tab <- calculate_icers(cost = df_inc_shift_05$tot_cost_health,
                            effect = df_inc_shift_05$tot_qalys,
                            strategies = v_scenarios)
icer_tab$Strategy <- factor(icer_tab$Strategy,
                            levels = v_scenarios)
icer_health <- format.qaly(icer_tab[2,6])
fig2a <- plot_icers_mod(icer_tab, label = "none") +
  coord_flip() + theme_bw() + ylab("Inc. Cost ($)") +
  xlab("Inc. Effectiveness (QALYs)") +
  annotate("text", x = 0.30, y = 2200, size = 3,
           label = icer_health, angle = 45)

icer_tab_soc <- calculate_icers(cost = df_inc_shift_05$tot_cost_societal,
                                effect = df_inc_shift_05$tot_qalys,
                                strategies = v_scenarios)
icer_tab_soc$Strategy <- factor(icer_tab_soc$Strategy,
                                levels = v_scenarios)
fig2b <- plot_icers_mod(icer_tab_soc, label = "none") +
  coord_flip() + theme_bw() + ylab("Inc. Cost ($)") +
  scale_y_continuous(labels = scales::comma) +
  xlab("Inc. Effectiveness (QALYs)")

fig_ice_health <- fig2a +
  labs(title = "Healthcare Sector", tag = "A")

fig_ice_soc    <- fig2b +
  labs(title = "Modified Healthcare Sector", tag = "B")

ice_combined <- ggarrange(fig_ice_health, fig_ice_soc, nrow = 1,
                          legend = "bottom", common.legend = TRUE)

ggsave(filename = paste0(serv_path, "ice_combined_", pdf_suffix),
       plot = ice_combined, width = 10, height = 6)
ggsave(filename = paste0(serv_path, "ice_combined_", png_suffix),
       plot = ice_combined, width = 10, height = 6)

## Figure 3: INMB by perspective
my_comparisons <- list(c("5%", "10%"),
                       c("5%", "15%"),
                       c("5%", "20%"))
fig3 <- ggplot(data = df_shift_05_all_long) +
  geom_boxplot(aes(x = scenario, y = value, fill = nmb)) +
  facet_grid(. ~ nmb) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  xlab("Percent Increase in Rate of Transplantation") +
  ylab("Incremental Net Monetary Benefit ($)") +
  scale_fill_nejm() +
  guides(fill = "none") +
  theme_bw() +
  scale_y_continuous(labels = scales::comma) +
  stat_compare_means(aes(x = scenario, y = value),
                     method = "anova")
ggsave(filename = paste0(serv_path, "inmb_perspective_", pdf_suffix),
       plot = fig3, height = 6, width = 8)
ggsave(filename = paste0(serv_path, "inmb_perspective_", png_suffix),
       plot = fig3, height = 6, width = 8)

## Figure 4: INMB by key subgroup
# Generate df, using inc shift 5%
df_fig4 <- df_shift_05_all %>%
  dplyr::mutate(inmb_health_65 = (tot_qalys_65 * wtp) - tot_cost_health_65) %>%
  dplyr::mutate(inmb_soc_65 = (tot_qalys_65 * wtp) - tot_cost_societal_65) %>%
  dplyr::mutate(inmb_health_70 = (tot_qalys_70 * wtp) - tot_cost_health_70) %>%
  dplyr::mutate(inmb_soc_70 = (tot_qalys_70 * wtp) - tot_cost_societal_70) %>%
  dplyr::mutate(inmb_health_75 = (tot_qalys_75 * wtp) - tot_cost_health_75) %>%
  dplyr::mutate(inmb_soc_75 = (tot_qalys_75 * wtp) - tot_cost_societal_75) %>%
  dplyr::mutate(inmb_health_diab_0 = (tot_qalys_diab_0 * wtp) - tot_cost_health_diab_0) %>%
  dplyr::mutate(inmb_soc_diab_0 = (tot_qalys_diab_0 * wtp) - tot_cost_societal_diab_0) %>%
  dplyr::mutate(inmb_health_diab_1 = (tot_qalys_diab_1 * wtp) - tot_cost_health_diab_1) %>%
  dplyr::mutate(inmb_soc_diab_1 = (tot_qalys_diab_1 * wtp) - tot_cost_societal_diab_1) %>%
  dplyr::mutate(inmb_health_white = (tot_qalys_white * wtp) - tot_cost_health_white) %>%
  dplyr::mutate(inmb_soc_white = (tot_qalys_white * wtp) - tot_cost_societal_white) %>%
  dplyr::mutate(inmb_health_black = (tot_qalys_black * wtp) - tot_cost_health_black) %>%
  dplyr::mutate(inmb_soc_black = (tot_qalys_black * wtp) - tot_cost_societal_black) %>%
  dplyr::mutate(inmb_health_hisp = (tot_qalys_hispanic * wtp) - tot_cost_health_hispanic) %>%
  dplyr::mutate(inmb_soc_hisp = (tot_qalys_hispanic * wtp) - tot_cost_societal_hispanic) %>%
  dplyr::mutate(inmb_health_other = (tot_qalys_other * wtp) - tot_cost_health_other) %>%
  dplyr::mutate(inmb_soc_other = (tot_qalys_other * wtp) - tot_cost_societal_other) %>%
  gather(inmb, value, inmb_health_65:inmb_soc_other) %>%
  dplyr::select(c(scenario, inmb, value)) %>%
  dplyr::mutate(subgrp = case_when(
    str_detect(inmb, "70") ~ "age",
    str_detect(inmb, "5") ~ "age",
    str_detect(inmb, "diab") ~ "diab",
    str_detect(inmb, "white") ~ "race",
    str_detect(inmb, "black") ~ "race",
    str_detect(inmb, "hisp") ~ "race",
    str_detect(inmb, "other") ~ "race",
  )) %>%
  dplyr::mutate(subgrp_spec = case_when(
    str_detect(inmb, "65") ~ "65-69",
    str_detect(inmb, "70") ~ "70-74",
    str_detect(inmb, "75") ~ "75+",
    str_detect(inmb, "diab_0") ~ "No",
    str_detect(inmb, "diab_1") ~ "Yes",
    str_detect(inmb, "white") ~ "NH White",
    str_detect(inmb, "black") ~ "NH Black",
    str_detect(inmb, "hisp") ~ "Hispanic",
    str_detect(inmb, "other") ~ "NH Other",
  )) %>%
  dplyr::mutate(pers = case_when(
    str_detect(inmb, "health") ~ "Healthcare Sector",
    str_detect(inmb, "soc") ~ "Modified Healthcare Sector",
  )) %>%
  filter(scenario == "25%")
l_df_fig4 <- split(df_fig4, df_fig4$subgrp)

my_comparisons <- list(c("65-69", "75+"))
fig4a <- ggplot(data = l_df_fig4[[1]],
                aes(x = subgrp_spec, y = value)) +
  geom_boxplot(aes(fill = pers)) +
  facet_grid(. ~ pers) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  xlab("Age at Listing") +
  scale_fill_nejm() +
  guides(fill = "none") +
  theme_bw() +
  scale_y_continuous(labels = scales::comma,
                     breaks = c(0, 40000, 80000, 120000),
                     limits = c(0, 140000)) +
  stat_compare_means(comparisons = my_comparisons,
                     aes(x = subgrp_spec, y = value,
                         label = after_stat(p.signif))) +
  theme(axis.title.y = element_blank())

my_comparisons <- list(c("Yes", "No"))
fig4b <- ggplot(data = l_df_fig4[[2]],
                aes(x = subgrp_spec, y = value)) +
  geom_boxplot(aes(fill = pers)) +
  facet_grid(. ~ pers) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  xlab("Diabetes History") +
  scale_fill_nejm() +
  guides(fill = "none") +
  theme_bw() +
  scale_y_continuous(labels = scales::comma,
                     breaks = c(0, 40000, 80000, 120000),
                     limits = c(0, 140000)) +
  stat_compare_means(comparisons = my_comparisons,
                     aes(x = subgrp_spec, y = value,
                         label = after_stat(p.signif))) +
  theme(axis.title.y = element_blank())

my_comparisons <- list(c("NH White", "NH Black"))
fig4c <- ggplot(data = l_df_fig4[[3]],
                aes(x = subgrp_spec, y = value)) +
  geom_boxplot(aes(fill = pers)) +
  facet_grid(. ~ pers) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  xlab("Race/Ethnicity") +
  scale_fill_nejm() +
  guides(fill = "none") +
  theme_bw() +
  scale_y_continuous(labels = scales::comma,
                     breaks = c(0, 40000, 80000, 120000),
                     limits = c(0, 140000)) +
  stat_compare_means(comparisons = my_comparisons,
                     aes(x = subgrp_spec, y = value,
                         label = after_stat(p.signif))) +
  theme(axis.title.y = element_blank())

gta <- ggplotGrob(fig4a)
gtb <- ggplotGrob(fig4b)
gtc <- ggplotGrob(fig4c)

new_width <- unit.pmax(gta$widths[2:3], gtb$widths[2:3], gtc$widths[2:3])

gta$widths[2:3] <- as.list(new_width)
gtb$widths[2:3] <- as.list(new_width)
gtc$widths[2:3] <- as.list(new_width)

fig4 <- grid.arrange(gta, gtb, gtc,
                     ncol = 3,
                     left = "Incremental Net Monetary Benefit ($)")

ggsave(filename = paste0(serv_path, "inmb_key_subgrp_", pdf_suffix),
       plot = fig4, width = 12, height = 4)
ggsave(filename = paste0(serv_path, "inmb_key_subgrp_", png_suffix),
       plot = fig4, width = 12, height = 4)


## INMB for each subgroup
# AGE
fig4_age <- ggplot(data = l_df_fig4[[1]],
                   aes(x = subgrp_spec, y = value)) +
  geom_boxplot(aes(fill = pers)) +
  facet_grid(. ~ pers) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  xlab("Age at Listing") +
  ylab("Incremental Net Monetary Benefit ($)") +
  scale_fill_nejm() +
  guides(fill = "none") +
  theme_bw() +
  scale_y_continuous(labels = scales::comma,
                     breaks = c(0, 40000, 80000, 120000),
                     limits = c(0, 140000)) +
  stat_compare_means(aes(x = subgrp_spec, y = value),
                     method = "anova", label.y = 10000)

my_comparisons <- list(c("Yes", "No"))
fig4_diab <- ggplot(data = l_df_fig4[[2]],
                aes(x = subgrp_spec, y = value)) +
  geom_boxplot(aes(fill = pers)) +
  facet_grid(. ~ pers) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  xlab("Diabetes History") +
  ylab("Incremental Net Monetary Benefit ($)") +
  scale_fill_nejm() +
  guides(fill = "none") +
  theme_bw() +
  scale_y_continuous(labels = scales::comma,
                     breaks = c(0, 40000, 80000, 120000),
                     limits = c(0, 140000)) +
  stat_compare_means(comparisons = my_comparisons,
                     aes(x = subgrp_spec, y = value,
                         label = after_stat(p.signif)))

my_comparisons <- list(c("NH White", "NH Black"),
                       c("NH White", "NH Other"),
                       c("NH White", "Hispanic"))
fig4_race <- ggplot(data = l_df_fig4[[3]],
                    aes(x = subgrp_spec, y = value)) +
  geom_boxplot(aes(fill = pers)) +
  facet_grid(. ~ pers) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  xlab("Race/Ethnicity") +
  ylab("Incremental Net Monetary Benefit ($)") +
  scale_fill_nejm() +
  guides(fill = "none") +
  theme_bw() +
  scale_y_continuous(labels = scales::comma,
                     breaks = c(0, 40000, 80000, 120000),
                     limits = c(0, 140000)) +
  stat_compare_means(method = "anova", label.y = 10000)+ # Add global p-value
  stat_compare_means(aes(label = after_stat(p.signif)),
                     method = "t.test", ref.group = "NH White")


ggsave(filename = paste0(out_path, "inmb_age.pdf"), plot = fig4_age,
       width = 8, height = 4)
ggsave(filename = paste0(out_path, "inmb_age.png"), plot = fig4_age,
       width = 8, height = 4)
ggsave(filename = paste0(out_path, "inmb_diab.pdf"), plot = fig4_diab,
       width = 8, height = 4)
ggsave(filename = paste0(out_path, "inmb_diab.png"), plot = fig4_diab,
       width = 8, height = 4)
ggsave(filename = paste0(out_path, "inmb_race.pdf"), plot = fig4_race,
       width = 8, height = 4)
ggsave(filename = paste0(out_path, "inmb_race.png"), plot = fig4_race,
       width = 8, height = 4)

## Figure 4: ICE Plane by key subgroup
# Overlay each subgroup on the same plot
# df_fig4 <- df_shift_05_all %>%
#   gather(icer, value, icer_health_diab_0:icer_societal_75) %>%
#   dplyr::select(c(scenario, icer, value)) %>%
#   dplyr::mutate(subgrp = case_when(
#     str_detect(icer, "5") ~ "age",
#     str_detect(icer, "7") ~ "age",
#     str_detect(icer, "diab") ~ "diab",
#     str_detect(icer, "white") ~ "race",
#     str_detect(icer, "black") ~ "race",
#   )) %>%
#   dplyr::mutate(subgrp_spec = case_when(
#     str_detect(icer, "65") ~ "65-69",
#     str_detect(icer, "70") ~ "70-74",
#     str_detect(icer, "75") ~ "75+",
#     str_detect(icer, "diab_0") ~ "No",
#     str_detect(icer, "diab_1") ~ "Yes",
#     str_detect(icer, "white") ~ "NH White",
#     str_detect(icer, "black") ~ "NH Black",
#     str_detect(icer, "hisp") ~ "Hispanic",
#     str_detect(icer, "other") ~ "NH Other",
#   )) %>%
#   dplyr::mutate(pers = case_when(
#     str_detect(icer, "health") ~ "Healthcare Sector",
#     str_detect(icer, "societal") ~ "Societal",
#   ))
# l_df_fig4 <- split(df_fig4, df_fig4$subgrp)

## AGE
icer_tab_65 <- calculate_icers(cost = df_inc_shift_05$tot_cost_health_65,
                               effect = df_inc_shift_05$tot_qalys_65,
                               strategies = v_scenarios)
icer_tab_65$subgrp <- "65-69"
icer_tab_70 <- calculate_icers(cost = df_inc_shift_05$tot_cost_health_70,
                               effect = df_inc_shift_05$tot_qalys_70,
                               strategies = v_scenarios)
icer_tab_70$subgrp <- "70-74"
icer_tab_75 <- calculate_icers(cost = df_inc_shift_05$tot_cost_health_75,
                               effect = df_inc_shift_05$tot_qalys_75,
                               strategies = v_scenarios)
icer_tab_75$subgrp <- "75+"
icer_tab_age_health <- rbind(icer_tab_65, icer_tab_70, icer_tab_75)
icer_tab_age_health$Strategy <- factor(icer_tab_age_health$Strategy,
                                levels = v_scenarios)
icer_tab_age_health$Perspective <- "Healthcare Sector"


icer_health_65 <- format.qaly(icer_tab_age_health[2, 6])
icer_health_70 <- format.qaly(icer_tab_age_health[8, 6])
icer_health_75 <- format.qaly(icer_tab_age_health[14, 6])

# fig4a <- plot_icers_mod_subgrp_3(icer_tab_age, label = "none") +
#   scale_y_continuous(labels = scales::comma, n.breaks = 6) +
#   coord_flip() + theme_bw() + ylab("Inc. Cost ($)") +
#   xlab("Inc. Effectiveness (QALYs)") +
#   annotate("text", x = 0.60, y = 0, size = 3, color = "#BC3C29FF",
#            label = icer_health_65, fontface = 2, hjust = 0) +
#   annotate("text", x = 0.57, y = 0, size = 3, color = "#0072B5FF",
#            label = icer_health_70, fontface = 2, hjust = 0) +
#   annotate("text", x = 0.54, y = 0, size = 3, color = "#E18727FF",
#            label = icer_health_75, fontface = 2, hjust = 0)

icer_tab_65 <- calculate_icers(cost = df_inc_shift_05$tot_cost_societal_65,
                               effect = df_inc_shift_05$tot_qalys_65,
                               strategies = v_scenarios)
icer_tab_65$subgrp <- "65-69"
icer_tab_70 <- calculate_icers(cost = df_inc_shift_05$tot_cost_societal_70,
                               effect = df_inc_shift_05$tot_qalys_70,
                               strategies = v_scenarios)
icer_tab_70$subgrp <- "70-74"
icer_tab_75 <- calculate_icers(cost = df_inc_shift_05$tot_cost_societal_75,
                               effect = df_inc_shift_05$tot_qalys_75,
                               strategies = v_scenarios)
icer_tab_75$subgrp <- "75+"
icer_tab_age_health_mod <- rbind(icer_tab_65, icer_tab_70, icer_tab_75)
icer_tab_age_health_mod$Strategy <- factor(icer_tab_age_health_mod$Strategy,
                                       levels = v_scenarios)
icer_tab_age_health_mod$Perspective <- "Modified Healthcare Sector"

icer_societal_65 <- format.qaly(icer_tab_age_health_mod[2, 6])
icer_societal_70 <- format.qaly(icer_tab_age_health_mod[8, 6])
icer_societal_75 <- format.qaly(icer_tab_age_health_mod[14, 6])

icer_tab_age_all <- rbind(icer_tab_age_health, icer_tab_age_health_mod)

dat_text <- data.frame(
  subgrp = c("65-69", "70-74", "75+", "65-69", "70-74", "75+"),
  Perspective = c("Healthcare Sector", "Healthcare Sector", "Healthcare Sector",
                  "Modified Healthcare Sector", "Modified Healthcare Sector",
                  "Modified Healthcare Sector"),
  label = c(icer_health_65, icer_health_70, icer_health_75, "", "", ""),
  x = c(1800, 1800, 1800, 0, 0, 0),
  y = c(0.32, 0.32, 0.32, 0, 0, 0)
)

fig4 <- plot_icers_mod_subgrp(icer_tab_age_all, label = "none") +
  geom_text(
    size    = 4,
    data    = dat_text,
    mapping = aes(x = x, y = y, label = label, angle = 30)
  ) +
    guides(shape = guide_legend("Strategy"),
           color = guide_legend("Age Group"),
           linetype = guide_legend(element_blank())) +
  theme(strip.text.x = element_text(size = 12),
        strip.text.y = element_text(size = 12))

# fig4b <- plot_icers_mod_subgrp_3(icer_tab_age, label = "none") +
#   scale_y_continuous(labels = scales::comma, n.breaks = 6) +
#   coord_flip() + theme_bw() + ylab("Inc. Cost ($)") +
#   xlab("Inc. Effectiveness (QALYs)")
#
# fig_ice_health <- fig4a +
#   labs(title = "Healthcare Sector Perspective", tag = "A") +
#   guides(shape = guide_legend("Strategy",
#                               theme = theme(legend.title = element_text(vjust = 0.85))),
#          color = guide_legend("Age Group"),
#          linetype = guide_legend(element_blank()))
#
# fig_ice_soc    <- fig4b +
#   labs(title = "Modified Healthcare Sector Perspective", tag = "B") +
#   guides(shape = guide_legend("Strategy",
#                               theme = theme(legend.title = element_text(vjust = 0.85))),
#          color = guide_legend("Age Group"),
#          linetype = guide_legend(element_blank()))
#
# ice_combined_age <- ggarrange(fig_ice_health, fig_ice_soc, nrow = 1,
#                               legend = "bottom", common.legend = TRUE)
# ice_combined_age

ggsave(filename = paste0(serv_path, "ice_combined_age_", pdf_suffix),
       plot = fig4, width = 12, height = 8)
ggsave(filename = paste0(serv_path, "ice_combined_age_", png_suffix),
       plot = fig4, width = 12, height = 8)

## Diabetes
icer_tab_diab_0 <- calculate_icers(cost = df_inc_shift_05$tot_cost_health_diab_0,
                               effect = df_inc_shift_05$tot_qalys_diab_0,
                               strategies = v_scenarios)
icer_tab_diab_0$subgrp <- "No"

icer_tab_diab_1 <- calculate_icers(cost = df_inc_shift_05$tot_cost_health_diab_1,
                                   effect = df_inc_shift_05$tot_qalys_diab_1,
                                   strategies = v_scenarios)
icer_tab_diab_1$subgrp <- "Yes"

icer_tab_diab_health <- rbind(icer_tab_diab_0, icer_tab_diab_1)
icer_tab_diab_health$Strategy <- factor(icer_tab_diab_health$Strategy,
                                 levels = v_scenarios)
icer_tab_diab_health$Perspective <- "Healthcare Sector"

icer_health_diab0 <- format.qaly(icer_tab_diab_health[2, 6])
icer_health_diab1 <- format.qaly(icer_tab_diab_health[8, 6])

# fig4a <- plot_icers_mod_subgrp_2(icer_tab_diab, label = "none") +
#   scale_y_continuous(labels = scales::comma, n.breaks = 7) +
#   coord_flip() + theme_bw() + ylab("Inc. Cost ($)") +
#   xlab("Inc. Effectiveness (QALYs)") +
#   annotate("text", x = 0.63, y = 0, size = 3, color = "#BC3C29FF",
#            label = icer_societal_diab0, fontface = 2, hjust = 0) +
#   annotate("text", x = 0.60, y = 0, size = 3, color = "#0072B5FF",
#            label = icer_societal_diab1, fontface = 2, hjust = 0)

icer_tab_diab_0 <- calculate_icers(cost = df_inc_shift_05$tot_cost_societal_diab_0,
                               effect = df_inc_shift_05$tot_qalys_diab_0,
                               strategies = v_scenarios)
icer_tab_diab_0$subgrp <- "No"
icer_tab_diab_1 <- calculate_icers(cost = df_inc_shift_05$tot_cost_societal_diab_1,
                               effect = df_inc_shift_05$tot_qalys_diab_1,
                               strategies = v_scenarios)
icer_tab_diab_1$subgrp <- "Yes"

icer_tab_diab_health_mod <- rbind(icer_tab_diab_0, icer_tab_diab_1)
icer_tab_diab_health_mod$Strategy <- factor(icer_tab_diab_health_mod$Strategy,
                                     levels = v_scenarios)
icer_tab_diab_health_mod$Perspective <- "Modified Healthcare Sector"
icer_tab_diab_all <- rbind(icer_tab_diab_health, icer_tab_diab_health_mod)

dat_text <- data.frame(
  subgrp = c("No", "Yes", "No", "Yes"),
  Perspective = c("Healthcare Sector", "Healthcare Sector",
                  "Modified Healthcare Sector", "Modified Healthcare Sector"),
  label = c(icer_health_diab0, icer_health_diab1, "", ""),
  x = c(1800, 2000, 0, 0),
  y = c(0.35, 0.26, 0, 0)
)

fig4_diab <- plot_icers_mod_subgrp(icer_tab_diab_all, label = "none") +
  geom_text(
    size    = 4,
    data    = dat_text,
    mapping = aes(x = x, y = y, label = label, angle = c(42, 30, 0, 0))
  ) +
  guides(shape = guide_legend("Strategy"),
         color = guide_legend("Diabetes History"),
         linetype = guide_legend(element_blank())) +
  theme(strip.text.x = element_text(size = 12),
        strip.text.y = element_text(size = 12))

# fig4b <- plot_icers_mod_subgrp_2(icer_tab_diab, label = "none") +
#   scale_y_continuous(labels = scales::comma, n.breaks = 7) +
#   coord_flip() + theme_bw() + ylab("Inc. Cost ($)") +
#   xlab("Inc. Effectiveness (QALYs)")
#
# fig_ice_health <- fig4a +
#   labs(title = "Healthcare Sector Perspective", tag = "A") +
#   guides(shape = guide_legend("Strategy",
#                               theme = theme(legend.title = element_text(vjust = 0.85))),
#          color = guide_legend("Diabetes History"),
#          linetype = guide_legend(element_blank()))
#
# fig_ice_soc    <- fig4b +
#   labs(title = "Modified Healthcare Sector Perspective", tag = "B")
#   guides(shape = guide_legend("Strategy",
#                               theme = theme(legend.title = element_text(vjust = 0.85))),
#          color = guide_legend("Diabetes History"),
#          linetype = guide_legend(element_blank()))
#
# ice_combined_diab <- ggarrange(fig_ice_health, fig_ice_soc, nrow = 1,
#                               legend = "bottom", common.legend = TRUE)
# ice_combined_diab

ggsave(filename = paste0(serv_path, "ice_combined_diab_", pdf_suffix),
       plot = fig4_diab, width = 12, height = 8)
ggsave(filename = paste0(serv_path, "ice_combined_diab_", png_suffix),
       plot = fig4_diab, width = 12, height = 8)

## RACE/ETHNICITY
icer_tab_white <- calculate_icers(cost = df_inc_shift_05$tot_cost_health_white,
                               effect = df_inc_shift_05$tot_qalys_white,
                               strategies = v_scenarios)
icer_tab_white$subgrp <- "NH White"
icer_tab_black <- calculate_icers(cost = df_inc_shift_05$tot_cost_health_black,
                               effect = df_inc_shift_05$tot_qalys_black,
                               strategies = v_scenarios)
icer_tab_black$subgrp <- "NH Black"
icer_tab_hispanic <- calculate_icers(cost = df_inc_shift_05$tot_cost_health_hispanic,
                               effect = df_inc_shift_05$tot_qalys_hispanic,
                               strategies = v_scenarios)
icer_tab_hispanic$subgrp <- "Hispanic"
icer_tab_other <- calculate_icers(cost = df_inc_shift_05$tot_cost_health_other,
                                     effect = df_inc_shift_05$tot_qalys_other,
                                     strategies = v_scenarios)
icer_tab_other$subgrp <- "NH Other"
icer_tab_race_health <- rbind(icer_tab_white, icer_tab_black,
                              icer_tab_hispanic, icer_tab_other)

icer_tab_race_health$Strategy <- factor(icer_tab_race_health$Strategy,
                                 levels = v_scenarios)
icer_tab_race_health$Perspective <- "Healthcare Sector"

icer_health_white <- format.qaly(icer_tab_race_health[2, 6])
icer_health_black <- format.qaly(icer_tab_race_health[8, 6])
icer_health_hisp  <- paste0(format.qaly(icer_tab_race_health[14, 6]), ", ",
                            format.qaly(icer_tab_race_health[15, 6]))
icer_health_other <- paste0(format.qaly(icer_tab_race_health[20, 6]), ", ",
                            format.qaly(icer_tab_race_health[21, 6]), ", ",
                            format.qaly(icer_tab_race_health[22, 6]))

icer_health_hisp_1 <- format.qaly(icer_tab_race_health[14, 6])
icer_health_hisp_2 <- format.qaly(icer_tab_race_health[15, 6])

icer_health_other_1 <- format.qaly(icer_tab_race_health[20, 6])
icer_health_other_2 <- format.qaly(icer_tab_race_health[21, 6])
icer_health_other_3 <- format.qaly(icer_tab_race_health[22, 6])

# fig4a <- plot_icers_mod_subgrp_4(icer_tab_race, label = "none") +
#   scale_y_continuous(labels = scales::comma, n.breaks = 7) +
#   coord_flip() + theme_bw() + ylab("Inc. Cost ($)") +
#   xlab("Inc. Effectiveness (QALYs)") +
#   annotate("text", x = 0.71, y = 0, size = 3, color = "#BC3C29FF",
#            label = icer_health_hisp, fontface = 2, hjust = 0) +
#   annotate("text", x = 0.68, y = 0, size = 3, color = "#0072B5FF",
#            label = icer_health_black, fontface = 2, hjust = 0) +
#   annotate("text", x = 0.65, y = 0, size = 3, color = "#E18727FF",
#            label = icer_health_other, fontface = 2, hjust = 0) +
#   annotate("text", x = 0.62, y = 0, size = 3, color = "#20854EFF",
#            label = icer_health_white, fontface = 2, hjust = 0)

icer_tab_white <- calculate_icers(cost = df_inc_shift_05$tot_cost_societal_white,
                                  effect = df_inc_shift_05$tot_qalys_white,
                                  strategies = v_scenarios)
icer_tab_white$subgrp <- "NH White"
icer_tab_black <- calculate_icers(cost = df_inc_shift_05$tot_cost_societal_black,
                                  effect = df_inc_shift_05$tot_qalys_black,
                                  strategies = v_scenarios)
icer_tab_black$subgrp <- "NH Black"
icer_tab_hispanic <- calculate_icers(cost = df_inc_shift_05$tot_cost_societal_hispanic,
                                     effect = df_inc_shift_05$tot_qalys_hispanic,
                                     strategies = v_scenarios)
icer_tab_hispanic$subgrp <- "Hispanic"
icer_tab_other <- calculate_icers(cost = df_inc_shift_05$tot_cost_societal_other,
                                  effect = df_inc_shift_05$tot_qalys_other,
                                  strategies = v_scenarios)
icer_tab_other$subgrp <- "NH Other"
icer_tab_race_health_mod <- rbind(icer_tab_white, icer_tab_black,
                              icer_tab_hispanic, icer_tab_other)

icer_tab_race_health_mod$Strategy <- factor(icer_tab_race_health_mod$Strategy,
                                        levels = v_scenarios)
icer_tab_race_health_mod$Perspective <- "Modified Healthcare Sector"

icer_tab_race_all <- rbind(icer_tab_race_health, icer_tab_race_health_mod)

dat_text <- data.frame(
  subgrp = c("NH White", "NH Black", "Hispanic", "Hispanic",
             "NH Other", "NH Other", "NH Other",
             "NH White", "NH Black", "Hispanic", "NH Other"),
  Perspective = c("Healthcare Sector", "Healthcare Sector", "Healthcare Sector",
                  "Healthcare Sector", "Healthcare Sector", "Healthcare Sector",
                  "Healthcare Sector",
                  "Modified Healthcare Sector", "Modified Healthcare Sector",
                  "Modified Healthcare Sector", "Modified Healthcare Sector"),
  label = c(icer_health_white, icer_health_black, icer_health_hisp_1,
            icer_health_hisp_2, icer_health_other_1, icer_health_other_2,
            icer_health_other_3, "", "", "", ""),
  x = c(2000, 2500, 2000, 4500, 200, 1100, 3000, 0, 0, 0, 0),
  y = c(0.3, 0.37, 0.31, 0.6, 0.1, 0.26, 0.54, 0, 0, 0, 0)
)

fig4_race <- plot_icers_mod_subgrp(icer_tab_race_all, label = "none") +
  geom_text(
    size    = 3,
    data    = dat_text,
    mapping = aes(x = x, y = y, label = label, angle = c(23, 23.5, 23, 22.5, 30,
                                                         29, 27.5, 0, 0, 0, 0))
  ) +
  guides(shape = guide_legend("Strategy"),
         color = guide_legend("Race/Ethnicity"),
         linetype = guide_legend(element_blank())) +
  theme(strip.text.x = element_text(size = 12),
        strip.text.y = element_text(size = 12))

# fig4b <- plot_icers_mod_subgrp_4(icer_tab_race, label = "none") +
#   scale_y_continuous(labels = scales::comma, n.breaks = 7) +
#   coord_flip() + theme_bw() + ylab("Inc. Cost ($)") +
#   xlab("Inc. Effectiveness (QALYs)")
#
# fig_ice_health <- fig4a +
#   labs(title = "Healthcare Sector Perspective", tag = "A") +
#   guides(shape = guide_legend("Strategy",
#                               theme = theme(legend.title = element_text(vjust = 0.85))),
#          color = guide_legend("Race/Ethnicity", ncol = 2,
#                               theme = theme(legend.title = element_text(vjust = 0.85))),
#          linetype = guide_legend(element_blank()))
#
# fig_ice_soc    <- fig4b +
#   labs(title = "Modified Healthcare Sector Perspective", tag = "B")
#   guides(shape = guide_legend("Strategy",
#                               theme = theme(legend.title = element_text(vjust = 0.85))),
#          color = guide_legend("Race/Ethnicity", ncol = 2,
#                               theme = theme(legend.title = element_text(vjust = 0.85))),
#          linetype = guide_legend(element_blank()))
#
# ice_combined_race <- ggarrange(fig_ice_health, fig_ice_soc, nrow = 1,
#                                legend = "bottom", common.legend = TRUE)
# ice_combined_race

ggsave(filename = paste0(serv_path, "ice_combined_race_", pdf_suffix),
       plot = fig4_race, width = 12, height = 9)
ggsave(filename = paste0(serv_path, "ice_combined_race_", png_suffix),
       plot = fig4_race, width = 12, height = 9)

## Figure 5: Scenario analysis for percentage shift in KDPI
# Generate dataset
df_25_shift_05 <- l_inc[[17]] %>%
  mutate(shift = "5%")
df_25_shift_10 <- l_inc[[18]] %>%
  mutate(shift = "10%")
df_25_shift_15 <- l_inc[[19]] %>%
  mutate(shift = "15%")
df_25_shift_20 <- l_inc[[20]] %>%
  mutate(shift = "20%")

df_25_shift_all <- rbind(df_25_shift_05, df_25_shift_10,
                         df_25_shift_15, df_25_shift_20)
df_25_shift_all$shift <- factor(df_25_shift_all$shift,
                                levels = c("5%", "10%",
                                           "15%", "20%"))

# INMB Figure
df_fig5 <- df_25_shift_all %>%
  gather(inmb, value, nmb_health_all:nmb_societal_all) %>%
  dplyr::mutate(pers = case_when(
    str_detect(inmb, "health") ~ "Healthcare Sector",
    str_detect(inmb, "societal") ~ "Modified Healthcare Sector",
  ))

fig5 <- ggplot(data = df_fig5) +
  geom_boxplot(aes(x = shift, y = value, fill = pers)) +
  facet_grid(. ~ pers) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  xlab("Percent of Lower Quality Donor Kidneys") +
  ylab("Incremental Net Monetary Benefit ($)") +
  scale_fill_nejm() +
  guides(fill = "none") +
  theme_bw() +
  scale_y_continuous(labels = scales::comma) +
  stat_compare_means(aes(x = shift, y = value),
                     method = "anova", label.y = 10000)
ggsave(filename = paste0(serv_path, "inmb_scenario_", pdf_suffix),
       plot = fig5, width = 8, height = 5)
ggsave(filename = paste0(serv_path, "inmb_scenario_", png_suffix),
       plot = fig5, width = 8, height = 5)

# ICE Plane Figure
icer_tab_5 <- calculate_icers(cost = df_inc_shift_05$tot_cost_health,
                              effect = df_inc_shift_05$tot_qalys,
                              strategies = v_scenarios)
icer_tab_5$subgrp <- "5%"
icer_tab_10 <- calculate_icers(cost = df_inc_shift_10$tot_cost_health,
                              effect = df_inc_shift_10$tot_qalys,
                              strategies = v_scenarios)
icer_tab_10$subgrp <- "10%"
icer_tab_15 <- calculate_icers(cost = df_inc_shift_15$tot_cost_health,
                               effect = df_inc_shift_15$tot_qalys,
                               strategies = v_scenarios)
icer_tab_15$subgrp <- "15%"

icer_tab_20 <- calculate_icers(cost = df_inc_shift_20$tot_cost_health,
                               effect = df_inc_shift_20$tot_qalys,
                               strategies = v_scenarios)
icer_tab_20$subgrp <- "20%"

icer_tab_shift_health <- rbind(icer_tab_5, icer_tab_10,
                               icer_tab_15, icer_tab_20)
icer_tab_shift_health$Strategy <- factor(icer_tab_shift_health$Strategy,
                                  levels = v_scenarios)
icer_tab_shift_health$subgrp <- factor(icer_tab_shift_health$subgrp,
                                levels = c("5%", "10%",
                                           "15%", "20%"))
icer_tab_shift_health$Perspective <- "Healthcare Sector"

icer_health_5  <- format.qaly(icer_tab_shift_health[2, 6])
icer_health_10 <- paste0(format.qaly(icer_tab_shift_health[8, 6]), ", ",
                         format.qaly(icer_tab_shift_health[9, 6]))
icer_health_15 <- format.qaly(icer_tab_shift_health[14, 6])
icer_health_20 <- format.qaly(icer_tab_shift_health[20, 6])

icer_health_10_1 <- format.qaly(icer_tab_shift_health[8, 6])
icer_health_10_2 <- format.qaly(icer_tab_shift_health[9, 6])

# fig5a <- plot_icers_mod_subgrp(icer_tab_shift, label = "none") +
#   scale_y_continuous(labels = scales::comma, n.breaks = 7) +
#   coord_flip() + theme_bw() + ylab("Inc. Cost ($)") +
#   xlab("Inc. Effectiveness (QALYs)") +
#   annotate("text", x = 0.56, y = 0, size = 3, color = "#BC3C29FF",
#          label = icer_health_5, fontface = 2, hjust = 0) +
#   annotate("text", x = 0.53, y = 0, size = 3, color = "#0072B5FF",
#            label = icer_health_10, fontface = 2, hjust = 0) +
#   annotate("text", x = 0.50, y = 0, size = 3, color = "#E18727FF",
#            label = icer_health_15, fontface = 2, hjust = 0) +
#   annotate("text", x = 0.47, y = 0, size = 3, color = "#20854EFF",
#            label = icer_health_20, fontface = 2, hjust = 0)


icer_tab_5 <- calculate_icers(cost = df_inc_shift_05$tot_cost_societal,
                              effect = df_inc_shift_05$tot_qalys,
                              strategies = v_scenarios)
icer_tab_5$subgrp <- "5%"
icer_tab_10 <- calculate_icers(cost = df_inc_shift_10$tot_cost_societal,
                               effect = df_inc_shift_10$tot_qalys,
                               strategies = v_scenarios)
icer_tab_10$subgrp <- "10%"
icer_tab_15 <- calculate_icers(cost = df_inc_shift_15$tot_cost_societal,
                               effect = df_inc_shift_15$tot_qalys,
                               strategies = v_scenarios)
icer_tab_15$subgrp <- "15%"

icer_tab_20 <- calculate_icers(cost = df_inc_shift_20$tot_cost_societal,
                               effect = df_inc_shift_20$tot_qalys,
                               strategies = v_scenarios)
icer_tab_20$subgrp <- "20%"

icer_tab_shift_health_mod <- rbind(icer_tab_5, icer_tab_10,
                               icer_tab_15, icer_tab_20)
icer_tab_shift_health_mod$Strategy <- factor(icer_tab_shift_health_mod$Strategy,
                                         levels = v_scenarios)
icer_tab_shift_health_mod$subgrp <- factor(icer_tab_shift_health_mod$subgrp,
                                       levels = c("5%", "10%",
                                                  "15%", "20%"))
icer_tab_shift_health_mod$Perspective <- "Modified Healthcare Sector"
icer_tab_shift_all <- rbind(icer_tab_shift_health, icer_tab_shift_health_mod)

# fig5b <- plot_icers_mod_subgrp(icer_tab_shift, label = "none") +
#   scale_y_continuous(labels = scales::comma, n.breaks = 7) +
#   coord_flip() + theme_bw() + ylab("Inc. Cost ($)") +
#   xlab("Inc. Effectiveness (QALYs)")
#
# fig_ice_health <- fig5a +
#   labs(title = "Healthcare Sector Perspective", tag = "A") +
#   guides(shape = guide_legend("Strategy",
#                               theme = theme(legend.title = element_text(vjust = 0.85))),
#          color = guide_legend(str_wrap("Percent of Lower Quality Donor Kidneys", 25), ncol = 2,
#                               theme = theme(legend.title = element_text(vjust = 0.85))),
#          linetype = guide_legend(element_blank()))
#
# fig_ice_soc    <- fig5b +
#   labs(title = "Modified Healthcare Sector Perspective", tag = "B")
#   guides(shape = guide_legend("Strategy",
#                               theme = theme(legend.title = element_text(vjust = 0.85))),
#          color = guide_legend(str_wrap("Percent of Lower Quality Donor Kidneys", 25), ncol = 2,
#                               theme = theme(legend.title = element_text (vjust = 0.85))),
#          linetype = guide_legend(element_blank()))
#
# ice_combined_shift <- ggarrange(fig_ice_health, fig_ice_soc, nrow = 1,
#                                legend = "bottom", common.legend = TRUE)
# ice_combined_shift

dat_text <- data.frame(
  subgrp = c("5%", "10%", "10%", "15%", "20%",
             "5%", "10%", "15%", "20%"),
  Perspective = c("Healthcare Sector", "Healthcare Sector", "Healthcare Sector",
                  "Healthcare Sector", "Healthcare Sector",
                  "Modified Healthcare Sector", "Modified Healthcare Sector",
                  "Modified Healthcare Sector", "Modified Healthcare Sector"),
  label = c(icer_health_5, icer_health_10_1, icer_health_10_2,
            icer_health_15, icer_health_20, "", "", "", ""),
  x = c(2000, 2000, 3600, 2000, 2000, 0, 0, 0, 0),
  y = c(0.3, 0.3, 0.51, 0.3, 0.3, 0, 0, 0, 0)
)

fig4_shift <- plot_icers_mod_subgrp(icer_tab_shift_all, label = "none") +
  geom_text(
    size    = 3,
    data    = dat_text,
    mapping = aes(x = x, y = y, label = label, angle = c(21.5, 21.5, 20, 22, 22,
                                                         0, 0, 0, 0))
  ) +
  guides(shape = guide_legend("Strategy"),
         color = guide_legend("KDPI Shift %"),
         linetype = guide_legend(element_blank())) +
  theme(strip.text.x = element_text(size = 12),
        strip.text.y = element_text(size = 12))


ggsave(filename = paste0(serv_path, "ice_combined_shift_", pdf_suffix),
       plot = fig4_shift, width = 10, height = 6)
ggsave(filename = paste0(serv_path, "ice_combined_shift_", png_suffix),
       plot = fig4_shift, width = 10, height = 6)


## Appendix Figures for Key Subgroups
# Generate df, using inc shift 5%
df_app_fig <- df_shift_05_all %>%
  dplyr::mutate(inmb_health_65 = (tot_qalys_65 * wtp) - tot_cost_health_65) %>%
  dplyr::mutate(inmb_soc_65 = (tot_qalys_65 * wtp) - tot_cost_societal_65) %>%
  dplyr::mutate(inmb_health_70 = (tot_qalys_70 * wtp) - tot_cost_health_70) %>%
  dplyr::mutate(inmb_soc_70 = (tot_qalys_70 * wtp) - tot_cost_societal_70) %>%
  dplyr::mutate(inmb_health_75 = (tot_qalys_75 * wtp) - tot_cost_health_75) %>%
  dplyr::mutate(inmb_soc_75 = (tot_qalys_75 * wtp) - tot_cost_societal_75) %>%
  dplyr::mutate(inmb_health_diab_0 = (tot_qalys_diab_0 * wtp) - tot_cost_health_diab_0) %>%
  dplyr::mutate(inmb_soc_diab_0 = (tot_qalys_diab_0 * wtp) - tot_cost_societal_diab_0) %>%
  dplyr::mutate(inmb_health_diab_1 = (tot_qalys_diab_1 * wtp) - tot_cost_health_diab_1) %>%
  dplyr::mutate(inmb_soc_diab_1 = (tot_qalys_diab_1 * wtp) - tot_cost_societal_diab_1) %>%
  dplyr::mutate(inmb_health_white = (tot_qalys_white * wtp) - tot_cost_health_white) %>%
  dplyr::mutate(inmb_soc_white = (tot_qalys_white * wtp) - tot_cost_societal_white) %>%
  dplyr::mutate(inmb_health_black = (tot_qalys_black * wtp) - tot_cost_health_black) %>%
  dplyr::mutate(inmb_soc_black = (tot_qalys_black * wtp) - tot_cost_societal_black) %>%
  dplyr::mutate(inmb_health_hispanic = (tot_qalys_hispanic * wtp) - tot_cost_health_hispanic) %>%
  dplyr::mutate(inmb_soc_hispanic = (tot_qalys_hispanic * wtp) - tot_cost_societal_hispanic) %>%
  dplyr::mutate(inmb_health_other = (tot_qalys_other * wtp) - tot_cost_health_other) %>%
  dplyr::mutate(inmb_soc_other = (tot_qalys_other * wtp) - tot_cost_societal_other) %>%
  gather(inmb, value, inmb_health_65:inmb_soc_other) %>%
  dplyr::select(c(scenario, inmb, value)) %>%
  dplyr::mutate(subgrp = case_when(
    str_detect(inmb, "5") ~ "age",
    str_detect(inmb, "0") ~ "age",
    str_detect(inmb, "diab") ~ "diab",
    str_detect(inmb, "white") ~ "race",
    str_detect(inmb, "black") ~ "race",
    str_detect(inmb, "hispanic") ~ "race",
    str_detect(inmb, "other") ~ "race",
  )) %>%
  dplyr::mutate(subgrp_spec = case_when(
    str_detect(inmb, "65") ~ "65-69",
    str_detect(inmb, "70") ~ "70",
    str_detect(inmb, "75") ~ "75+",
    str_detect(inmb, "diab_0") ~ "No",
    str_detect(inmb, "diab_1") ~ "Yes",
    str_detect(inmb, "white") ~ "NH White",
    str_detect(inmb, "black") ~ "NH Black",
    str_detect(inmb, "hispanic") ~ "Hispanic",
    str_detect(inmb, "other") ~ "NH Other",
  )) %>%
  dplyr::mutate(pers = case_when(
    str_detect(inmb, "health") ~ "Healthcare Sector",
    str_detect(inmb, "soc") ~ "Societal",
  )) %>%
  filter(scenario == "25%")
l_df_app_fig <- split(df_app_fig, df_app_fig$subgrp)



## Table 2: Event Counts
l_df_table2 <- list(base_df_pp, df_05, df_10, df_15, df_20, df_25)
m_counts <- matrix(nrow = 8, ncol = 18)
colnames(m_counts) <- c("base_mean", "base_min", "base_max",
                        "5_mean",  "5_min", "5_max",
                        "10_mean", "10_min", "10_max",
                        "15_mean", "15_min", "15_max",
                        "20_mean", "20_min", "20_max",
                        "25_mean", "25_min", "25_max")
x <- 1
for (df in l_df_table2) {
  y <- 0
  for (i in 32:37) {
    y <- y + 1
    m_counts[y, x]     <- round(mean(df[, i]), 0)
    m_counts[y, x + 1] <- round(min(df[, i]), 0)
    m_counts[y, x + 2] <- round(max(df[, i]), 0)
  }
  for (i in 92:93) {
    y <- y + 1
    m_counts[y, x]     <- round(mean(df[, i]), 0)
    m_counts[y, x + 1] <- round(min(df[, i]), 0)
    m_counts[y, x + 2] <- round(max(df[, i]), 0)
  }
  x <- x + 3
}

