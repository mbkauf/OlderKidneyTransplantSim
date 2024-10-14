library(readxl)
# library(haven)
library(data.table)

##### file paths #####
# in_path  <- "//tsclient/Documents/Simulation Model/Model Inputs/"
# out_path <- "//tsclient/Documents/Simulation Model/Model Outputs/"

in_path  <- "~/Documents/R/OlderKidneyTransplantSim/data/Model Inputs/"
# out_path <- "~/Documents/Simulation Model/Model Outputs/"

# if (Sys.getenv("LMOD_SYSHOST") == "sherlock") {
#   in_path  <- "/home/users/mbkauf/OlderKidneyTransplantSim/data/Model Inputs/"
# }

# in_path  <- "/home/users/mbkauf/OlderKidneyTransplantSim/data/Model Inputs/"

##### Data set #####
sample_df  <- as.matrix(fread(file = paste0(in_path, "sample_df_2.csv")))[, -1]

## Discount Rate
disc <- ((1 + 0.03)^(1 / 12)) - 1  # Monthly Discount Rate (3% Annual)

## Model Coefficients
ddtx_coef_trim       <- as.matrix(fread(paste0(in_path, "ddtx_coef_trim.csv")))[, -1]
ldtx_coef_trim       <- as.matrix(fread(paste0(in_path, "ldtx_coef_trim.csv")))[, -1]
mort_coef_trim       <- as.matrix(fread(paste0(in_path, "mort_coef_trim.csv")))[, -1]
remove_coef_trim     <- as.matrix(fread(paste0(in_path, "remove_coef_trim.csv")))[, -1]
ddtx_shape_trim      <- as.matrix(fread(paste0(in_path, "ddtx_shape_trim.csv")))[, -1]
ldtx_shape_trim      <- as.matrix(fread(paste0(in_path, "ldtx_shape_trim.csv")))[, -1]
mort_shape_trim      <- as.matrix(fread(paste0(in_path, "mort_shape_trim.csv")))[, -1]
remove_shape_trim    <- as.matrix(fread(paste0(in_path, "remove_shape_trim.csv")))[, -1]
mlogit_coef_trim     <- as.matrix(fread(paste0(in_path, "mlogit_coef_trim.csv")))
gs_mort_coef_trim    <- as.matrix(fread(paste0(in_path, "gs_mort_coef_trim.csv")))[, -19]
gl_coef_trim         <- as.matrix(fread(paste0(in_path, "gl_coef_trim.csv")))[, -20]
dial_mort_coef_trim  <- as.matrix(fread(paste0(in_path, "dial_mort_coef_trim.csv")))[, -21]
gl_mort_coef_trim    <- as.matrix(fread(paste0(in_path, "gl_mort_coef_trim.csv")))[, -21]
gs_mort_shape_trim   <- as.matrix(fread(paste0(in_path, "gs_mort_shape_trim.csv")))[, -1]
gl_shape_trim        <- as.matrix(fread(paste0(in_path, "gl_shape_trim.csv")))[, -1]
dial_mort_shape_trim <- as.matrix(fread(paste0(in_path, "dial_mort_shape_trim.csv")))[, -1]

generate_sse <- function(beta, mu_target, sd_target) {
  alpha <- mu_target * beta
  v_samples <- rgamma(n = 100000, shape = alpha, rate = beta)
  mean_samples <- mean(v_samples)
  sd_samples   <- sd(v_samples)

  out <- (mu_target - mean_samples)^2 + (sd_target - sd_samples)^2

  return(out)
}

optim_cost_params <- function(mu, sd, n_samples = 10000) {
  optim_result <- optim(par = 0.5, fn = generate_sse,
                        mu_target = mu,
                        sd_target = sd,
                        method = "L-BFGS-B",
                        lower = 0.0001,
                        upper = 1)
  beta  <- optim_result$par
  alpha <- mu * beta

  return(c(alpha = alpha, beta = beta))
}

get_psa_params_costs <- function(n_samples = 10000, params) {
  ## Load parameters
  mu_dial_65        <- params$v_c_dialysis[1]
  mu_dial_70        <- params$v_c_dialysis[2]
  mu_dial_75        <- params$v_c_dialysis[3]
  mu_dial_80        <- params$v_c_dialysis[4]
  mu_dial_85        <- params$v_c_dialysis[5]
  sd_dial           <- params$c_dialysis_sd
  mu_posttx_65      <- params$v_c_posttx[1]
  mu_posttx_70      <- params$v_c_posttx[2]
  mu_posttx_75      <- params$v_c_posttx[3]
  mu_posttx_80      <- params$v_c_posttx[4]
  mu_posttx_85      <- params$v_c_posttx[5]
  sd_posttx         <- params$c_posttx_sd
  mu_ddtx           <- params$c_ddtx
  sd_ddtx           <- params$c_ddtx_sd
  mu_oacc           <- params$c_oacc
  sd_oacc           <- params$c_oacc_sd
  mu_ldtx           <- params$c_ldtx
  sd_ldtx           <- params$c_ldtx_sd
  mu_dgf            <- params$c_dgf
  sd_dgf            <- params$c_dgf_sd
  mu_graftloss      <- params$c_graftloss
  sd_graftloss      <- params$c_graftloss_sd
  mu_giver_dialysis <- params$c_giver_dialysis
  sd_giver_dialysis <- params$c_giver_dialysis_sd
  mu_giver_posttx   <- params$c_giver_posttx
  sd_giver_posttx   <- params$c_giver_posttx_sd

  ## Dialysis costs
  # 65-69
  c_dial_65_params <- optim_cost_params(mu = mu_dial_65, sd = sd_dial)
  # 70-74
  c_dial_70_params <- optim_cost_params(mu = mu_dial_70, sd = sd_dial)
  # 75-79
  c_dial_75_params <- optim_cost_params(mu = mu_dial_75, sd = sd_dial)
  # 80-84
  c_dial_80_params <- optim_cost_params(mu = mu_dial_80, sd = sd_dial)
  # 85+
  c_dial_85_params <- optim_cost_params(mu = mu_dial_85, sd = sd_dial)

  ## Post-transplant costs
  # 65-69
  c_posttx_65_params <- optim_cost_params(mu = mu_posttx_65, sd = sd_posttx)
  # 70-74
  c_posttx_70_params <- optim_cost_params(mu = mu_posttx_70, sd = sd_posttx)
  # 75-79
  c_posttx_75_params <- optim_cost_params(mu = mu_posttx_75, sd = sd_posttx)
  # 80-84
  c_posttx_80_params <- optim_cost_params(mu = mu_posttx_80, sd = sd_posttx)
  # 85+
  c_posttx_85_params <- optim_cost_params(mu = mu_posttx_85, sd = sd_posttx)

  ## Deceased donor tx cost
  c_ddtx_params <- optim_cost_params(mu = mu_ddtx, sd = sd_ddtx)
  ## OACC costs
  c_oacc_params <- optim_cost_params(mu = mu_oacc, sd = sd_oacc)
  ## Living donor tx cost
  c_ldtx_params <- optim_cost_params(mu = mu_ldtx, sd = sd_ldtx)
  ## Delayed graft function cost
  c_dgf_params  <- optim_cost_params(mu = mu_dgf, sd = sd_dgf)
  ## Graft loss costs
  c_graftloss_params <- optim_cost_params(mu = mu_graftloss, sd = sd_graftloss)
  ## Caregiver time dialysis
  c_giver_dialysis_params <- optim_cost_params(mu = mu_giver_dialysis,
                                               sd = sd_giver_dialysis)
  ## Caregiver time post-tx
  c_giver_posttx_params <- optim_cost_params(mu = mu_giver_posttx,
                                             sd = sd_giver_posttx)

  return(list(c_dial_65_params = c_dial_65_params,
              c_dial_70_params = c_dial_70_params,
              c_dial_75_params = c_dial_75_params,
              c_dial_80_params = c_dial_80_params,
              c_dial_85_params = c_dial_85_params,
              c_posttx_65_params = c_posttx_65_params,
              c_posttx_70_params = c_posttx_70_params,
              c_posttx_75_params = c_posttx_75_params,
              c_posttx_80_params = c_posttx_80_params,
              c_posttx_85_params = c_posttx_85_params,
              c_ddtx_params = c_ddtx_params,
              c_oacc_params = c_oacc_params,
              c_ldtx_params = c_ldtx_params,
              c_dgf_params = c_dgf_params,
              c_graftloss_params = c_graftloss_params,
              c_giver_dialysis_params = c_giver_dialysis_params,
              c_giver_posttx_params = c_giver_posttx_params))
}

### Cost and QALY inputs
l_cost_params <- list(
  # monthly dialysis cost (65-69, 70-74, 75-79, 80-84, 85+)
  v_c_dialysis = c(8791, 8675, 8553, 8672, 8387),
  c_dialysis_sd = 153,
  # monthly post-tx cost (65-69, 70-74, 75-79, 80-84, 85+)
  v_c_posttx = c(3520, 3498, 3604, 3419, 3182),
  c_posttx_sd = 161,
  c_ddtx = 114688, # one-time deceased donor transplant cost (was 83536)
  c_ddtx_sd = 3475, # (was 2531)
  c_oacc = 112318, # one-time OACC charges
  c_oacc_sd = 24148,
  c_ldtx = 111976, # one-time living donor transplant cost (was 70384)
  c_ldtx_sd = 3393, # (was 2133)
  c_dgf = 41592, # cost of experiencing delayed graft function
  c_dgf_sd = 12440,
  c_graftloss = 85989, # one-time graft loss cost
  c_graftloss_sd = 19257,
  c_gsdeath = 0,  # one-time death with function cost (place-holder)
  c_gldeath = 0,  # one-time death with function cost (place-holder)
  c_wldeath = 0,  # one-time death on wait list (place-holder)
  # 20% increase in spending compared to those still on the list
  rr_removal = 1.2,
  c_giver_dialysis = 5867,  # monthly caregiver time costs, dialysis
  c_giver_dialysis_sd = 4444,
  c_giver_posttx = 1141,  # monthly caregiver time costs, post-transplant
  c_giver_posttx_sd = 859,
  c_patient_dialysis = 1611,  # monthly patient time costs, dialysis
  # monthly patient time costs, post-transplant
  v_c_patient_posttx = c(372, 186, 62),
  v_c_giver_posttx = c(1355, 677, 276)
)

l_cost_dist_params <- get_psa_params_costs(params = l_cost_params)
l_cost_params <- c(l_cost_params, l_cost_dist_params)

l_qaly_params <- list(
  ## QALY inputs (deterministic)
  # QALYs for males on dialysis by age (65-69, 70-79, 80+)
  v_u_dialysis_m = c(0.46, 0.42, 0.40),
  # QALYs for females on dialysis by age (65-69, 70-79, 80+)
  v_u_dialysis_f = c(0.43, 0.39, 0.34),
  # QALYs for males 0-3 months post-tx by age (65-69, 70-79, 80+)
  v_u_posttx_03_m = c(0.58, 0.54, 0.52),
  # QALYs for females 0-3 months post-tx by age (65-69, 70-79, 80+)
  v_u_posttx_03_f = c(0.55, 0.51, 0.46),
  # QALYs for males 4-8 months post-tx by age (65-69, 70-79, 80+)
  v_u_posttx_48_m = c(0.66, 0.63, 0.61),
  # QALYs for females 4-8 months post-tx by age (65-69, 70-79, 80+)
  v_u_posttx_48_f = c(0.64, 0.60, 0.55),
  # QALYs for males 9-12 months post-tx by age (65-69, 70-79, 80+)
  v_u_posttx_912_m = c(0.60, 0.56, 0.54),
  # QALYs for females 9-12 months post-tx by age (65-69, 70-79, 80+)
  v_u_posttx_912_f = c(0.57, 0.53, 0.48),
  # QALYs for males 13+ months post-tx by age (65-69, 70-79, 80+)
  v_u_posttx_13_m = c(0.59, 0.55, 0.53),
  # QALYs for females 13+ months post-tx by age (65-69, 70-79, 80+)
  v_u_posttx_13_f = c(0.56, 0.52, 0.47),
  u_gl = 0.20,  # disutility of graft loss event
  u_dgf = 0.20,  # disutility of delayed graft function
  u_death = 0,  # utility after death
  n_qaly = 727
)

l_qaly_dist_params <- list(
  ## QALY inputs (probabilistic)
  # alpha parameter for males on dialysis by age (65-69, 70-79, 80+)
  v_u_dialysis_m_alpha = l_qaly_params$v_u_dialysis_m * l_qaly_params$n_qaly,
  # beta parameter for males on dialysis by age (65-69, 70-79, 80+)
  v_u_dialysis_m_beta  = (1 - l_qaly_params$v_u_dialysis_m) * l_qaly_params$n_qaly,
  # alpha parameter for females on dialysis by age (65-69, 70-79, 80+)
  v_u_dialysis_f_alpha = l_qaly_params$v_u_dialysis_f * l_qaly_params$n_qaly,
  # beta parameter for females on dialysis by age (65-69, 70-79, 80+)
  v_u_dialysis_f_beta  = (1 - l_qaly_params$v_u_dialysis_f) * l_qaly_params$n_qaly,
  # alpha parameter for males 0-3 months post-tx by age (65-69, 70-79, 80+)
  v_u_posttx_03_m_alpha = l_qaly_params$v_u_posttx_03_m * l_qaly_params$n_qaly,
  # beta parameter for males 0-3 months post-tx by age (65-69, 70-79, 80+)
  v_u_posttx_03_m_beta  = (1 - l_qaly_params$v_u_posttx_03_m) * l_qaly_params$n_qaly,
  # alpha parameter for females 0-3 months post-tx by age (65-69, 70-79, 80+)
  v_u_posttx_03_f_alpha = l_qaly_params$v_u_posttx_03_f * l_qaly_params$n_qaly,
  # beta parameter for females 0-3 months post-tx by age (65-69, 70-79, 80+)
  v_u_posttx_03_f_beta  = (1 - l_qaly_params$v_u_posttx_03_f) * l_qaly_params$n_qaly,
  # alpha parameter for males 4-8 months post-tx by age (65-69, 70-79, 80+)
  v_u_posttx_48_m_alpha = l_qaly_params$v_u_posttx_48_m * l_qaly_params$n_qaly,
  # beta parameter for males 4-8 months post-tx by age (65-69, 70-79, 80+)
  v_u_posttx_48_m_beta  = (1 - l_qaly_params$v_u_posttx_48_m) * l_qaly_params$n_qaly,
  # alpha parameter for females 4-8 months post-tx by age (65-69, 70-79, 80+)
  v_u_posttx_48_f_alpha = l_qaly_params$v_u_posttx_48_f * l_qaly_params$n_qaly,
  # beta parameter for females 4-8 months post-tx by age (65-69, 70-79, 80+)
  v_u_posttx_48_f_beta  = (1 - l_qaly_params$v_u_posttx_48_f) * l_qaly_params$n_qaly,
  # alpha parameter for males 9-12 months post-tx by age (65-69, 70-79, 80+)
  v_u_posttx_912_m_alpha = l_qaly_params$v_u_posttx_912_m * l_qaly_params$n_qaly,
  # beta parameter for males 9-12 months post-tx by age (65-69, 70-79, 80+)
  v_u_posttx_912_m_beta  = (1 - l_qaly_params$v_u_posttx_912_m) * l_qaly_params$n_qaly,
  # alpha parameter for females 9-12 months post-tx by age (65-69, 70-79, 80+)
  v_u_posttx_912_f_alpha = l_qaly_params$v_u_posttx_912_f * l_qaly_params$n_qaly,
  # beta parameter for females 9-12 months post-tx by age (65-69, 70-79, 80+)
  v_u_posttx_912_f_beta  = (1 - l_qaly_params$v_u_posttx_912_f) * l_qaly_params$n_qaly,
  # alpha parameter for males 13+ months post-tx by age (65-69, 70-79, 80+)
  v_u_posttx_13_m_alpha = l_qaly_params$v_u_posttx_13_m * l_qaly_params$n_qaly,
  # beta parameter for males 13+ months post-tx by age (65-69, 70-79, 80+)
  v_u_posttx_13_m_beta  = (1 - l_qaly_params$v_u_posttx_13_m) * l_qaly_params$n_qaly,
  # alpha parameter for females 13+ months post-tx by age (65-69, 70-79, 80+)
  v_u_posttx_13_f_alpha = l_qaly_params$v_u_posttx_13_f * l_qaly_params$n_qaly,
  # beta parameter for females 13+ months post-tx by age (65-69, 70-79, 80+)
  v_u_posttx_13_f_beta  = (1 - l_qaly_params$v_u_posttx_13_f) * l_qaly_params$n_qaly
)

l_qaly_params <- c(l_qaly_params, l_qaly_dist_params)

# score  <- read.csv(paste0(in_path, "df_100.csv"))[,32]

# #### Equation Coefficients #####
# ### Wait List Equation Variables ###
# # Deceased Donor Tx (AFT Weibull)
# ddtx_coef_w   <- as.matrix(read_excel(paste0(in_path, "list_results_input.xlsx"),
#                                       range = "B5:B30", col_names = FALSE))
# ddtx_p_w      <- read_excel(paste0(in_path, "list_results_input.xlsx"),
#                             range = "B31", col_names = FALSE)[[1]]
# ddtx_cov <- as.matrix(read_excel(paste0(in_path, "ddtx_results.xlsx"),
#                                  sheet = "cov", col_names = FALSE))
#
# # Living Donor Tx (PH Gompertz)
# ldtx_coef_g   <- as.matrix(read_excel(paste0(in_path, "list_results_input.xlsx"),
#                                       range = "C5:C30", col_names = FALSE))
# ldtx_gamma_g  <- read_excel(paste0(in_path, "list_results_input.xlsx"),
#                             range = "C31", col_names = FALSE)[[1]]
# ldtx_cov <- as.matrix(read_excel(paste0(in_path, "ldtx_results.xlsx"),
#                                  sheet = "cov", col_names = FALSE))
#
# # Wait list mortality (AFT Weibull)
# mort_coef_w   <- as.matrix(read_excel(paste0(in_path, "list_results_input.xlsx"),
#                                       range = "D5:D30", col_names = FALSE))
# mort_p_w  <- read_excel(paste0(in_path, "list_results_input.xlsx"),
#                         range = "D31", col_names = FALSE)[[1]]
# mort_cov <- as.matrix(read_excel(paste0(in_path, "mort_results.xlsx"),
#                                  sheet = "cov", col_names = FALSE))
#
# # Other wait list removals (AFT Weibull)
# remove_coef_w   <- as.matrix(read_excel(paste0(in_path, "list_results_input.xlsx"),
#                                         range = "E5:E30", col_names = FALSE))
# remove_p_w      <- read_excel(paste0(in_path, "list_results_input.xlsx"),
#                               range = "E31", col_names = FALSE)[[1]]
# remove_cov <- as.matrix(read_excel(paste0(in_path, "remove_results.xlsx"),
#                                    sheet = "cov", col_names = FALSE))
#
# # 30-Day Outcomes (w/o region)
# # Base Outcome: Graft Success
# gl30_coef_ml_noreg   <- as.matrix(read_excel(paste0(in_path, "tx_results_input_no_reg.xlsx"),
#                                              range = "C5:C22", col_names = FALSE))
# dgf_coef_ml_noreg    <- as.matrix(read_excel(paste0(in_path, "tx_results_input_no_reg.xlsx"),
#                                              range = "D5:D22", col_names = FALSE))
# mort30_coef_ml_noreg <- as.matrix(read_excel(paste0(in_path, "tx_results_input_no_reg.xlsx"),
#                                              range = "E5:E22", col_names = FALSE))
# mlogit_cov <- as.matrix(read_excel(paste0(in_path, "mlogit_results.xlsx"),
#                                    sheet = "cov", col_names = FALSE))
#
# # Post-Tx: Graft Loss (Gompertz) (w/o region)
# gl_coef_g_noreg <- as.matrix(read_excel(paste0(in_path, "gl_results_input_no_reg.xlsx"),
#                                         range = "B5:B23", col_names = FALSE))
# gl_gamma_g_noreg <- read_excel(paste0(in_path, "gl_results_input_no_reg.xlsx"),
#                                range = "B24", col_names = FALSE)[[1]]
# gl_cov <- as.matrix(read_excel(paste0(in_path, "gl_results.xlsx"),
#                                sheet = "cov", col_names = FALSE))
#
# # Post-Tx: Graft Loss (Weibull) (w/o region)
# gl_coef_w_noreg <- as.matrix(read_excel(paste0(in_path, "gl_results_input_no_reg.xlsx"),
#                                         range = "C5:C23", col_names = FALSE))
# gl_p_w_noreg <- read_excel(paste0(in_path, "gl_results_input_no_reg.xlsx"),
#                            range = "C24", col_names = FALSE)[[1]]
#
# # Post-Tx: Death with Function (Gompertz) (w/o region)
# gs_mort_coef_g_noreg <- as.matrix(read_excel(paste0(in_path, "gs_mort_results_input_no_reg.xlsx"),
#                                              range = "B5:B22", col_names = FALSE))
# gs_mort_gamma_g_noreg <- read_excel(paste0(in_path, "gs_mort_results_input_no_reg.xlsx"),
#                                     range = "B23", col_names = FALSE)[[1]]
# gs_mort_cov <- as.matrix(read_excel(paste0(in_path, "gs_mort_results.xlsx"),
#                                     sheet = "cov", col_names = FALSE))
#
# # w/o region
# # Post-Tx: Death after graft loss (put months2gl variable last) (Weibull)
# dial_mort_coef_w_noreg <- as.matrix(read_excel(paste0(in_path, "gl_mort_results_input_no_reg.xlsx"),
#                                                range = "B5:B24", col_names = FALSE))
# dial_mort_lnp_w_noreg <- read_excel(paste0(in_path, "gl_mort_results_input_no_reg.xlsx"),
#                                     range = "B25", col_names = FALSE)[[1]]
# dial_mort_cov <- as.matrix(read_excel(paste0(in_path, "gl_mort_results.xlsx"),
#                                       sheet = "cov", col_names = FALSE))
#
# ### Same day graft loss - mortality
# gl_mort_coef_logit_noreg <- as.matrix(read_excel(paste0(in_path,
#                                                         "gl_mort_results_logit_input_no_reg_rev.xlsx"),
#                                                  range = "B5:B24", col_names = FALSE))
# gl_mort_cov <- as.matrix(read_excel(paste0(in_path, "gl_mort_results_logit_rev.xlsx"),
#                                     sheet = "cov", col_names = FALSE))
#

# ## Function to load parameters
# load_cea_parameters <- function() {
#   v_cost_params <- list(
#     v_c_dialysis = c(8791, 8675, 8553, 8672, 8387), # monthly dialysis cost (65-69, 70-74, 75-79, 80-84, 85+)
#     v_c_posttx = c(3520, 3498, 3604, 3419, 3182), # monthly post-tx cost (65-69, 70-74, 75-79, 80-84, 85+)
#     c_ddtx = 83536, # one-time deceased donor transplant cost
#     c_oacc = 112318, # one-time OACC charges
#     c_ldtx = 70384, # one-time living donor transplant cost (place-holder)
#     c_dgf = 41592, # cost of experiencing delayed graft function
#     c_graftloss = 85989, # one-time graft loss cost (place-holder)
#     c_gsdeath = 0,  # one-time death with function cost (place-holder)
#     c_gldeath = 0,  # one-time death with function cost (place-holder)
#     c_wldeath = 0,  # one-time death on wait list (place-holder)
#     rr_removal = 1.2, # 20% increase in spending compared to those still on the list
#     c_giver_dialysis = 5867,  # monthly caregiver time costs, dialysis
#     c_giver_posttx = 1141,  # monthly caregiver time costs, post-transplant
#     c_patient_dialysis = 1611,  # monthly patient time costs, dialysis
#     v_c_patient_posttx = c(372, 186, 62),  # monthly patient time costs, post-transplant
#     v_c_dialysis_shape = c(8791, 8675, 8553, 8672, 8387) * 0.01, # shape parameter for monthly dialysis cost (65-69, 70-74, 75-79, 80-84, 85+)
#     v_c_dialysis_rate = rep(0.01, 4), # shape parameter for monthly dialysis cost (65-69, 70-74, 75-79, 80-84, 85+)
#     v_c_posttx_shape = c(3520, 3498, 3604, 3419, 3182) * 0.01, # shape parameter for monthly dialysis cost (65-69, 70-74, 75-79, 80-84, 85+)
#     v_c_posttx_rate = rep(0.01, 4), # shape parameter for monthly dialysis cost (65-69, 70-74, 75-79, 80-84, 85+)
#     c_ddtx_shape = 83536 * 0.001,
#     c_ddtx_rate = 0.001

#   )

#   v_qaly_params <- list(
#     v_u_dialysis_m = c(0.46, 0.42, 0.40), # QALYs for males on dialysis by age (65-69, 70-79, 80+)
#     v_u_dialysis_f = c(0.43, 0.39, 0.34), # QALYs for females on dialysis by age (65-69, 70-79, 80+)
#     v_u_posttx_03_m = c(0.58, 0.54, 0.52), # QALYs for males 0-3 months post-tx by age (65-69, 70-79, 80+)
#     v_u_posttx_03_f = c(0.55, 0.51, 0.46), # QALYs for females 0-3 months post-tx by age (65-69, 70-79, 80+)
#     v_u_posttx_48_m = c(0.66, 0.63, 0.61), # QALYs for males 4-8 months post-tx by age (65-69, 70-79, 80+)
#     v_u_posttx_48_f = c(0.64, 0.60, 0.55), # QALYs for females 4-8 months post-tx by age (65-69, 70-79, 80+)
#     v_u_posttx_912_m = c(0.60, 0.56, 0.54), # QALYs for males 9-12 months post-tx by age (65-69, 70-79, 80+)
#     v_u_posttx_912_f = c(0.57, 0.53, 0.48), # QALYs for females 9-12 months post-tx by age (65-69, 70-79, 80+)
#     v_u_posttx_13_m = c(0.59, 0.55, 0.53), # QALYs for males 13+ months post-tx by age (65-69, 70-79, 80+)
#     v_u_posttx_13_f = c(0.56, 0.52, 0.47), # QALYs for females 13+ months post-tx by age (65-69, 70-79, 80+)
#     u_gl = 0.20, # disutility of graft loss event
#     u_dgf = 0.20, # disutility of delayed graft function
#     u_death = 0,  # utility after death
#     v_u_dialysis_m_alpha = c(9.527, 8.562, 8.027), # alpha parameter for males on dialysis by age (65-69, 70-79, 80+)
#     v_u_dialysis_m_beta  = c(11.305, 11.858, 12.074), # beta parameter for males on dialysis by age (65-69, 70-79, 80+)
#     v_u_dialysis_f_alpha = c(8.798, 7.726, 6.409), # alpha parameter for females on dialysis by age (65-69, 70-79, 80+)
#     v_u_dialysis_f_beta  = c(11.743, 12.170, 12.368), # beta parameter for females on dialysis by age (65-69, 70-79, 80+)
#     v_u_posttx_03_m_alpha = c(54.481, 51.877, 50.216), # alpha parameter for males 0-3 months post-tx by age (65-69, 70-79, 80+)
#     v_u_posttx_03_m_beta  = c(39.551, 43.945, 46.096), # beta parameter for males 0-3 months post-tx by age (65-69, 70-79, 80+)
#     v_u_posttx_03_f_alpha = c(52.562, 49.226, 44.470), # alpha parameter for females 0-3 months post-tx by age (65-69, 70-79, 80+)
#     v_u_posttx_03_f_beta  = c(42.938, 47.222, 51.497), # beta parameter for females 0-3 months post-tx by age (65-69, 70-79, 80+)
#     v_u_posttx_48_m_alpha = c(28.671, 28.398, 28.048), # alpha parameter for males 4-8 months post-tx by age (65-69, 70-79, 80+)
#     v_u_posttx_48_m_beta  = c(14.501, 16.958, 18.228), # beta parameter for males 4-8 months post-tx by age (65-69, 70-79, 80+)
#     v_u_posttx_48_f_alpha = c(28.510, 27.800, 26.335), # alpha parameter for females 4-8 months post-tx by age (65-69, 70-79, 80+)
#     v_u_posttx_48_f_beta  = c(16.380, 18.914, 21.712), # beta parameter for females 4-8 months post-tx by age (65-69, 70-79, 80+)
#     v_u_posttx_912_m_alpha = c(40.144, 38.499, 37.405), # alpha parameter for males 9-12 months post-tx by age (65-69, 70-79, 80+)
#     v_u_posttx_912_m_beta  = c(27.297, 30.589, 32.219), # beta parameter for males 9-12 months post-tx by age (65-69, 70-79, 80+)
#     v_u_posttx_912_f_alpha = c(38.942, 36.744, 33.497), # alpha parameter for females 9-12 months post-tx by age (65-69, 70-79, 80+)
#     v_u_posttx_912_f_beta  = c(29.830, 33.079, 36.399), # beta parameter for females 9-12 months post-tx by age (65-69, 70-79, 80+)
#     v_u_posttx_13_m_alpha = c(65.472, 62.530, 60.625), # alpha parameter for males 13+ months post-tx by age (65-69, 70-79, 80+)
#     v_u_posttx_13_m_beta  = c(46.205, 51.522, 54.138), # beta parameter for males 13+ months post-tx by age (65-69, 70-79, 80+)
#     v_u_posttx_13_f_alpha = c(63.311, 59.483, 53.950), # alpha parameter for females 13+ months post-tx by age (65-69, 70-79, 80+)
#     v_u_posttx_13_f_beta  = c(50.300, 55.513, 60.773) # beta parameter for females 13+ months post-tx by age (65-69, 70-79, 80+)
#   )

#   return(list(v_cost_params = v_cost_params,
#               v_qaly_params = v_qaly_params))
# }
