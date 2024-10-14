library(fitdistrplus)
library(MASS)
library(copula)
library(psych)
library(dtplyr)
library(tidyr)
library(survival)
library(dplyr)

## Function to run waitlist simulation
list.simulation <- function(seed,
                            data = sample_df,
                            ddtx_coef, ddtx_p,
                            ldtx_coef, ldtx_gamma,
                            mort_coef, mort_shape,
                            remove_coef, remove_p,
                            ddtx_rate = 0) {
  ### Function Inputs
  # seed: used to replicate results
  # n: sample size used for simulation
  # t: number of months to run simulation over
  ### Create bootstrapped data set
  set.seed(seed)
  num_patients <- nrow(data)

  ### Set up wait list equations
  # Deceased Donor Tx (AFT Weibull)
  ddtx_lambda <- exp(-exp(ddtx_p) * data %*% ddtx_coef)

  # Living Donor Tx (PH Gompertz)
  ldtx_lambda <- exp(data %*% ldtx_coef)

  # Wait list mortality (AFT Weibull)
  mort_lambda <- exp(-exp(mort_shape) * data %*% mort_coef)

  # Other wait list removal (AFT Weibull)
  remove_lambda <- exp(-exp(remove_p) * data %*% remove_coef)

  ### Initialize vectors and matrices
  event_results   <- rep(0, num_patients)

  ddtx_results    <- rep(0, num_patients)

  ldtx_results    <- rep(0, num_patients)

  mort_results    <- rep(0, num_patients)

  remove_results  <- rep(0, num_patients)

  m2event         <- rep(0, num_patients)
  events          <- rep(0, num_patients)

  m2age100        <- (rep(35, num_patients) - data[,1]) * 12

  c.mat           <- matrix(data = 0, nrow = num_patients, ncol = 420) # cost matrix
  u.mat           <- matrix(data = 0, nrow = num_patients, ncol = 420) # QALY matrix

  # Generate vector of months until candidate is 65
  # m2_65           <- boot_data[,1]
  # m2_65           <- replace(m2_65, m2_65>0, 0)
  # u65             <- replace(m2_65, m2_65<0, 1)
  # m2_65           <- (m2_65+0.5)*-12
  # m2_65           <- m2_65 * u65

  ### Run Simulation
  i <- 0
  while (all(event_results == 1) == FALSE) {
    i = i + 1
    ### Deceased Donor Tx?
    # Generate vector of ddtx probabilities
    ddtx_haz  <- (exp(ddtx_p) * ddtx_lambda * i^(exp(ddtx_p) - 1)) *
                   (1 + ddtx_rate)
    ddtx_prob <- 1 - exp(-ddtx_haz)

    # Determine who has event by using random uniform draws and determine who
    # had a new event this period
    # Vector of random probabilities from uniform
    rand_ddtx_surv <- runif(num_patients, min = 0, max = 1)
    # vector of those who had event in period i
    y              <- as.integer(rand_ddtx_surv < ddtx_prob)
    # vector of new events
    new_ddtx       <- as.integer(y == 1 & pmax(y, ddtx_results) > ddtx_results)

    # Check check that they are still at risk for event
    ddtx_results <- pmax(y, ddtx_results)
    match        <- as.integer(ddtx_results == pmax(ldtx_results, mort_results,
                                                    remove_results) &
                                 pmax(ldtx_results, mort_results, remove_results)==1)
    ddtx_results <- ddtx_results - match  # remove those who were not at risk for ddtx
    new_ddtx     <- new_ddtx - match  # final vector of new ddtx events

    ### Living Donor Tx?
    # Generate vector of ldtx probabilities
    ldtx_prob <- 1 - exp(-(ldtx_lambda * exp(ldtx_gamma*i)))  # probability vector

    # Determine who has event by using random uniform draws and determine who
    # had a new event this period
    rand_ldtx_surv <- runif(num_patients, min = 0, max = 1)  # random number vector
    x              <- as.integer(rand_ldtx_surv < ldtx_prob)  # vector of possible events
    new_ldtx       <- as.integer(x==1 & pmax(x, ldtx_results) > ldtx_results)

    # Check check that they are still at risk for event
    ldtx_results <- pmax(x, ldtx_results)
    match        <- as.integer(ldtx_results == pmax(ddtx_results, mort_results,
                                                    remove_results) &
                                 pmax(ddtx_results, mort_results, remove_results) == 1)
    ldtx_results <- ldtx_results - match
    new_ldtx     <- new_ldtx - match

    ### Wait list mortality?
    mort_prob <- 1 - exp(-(exp(mort_shape) * mort_lambda * i^(exp(mort_shape) - 1)))

    # Determine who has event by using random uniform draws and determine who
    # had a new event this period
    rand_mort_surv <- runif(num_patients, min = 0, max = 1)  # random number vector
    x              <- as.integer(rand_mort_surv < mort_prob)  # vector of possible events
    age100         <- as.integer(i >= m2age100)      # force death at age 100
    x              <- pmax(x, age100)   # new vector of possible events
    new_mort       <- as.integer(x == 1 & pmax(x, mort_results) > mort_results)

    # Check check that they are still at risk for event
    mort_results <- pmax(x, mort_results)
    match        <- as.integer(mort_results == pmax(ddtx_results, ldtx_results,
                                                    remove_results) &
                                 pmax(ddtx_results, ldtx_results, remove_results)==1)
    mort_results <- mort_results - match
    new_mort     <- new_mort - match

    ### Other List Removal?
    # Generate vector of list removal probabilities
    remove_prob <- 1 - exp(-(exp(remove_p) * remove_lambda * i^(exp(remove_p) - 1)))

    # Determine who has event by using random uniform draws and determine who
    # had a new event this period
    rand_remove_surv <- runif(num_patients, min = 0, max = 1)  # Vector of random probabilities from uniform
    y                <- as.integer(rand_remove_surv < remove_prob)  # vector of those who had event in period i
    new_remove       <- as.integer(y==1 & pmax(y, remove_results) > remove_results)

    # Check check that they are still at risk for event
    remove_results <- pmax(y, remove_results)
    match          <- as.integer(remove_results == pmax(ddtx_results,
                                                        ldtx_results, mort_results)
                                 & pmax(ddtx_results, ldtx_results, mort_results)==1)
    remove_results <- remove_results - match
    new_remove     <- new_remove - match

    # Create vector the tracks cumulative removal events, create patient trace
    # remove_trace[,i] <- remove_results

    ### Time to any event
    event_results <- pmax(ddtx_results, ldtx_results,
                          mort_results, remove_results)
    m2event       <- m2event + pmax(new_ddtx, new_ldtx,
                                    new_mort, new_remove)*i

    ### Update cost and QALY matrix
    # c.mat[,i][event_results==0] <- c_dialysis
    # c.mat[,i][new_ddtx==1]      <- c_ddtx
    # c.mat[,i][new_ldtx==1]      <- c_ldtx
    # c.mat[,i][new_mort==1]      <- c_wldeath
    # c.mat[,i][new_remove==1]    <- c_remove
    #
    # u.mat[,i][event_results==0] <- u_dialysis
  }

  ### Build dataset for analysis
  # results_sim <- as.data.frame(cbind(data, event_results, ddtx_results,
  #                                    ldtx_results, mort_results, remove_results,
  #                                    m2event))
  results_sim <- cbind(data, event_results, ddtx_results,
                       ldtx_results, mort_results, remove_results,
                       m2event)

  results_sim <- lazy_dt(results_sim) %>%
    rename(
      event = event_results,
      ddtx = ddtx_results,
      ldtx = ldtx_results,
      mort = mort_results,
      remove = remove_results,
      months_to_event = m2event
    ) %>%
    dplyr::mutate(months_to_event = ifelse(event == 0, i,
                                           months_to_event)) %>%
    dplyr::mutate(race = ifelse(black == 0 & hispanic == 0 &
                                  other_race == 0, 0,
                                ifelse(black == 1 & hispanic == 0 &
                                         other_race == 0, 1,
                                ifelse(black==0 & hispanic==1 & other_race==0, 2, 3)))) %>%
    mutate(race = factor(race, levels = c(0, 1, 2, 3),
                         labels = c("White", "Black",
                                    "Hispanic", "Other"))) %>%
    as_tibble()

  return(results_sim)
}

## Function to run post-tx simulation
run.simulation <- function(seed, data = list_sim_df,
                           gl30_coef,
                           dgf_coef,
                           mort30_coef,
                           gl_coef,
                           gl_shape,
                           gs_mort_coef,
                           gs_mort_shape,
                           dial_mort_coef,
                           dial_mort_shape,
                           gl_mort_coef,
                           france = FALSE,
                           kdpi_pct_shift = 0,
                           kdpi_max_shift = 10,
                           gl_mort_reduction = 1) {
  ### Function Inputs
  # seed: used to replicate results
  # n: sample size used for simulation
  # t: number of months to run simulation over
  ### Create bootstrapped data set
  set.seed(seed)
  # months <- t

  ### Start first 30 days
  ### Make the following changes to the data used for the next part
  # keep only deceased donor transplants
  # update to age at transplant
  # update to time on dialysis prior to transplant
  # update list year to transplant year
  # draw from joint distribution of cold ischemia time and KDPI
  # reorder variables to match order for transplant outcomes
  '"
  data <- cbind(data, Z)
  names(data)[ncol(data)] <- "kdpi"
  names(data)[ncol(data)-1] <- "rec_cold_isch_tm"

  ddtx            <- data$ddtx
  ldtx            <- data$ldtx
  mort            <- data$mort
  remove          <- data$remove
  months_to_event <- data$months_to_event

  optn_reg2  <- data$optn_reg2
  optn_reg3  <- data$optn_reg3
  optn_reg4  <- data$optn_reg4
  optn_reg5  <- data$optn_reg5
  optn_reg6  <- data$optn_reg6
  optn_reg7  <- data$optn_reg7
  optn_reg8  <- data$optn_reg8
  optn_reg9  <- data$optn_reg9
  optn_reg10 <- data$optn_reg10
  optn_reg11 <- data$optn_reg11

  data <- data %>%
    mutate(rec_age_at_tx_c = can_age_at_listing_c + (months_to_event/12)) %>%
    mutate(years_dial = baseline_yrs_dial + (months_to_event/12)) %>%
    mutate(tx_year_c = list_year_c + (months_to_event/12)) %>%
    mutate(tx_year_c = ifelse(tx_year_c > 10, 10, tx_year_c)) %>%
    mutate(constant = 1) %>%
    select(rec_age_at_tx_c,	male,	black, hispanic, other_race,
          years_dial, canhx_cpra, kdpi, can_blood_ab, can_blood_b,
          can_blood_o, rec_cold_isch_tm, tx_year_c, diab_stat,
          copd, pvd,	ang_cad, constant)


  data <- as.matrix(data)
  num_patients <- nrow(data)

  "'
  no_ddtx <- data %>%
    dplyr::filter(ddtx==0)

  posttx_df <- data %>%
    dplyr::filter(ddtx==1)

  colvars = names(posttx_df)
  colvars = colnames(posttx_df)
  start_loc = match("optn_reg2",colvars)
  end_loc = match("optn_reg11",colvars)
  optn <- posttx_df[,start_loc:end_loc]

  if(france == FALSE){
    invisible(Z <- rand_kdpi_cold(num = nrow(posttx_df),
                                  pct_shift = kdpi_pct_shift,
                                  max_shift = kdpi_max_shift))
  } else {
    invisible(Z <- rand_kdpi_cold(num = nrow(posttx_df), france = TRUE,
                                  pct_shift = kdpi_pct_shift,
                                  max_shift = kdpi_max_shift))
  }

  # posttx_df <- cbind(posttx_df, Z)
  # names(posttx_df)[ncol(posttx_df)] <- "kdpi"
  # names(posttx_df)[ncol(posttx_df)-1] <- "rec_cold_isch_tm"

  posttx_df <- cbind(posttx_df, Z)
  names(posttx_df)[ncol(posttx_df)] <- "kdpi"
  names(posttx_df)[ncol(posttx_df)-1] <- "rec_cold_isch_tm"

  # data <- as.data.frame(data)

  ddtx_df <- posttx_df

  posttx_df <- lazy_dt(posttx_df) %>%
    dplyr::mutate(rec_age_at_tx_c = can_age_at_listing_c + (months_to_event/12)) %>%
    dplyr::mutate(years_dial = baseline_yrs_dial + (months_to_event/12)) %>%
    dplyr::mutate(tx_year_c = list_year_c + (months_to_event/12)) %>%
    dplyr::mutate(tx_year_c = ifelse(tx_year_c > 10, 10, tx_year_c)) %>%
    dplyr::mutate(constant = 1) %>%
    dplyr::select(rec_age_at_tx_c,	male,	black, hispanic, other_race,
                  years_dial, canhx_cpra, kdpi, can_blood_ab, can_blood_b,
                  can_blood_o, rec_cold_isch_tm, tx_year_c, diab_stat,
                  copd, pvd,	ang_cad, constant) %>%
    as_tibble()
  posttx_df <- as.matrix(posttx_df)
  num_patients <- nrow(posttx_df)

  ### Which 30-Day event? (Multinomial Logit)
  # Base Outcome: Graft Success
  gl30_exb   <- exp(posttx_df %*% gl30_coef)
  dgf_exb    <- exp(posttx_df %*% dgf_coef)
  mort30_exb <- exp(posttx_df %*% mort30_coef)

  gl30_prob     <- gl30_exb / (1 + gl30_exb + dgf_exb + mort30_exb)
  dgf_prob      <- dgf_exb / (1 + gl30_exb + dgf_exb + mort30_exb)
  mort30_prob   <- mort30_exb / (1 + gl30_exb + dgf_exb + mort30_exb)
  success_prob  <- 1 / (1 + gl30_exb + dgf_exb + mort30_exb)

  gl30_range    <- gl30_prob
  dgf_range     <- gl30_prob + dgf_prob
  mort30_range  <- dgf_prob + mort30_prob

  dgf_results     <- rep(0, num_patients)
  mort30_results  <- rep(0, num_patients)
  success_results <- rep(0, num_patients)
  tx_outcome      <- rep(0, num_patients)

  rand_num_v  <- runif(num_patients, min = 0, max = 1)  # Vector of random probabilities from uniform

  gl30_results    <- as.integer(rand_num_v <= gl30_range)
  gl_results      <- as.integer(rand_num_v <= gl30_range)
  dgf_results     <- as.integer((rand_num_v > gl30_range) &
                                  (rand_num_v <= dgf_range))
  mort30_results  <- as.integer((rand_num_v > dgf_range) &
                                  (rand_num_v <= mort30_range))
  success_results <- as.integer(rand_num_v > mort30_range)

  tx_outcome <- tx_outcome + gl30_results + (dgf_results*2) + (mort30_results*3)

  data_gl <- as.data.frame(cbind(posttx_df, dgf_results))


  data_gl_mort <- lazy_dt(data_gl) %>%
    mutate(m2gl = rep(0, nrow(data_gl))) %>%
    dplyr::select(rec_age_at_tx_c,	male,	black, hispanic, other_race,
           years_dial, canhx_cpra, kdpi, can_blood_ab, can_blood_b,
           can_blood_o, rec_cold_isch_tm, tx_year_c, diab_stat, copd,
           pvd, ang_cad, dgf_results, m2gl, constant) %>%
    as_tibble()

  data_gs_mort <- lazy_dt(data_gl) %>%
    dplyr::select(rec_age_at_tx_c,	male,	black, hispanic, other_race,
           years_dial, canhx_cpra, kdpi, can_blood_ab, can_blood_b,
           can_blood_o, rec_cold_isch_tm, diab_stat, copd,
           pvd, ang_cad, dgf_results, constant) %>%
    as_tibble()
  data_gl <- lazy_dt(data_gl) %>%
    dplyr::select(rec_age_at_tx_c,	male,	black, hispanic, other_race,
           years_dial, canhx_cpra, kdpi, can_blood_ab, can_blood_b,
           can_blood_o, rec_cold_isch_tm, tx_year_c, diab_stat, copd,
           pvd, ang_cad, dgf_results, constant) %>%
    as_tibble()

  data_gl_mort <- as.matrix(data_gl_mort)
  data_gs_mort <- as.matrix(data_gs_mort)
  data_gl <- as.matrix(data_gl)

  # build.data(n=num_patients)
  ### Graft loss (PH Gompertz)
  gl_lambda <- exp(data_gl %*% gl_coef)

  ### Death w/ functioning graft (PH Gompertz)
  gs_mort_lambda <- exp(data_gs_mort %*% gs_mort_coef)

  ### Initialize vectors and matrices
  gs_mort_results <- rep(0, num_patients)
  gl_results      <- rep(0, num_patients)
  gl_mort_results <- rep(0, num_patients)

  m2gl      <- rep(0, num_patients)
  m2gs_mort <- rep(0, num_patients)
  m2gl_mort <- rep(0, num_patients)

  death    <- rep(0, num_patients)
  m2age100 <- (rep(35, num_patients) - data_gs_mort[, 1]) * 12

  ### Run Simulation
  i <- 0
  while (all(death == 1) == FALSE) {
    i = i + 1
    ### First we will see who has died w/ functioning graft
    gs_mort_haz <- gs_mort_lambda * exp(gs_mort_shape*i) # vector of hazards
    gs_mort_prob <- 1 - exp(-gs_mort_haz)  # Graft loss probability vector
    rand_gs_mort_surv <- runif(num_patients, min = 0, max = 1)  # random number vector
    x <- as.integer(rand_gs_mort_surv < gs_mort_prob)  # vector of possible events
    age100         <- as.integer(i >= m2age100)      # force death at age 100
    x              <- pmax(x, age100)   # new vector of possible events
    new_gs_mort <- as.integer(x==1 & pmax(x, gs_mort_results) > gs_mort_results)
    gs_mort_results <- pmax(x, gs_mort_results)   # create vector that combines those with previous deaths to new deaths
    match <- as.integer(gs_mort_results == 1 &
                          (gl_results==1 | mort30_results==1 | gl30_results==1))
    gs_mort_results <- gs_mort_results - match
    new_gs_mort <- new_gs_mort - match
    m2gs_mort <- m2gs_mort + new_gs_mort*i
    # gs_mort_trace[,i] <- gs_mort_results

    ### Next we see who died after graft failure
    data_gl_mort[,ncol(data_gl_mort) - 1] <- m2gl # Update data set
    dial_mort_lambda <- exp(-exp(dial_mort_shape) * data_gl_mort %*% dial_mort_coef) # Update lambda vector
    dial_mort_haz <- exp(dial_mort_shape) * dial_mort_lambda * i^(exp(dial_mort_shape) - 1)  # Graft loss hazard vector
    dial_mort_prob <- 1 - exp(-dial_mort_haz)  # Graft loss probability vector
    rand_dial_surv <- runif(num_patients, min = 0, max = 1)  # Vector of random probabilities from uniform
    z <- as.integer(rand_dial_surv < dial_mort_prob)  # vector of those who had event in period i
    age100         <- as.integer(i >= m2age100)      # force death at age 100
    z              <- pmax(z, age100)   # new vector of possible events
    new_gl_mort <- as.integer(z==1 & pmax(z, gl_mort_results) > gl_mort_results)
    gl_mort_results <- pmax(z, gl_mort_results)  # create vector that combines those with previous graft loss to new graft loss
    gl_mort_results <- gl_results * gl_mort_results
    new_gl_mort <- gl_results * new_gl_mort
    m2gl_mort <- m2gl_mort + new_gl_mort*i

    gl_haz <- gl_lambda * exp(gl_shape*i) # vector of hazards (Gompertz)
    gl_prob <- 1 - exp(-gl_haz)  # Graft loss probability vector

    rand_gl_surv <- runif(num_patients, min = 0, max = 1)  # Vector of random probabilities from uniform
    y <- as.integer(rand_gl_surv < gl_prob)  # vector of those who had event in period i
    new_gl <- as.integer(y==1 & pmax(y, gl_results) > gl_results)
    gl_results <- pmax(y, gl_results)  # create vector that combines those with previous graft loss to new graft loss
    match <- as.integer(gl_results==1 &
                          (gs_mort_results==1 | mort30_results==1))
    gl_results <- gl_results - match
    new_gl <- new_gl - match
    m2gl <- m2gl + new_gl*i


    ### Day of graft loss death
    gl_mort_prob <- exp(data_gl_mort %*% gl_mort_coef) /
      (1 + exp(data_gl_mort %*% gl_mort_coef))
    gl_mort_prob <- gl_mort_prob * gl_mort_reduction
    rand_gl_mort_surv <- runif(num_patients, min = 0, max = 1)  # Vector of random probabilities from uniform
    a <- as.integer(rand_gl_mort_surv < gl_mort_prob)  # vector of those who had event in period i
    new_gl_mort <- as.integer(a==1 & pmax((new_gl *a), gl_mort_results) > gl_mort_results)
    gl_mort_results <- pmax(new_gl * a, gl_mort_results)
    m2gl_mort <- m2gl_mort + new_gl_mort*i

    # Is anyone still alive?
    death = pmax(gs_mort_results, gl_mort_results, age100)
  }

  ### Build dataset for analysis
  # list_vars_tx = c("graft_loss", "gs_death", "gl_death", "gl_time",
  #                  "gs_death_time", "gl_death_time", "dgf",
  #                  "kdpi", "rec_cold_isch_tm")
  #
  # sim_df <- as.data.frame(cbind(data_gl, ddtx, ldtx, mort, remove,
  #                                months_to_event, gl_results, gs_mort_results,
  #                                gl_mort_results, m2gl, m2gs_mort, m2gl_mort,
  #                                optn_reg2, optn_reg3, optn_reg4, optn_reg5,
  #                                optn_reg6, optn_reg7, optn_reg8, optn_reg9,
  #                                optn_reg10, optn_reg11, tx_outcome))
  #
  # sim_df <- sim_df %>%
  #   rename(
  #     graft_loss = gl_results,
  #     gs_death = gs_mort_results,
  #     gl_death = gl_mort_results,
  #     gl_time = m2gl,
  #     gs_death_time = m2gs_mort,
  #     gl_death_time = m2gl_mort,
  #     dgf = dgf_results
  #   ) %>%
  #   mutate_at(.vars = list_vars_tx, .funs = list(~ifelse(ddtx==1, ., NA))) %>%
  #   mutate(gl_time = ifelse(graft_loss==0 & gs_death==0,
  #                           i, gl_time)) %>%
  #   mutate(gl_time = ifelse(graft_loss==0 & gs_death==1,
  #                           gs_death_time, gl_time)) %>%
  #   mutate(gs_death_time = ifelse(gs_death==0 & graft_loss==0,
  #                                 i, gs_death_time)) %>%
  #   mutate(gs_death_time = ifelse(gs_death==0 & graft_loss==1,
  #                                 gl_time, gs_death_time)) %>%
  #   mutate(gl_death_time = gl_death_time - gl_time) %>%
  #   mutate(gl_death_time = ifelse(gl_death==0 & graft_loss==1,
  #                                 -gl_time+i, gl_death_time)) %>%
  #   mutate(gl_death_time = ifelse(gl_death_time < 0, NA, gl_death_time)) %>%
  #   mutate(race = ifelse(black==0 & hispanic==0 & other_race==0, 0,
  #                        ifelse(black==1 & hispanic==0 & other_race==0, 1,
  #                               ifelse(black==0 & hispanic==1 & other_race==0, 2, 3)))) %>%
  #   # mutate(sim = rep(1, num_patients)) %>%
  #   mutate(race = factor(race, levels = c(0, 1, 2, 3),
  #                        labels = c("White", "Black",
  #                                   "Hispanic", "Other"))) %>%
  #   #  mutate(sim = factor(sim, levels = c(0, 1),
  #   #                      labels = c("Observed", "Simulated"))) %>%
  #   mutate(tx_outcome = factor(tx_outcome, levels = c(0, 1, 2, 3),
  #                              labels = c("GS", "GL",
  #                                         "DGF", "Mort")))
  #
  # return(sim_df)

  list_vars_tx = c("gl_results", "gs_mort_results", "gl_mort_results", "m2gl",
                   "m2gs_mort", "m2gl_mort", "dgf_results", "tx_outcome",
                   "kdpi", "rec_cold_isch_tm")

  no_ddtx <- cbind(no_ddtx, setNames(lapply(list_vars_tx,
                                            function(x) x=NA), list_vars_tx))

  sim_df <- as.data.frame(cbind(ddtx_df, gl_results, gs_mort_results,
                                gl_mort_results, m2gl, m2gs_mort, m2gl_mort,
                                dgf_results, tx_outcome))
  all_df <- rbind(sim_df, no_ddtx)
  all_df <- lazy_dt(all_df) %>%
    rename(
      graft_loss = gl_results,
      gs_death = gs_mort_results,
      gl_death = gl_mort_results,
      gl_time = m2gl,
      gs_death_time = m2gs_mort,
      gl_death_time = m2gl_mort,
      dgf = dgf_results
    ) %>%
    mutate(gl_time = ifelse(graft_loss==0 & gs_death==0,
                            i, gl_time)) %>%
    mutate(gl_time = ifelse(graft_loss==0 & gs_death==1,
                            gs_death_time, gl_time)) %>%
    mutate(gs_death_time = ifelse(gs_death==0 & graft_loss==0,
                                  i, gs_death_time)) %>%
    mutate(gs_death_time = ifelse(gs_death==0 & graft_loss==1,
                                  gl_time, gs_death_time)) %>%
    mutate(gl_death_time = gl_death_time - gl_time) %>%
    mutate(gl_death_time = ifelse(gl_death==0 & graft_loss==1,
                                  -gl_time+i, gl_death_time)) %>%
    mutate(gl_death_time = ifelse(gl_death_time < 0, NaN, gl_death_time)) %>%
    mutate(race = ifelse(black==0 & hispanic==0 & other_race==0, 0,
                         ifelse(black==1 & hispanic==0 & other_race==0, 1,
                                ifelse(black==0 & hispanic==1 & other_race==0, 2, 3)))) %>%
    as_tibble()

  return(all_df)
}

## Function to generate random values for KDPI and cold ischemia time
rand_kdpi_cold <- function(num, france=FALSE, shift=TRUE,
                           pct_shift = 0, max_shift = 10) {
  if (france == TRUE) {
    kdpi_shape1 = 1.045481
    kdpi_shape2 = 0.697808
  } else {
    kdpi_shape1 = 1.458402
    kdpi_shape2 = 1.033602
  }

  cor <- 0.130092
  myCop <- copula::normalCopula(param=cor, dim = 2, dispstr = "un")
  myMvd <- copula::mvdc(copula=myCop, margins=c("gamma", "beta"),
                        paramMargins=list(list(shape = 3.949532,
                                               rate = 0.2128366),
                                          list(shape1 = kdpi_shape1,
                                               shape2 = kdpi_shape2)))
  Z <- copula::rMvdc(num, myMvd)
  colnames(Z) <- c("rec_cold_isch_tm", "kdpi")
  Z[,2] <- Z[,2] * 100
  Z[,2] <- round(Z[,2],0)
  # return(Z)

  # Test discrete/continuous copula
  # kdpi_counts <- c(4, 74, 88, 118, 111, 91, 102, 111, 117, 113, 112, 99, 112,
  #                  81, 103, 142, 106, 161, 155, 199, 200, 173, 222, 256, 208,
  #                  169, 168, 233, 234, 188, 177, 212, 173, 231, 177, 173, 229,
  #                  283, 231, 207, 266, 239, 245, 212, 256, 234, 212, 254, 279,
  #                  220, 184, 284, 225, 243, 228, 225, 264, 287, 196, 289, 356,
  #                  287, 336, 260, 237, 266, 285, 292, 277, 289, 341, 354, 316,
  #                  281, 348, 334, 286, 369, 315, 298, 316, 307, 349, 302, 336,
  #                  350, 306, 297, 355, 395, 344, 391, 300, 361, 326, 327, 300,
  #                  229, 270, 187, 131)
  # cum <- c(0,cumsum(kdpi_counts))
  # kdpi_counts <- c(kdpi_counts, 0)
  # m_kdpi_w <- cbind(seq(from = 0, to = 101), kdpi_counts, cum)
  # colnames(m_kdpi_w) <- c("kdpi", "w", "cum")
  # m_kdpi_w[,2:3] <- m_kdpi_w[,2:3] / sum(m_kdpi_w[,2])
  #
  # r     <- cor
  # sigma <- matrix(c(1, r, r, 1), ncol = 2)
  # s     <- chol(sigma)
  # z     <- s %*% matrix(rnorm(num*2), nrow=2)
  # u     <- pnorm(z)
  #
  # v_cold <- qgamma(p = u[1,], shape = 3.949532, rate = 0.2128366)
  # values <- u[2,]
  # v_kdpi <- m_kdpi_w[,1][findInterval(x = values, vec = m_kdpi_w[,3],
  #                                     all.inside = T)]
  #
  # test <- cbind(v_cold, v_kdpi)
  # colnames(test) <- c("rec_cold_isch_tm", "kdpi")

  if (shift == TRUE) {
    pct_0   <- 1 - pct_shift
    v_shift <- sample(x = seq(0, max_shift), size = num, replace = TRUE,
                      prob = c(pct_0, rep((pct_shift/max_shift), max_shift)))
    Z[,2] <- Z[,2] + v_shift
    Z[,2] <- sapply(Z[,2], FUN = function(x) min(x, 100))
  }

  return(Z)
}

## Function to set up a normal copula allowing for custom margins (continous or discrete)
custom_copula <- function(sigma, n.params, n.draws,
                          v_dist) {
  s <- chol(sigma)
  z <- s %*% matrix(rnorm(num*n.params), nrow = n.params)
  u <- pnorm(z)

  # for(i in 1:length(v_dist)){
  #   if(v_dist %in% c("gamma",  "beta", "normal")) {
  #     if(v_dist=="gamma")
  #     v_tmp
  #   }
  # }


}


## Function to run code with output suppressed
hush=function(code){
  sink("NUL") # use /dev/null in UNIX
  tmp = code
  sink()
  return(tmp)
}

## Function to vary parameters for survival equations
chol_coef_surv <- function(cov_mat, v_coef, shape) {
  chol_mat <- chol(cov_mat)
  v_in <- qnorm(runif(nrow(cov_mat), 0, 1), 0, 1)  # vector of random inverse normal
  v_tz <- rep(0, nrow(cov_mat))

  for (i in 1:nrow(cov_mat)) {
    v_tz[i] <- v_in %*% chol_mat[i,]
  }
  shape_tz <- v_tz[nrow(cov_mat)]
  coef_tz <- v_tz[1: nrow(cov_mat)-1]
  coef_psa <<- v_coef + coef_tz
  shape_psa <<- shape + shape_tz
}

## Function to vary parameters for non-survival equations
chol_coef_logit <- function(cov_mat, v_coef) {
  chol_mat <- chol(cov_mat)
  v_in <- qnorm(runif(nrow(cov_mat), 0, 1), 0, 1)  # vector of random inverse normal
  v_tz <- rep(0, nrow(cov_mat))

  for (i in 1:nrow(cov_mat)) {
    v_tz[i] <- v_in %*% chol_mat[i,]
  }
  coef_tz <- v_tz
  coef_psa <<- v_coef + coef_tz
}

## Function to get probabilistically varied coefficients for waitlist equations
psa_list <- function(ddtx_cov_mat = ddtx_cov, ddtx_v_coef = ddtx_coef_w,
                     ddtx_shape = ddtx_p_w, ldtx_cov_mat = ldtx_cov,
                     ldtx_v_coef = ldtx_coef_g, ldtx_shape = ldtx_gamma_g,
                     mort_cov_mat = mort_cov, mort_v_coef = mort_coef_w,
                     mort_shape = mort_p_w, remove_cov_mat = remove_cov,
                     remove_v_coef = remove_coef_w, remove_shape = remove_p_w) {
  chol_coef_surv(cov_mat = ddtx_cov_mat, v_coef = ddtx_v_coef, shape = ddtx_shape)
  ddtx_coef_psa <- coef_psa
  ddtx_shape_psa <- shape_psa

  chol_coef_surv(cov_mat = ldtx_cov_mat, v_coef = ldtx_v_coef, shape = ldtx_shape)
  ldtx_coef_psa <- coef_psa
  ldtx_shape_psa <- shape_psa

  chol_coef_surv(cov_mat = mort_cov_mat, v_coef = mort_v_coef, shape = mort_shape)
  mort_coef_psa <- coef_psa
  mort_shape_psa <- shape_psa

  chol_coef_surv(cov_mat = remove_cov_mat, v_coef = remove_v_coef, shape = remove_shape)
  remove_coef_psa <- coef_psa
  remove_shape_psa <- shape_psa

  return(list(ddtx_coef_psa = ddtx_coef_psa, ddtx_shape_psa = ddtx_shape_psa,
              ldtx_coef_psa = ldtx_coef_psa, ldtx_shape_psa = ldtx_shape_psa,
              mort_coef_psa = mort_coef_psa, mort_shape_psa = mort_shape_psa,
              remove_coef_psa = remove_coef_psa, remove_shape_psa = remove_shape_psa))
}

## Function to get probabilistically varied coefficients for post-tx equations
psa_posttx <- function(gl_cov_mat = gl_cov,
                       gl_v_coef = gl_coef_g_noreg,
                       gl_shape = gl_gamma_g_noreg,
                       gs_mort_cov_mat = gs_mort_cov,
                       gs_mort_v_coef = gs_mort_coef_g_noreg,
                       gs_mort_shape = gs_mort_gamma_g_noreg,
                       dial_mort_cov_mat = dial_mort_cov,
                       dial_mort_v_coef = dial_mort_coef_w_noreg,
                       dial_mort_shape = dial_mort_lnp_w_noreg,
                       gl_mort_cov_mat = gl_mort_cov,
                       gl_mort_v_coef = gl_mort_coef_logit_noreg,
                       mlogit_cov_mat = mlogit_cov,
                       mlogit_v_coef = rbind(gl30_coef_ml_noreg,
                                             dgf_coef_ml_noreg,
                                             mort30_coef_ml_noreg)) {
  chol_coef_surv(cov_mat = gl_cov_mat, v_coef = gl_v_coef, shape = gl_shape)
  gl_coef_psa <- coef_psa
  gl_shape_psa <- shape_psa

  chol_coef_surv(cov_mat = gs_mort_cov_mat, v_coef = gs_mort_v_coef,
                 shape = gs_mort_shape)
  gs_mort_coef_psa <- coef_psa
  gs_mort_shape_psa <- shape_psa

  chol_coef_surv(cov_mat = dial_mort_cov_mat, v_coef = dial_mort_v_coef,
                 shape = dial_mort_shape)
  dial_mort_coef_psa <- coef_psa
  dial_mort_shape_psa <- shape_psa

  chol_coef_logit(cov_mat = gl_mort_cov_mat, v_coef = gl_mort_v_coef)
  gl_mort_coef_psa <- coef_psa

  chol_coef_logit(cov_mat = mlogit_cov_mat, v_coef = mlogit_v_coef)
  mlogit_coef_psa <- coef_psa
  gl30_coef_psa   <- coef_psa[1:18]
  dgf_coef_psa    <- coef_psa[19:36]
  mort30_coef_psa <- coef_psa[37:54]

  return(list(gl_coef_psa = gl_coef_psa, gl_shape_psa = gl_shape_psa,
              gs_mort_coef_psa = gs_mort_coef_psa,
              gs_mort_shape_psa = gs_mort_shape_psa,
              dial_mort_coef_psa = dial_mort_coef_psa,
              dial_mort_shape_psa = dial_mort_shape_psa,
              gl_mort_coef_psa = gl_mort_coef_psa,
              mlogit_coef_psa = mlogit_coef_psa))
}

## Function to clean data between SIR runs
clean_SIR <- function() {
  sim_df <<- sim_df %>%
    dplyr::select(-c(constant, rec_cold_isch_tm)) %>%
    rename(months_to_gl = gl_time) %>%
    rename(months_to_gs_mort = gs_death_time) %>%
    rename(months_to_gl_mort = gl_death_time) %>%
    mutate(can_age_at_listing_c = rec_age_at_tx_c - (months_to_event/12)) %>%
    mutate(baseline_yrs_dial = years_dial - (months_to_event/12)) %>%
    mutate(list_year_c = tx_year_c - (months_to_event/12)) %>%
    mutate(start_time = ifelse(can_age_at_listing_c >=65, 0,
                               abs(can_age_at_listing_c)-0.5)) %>%
    mutate(months_to_event = ifelse(months_to_event < 0, 0, months_to_event)) %>%
    mutate(ddtx = ifelse(ddtx==1 & months_to_event>240 & sim=="Simulated",
                         0, ddtx)) %>%
    mutate(ldtx = ifelse(ldtx==1 & months_to_event>240 & sim=="Simulated",
                         0, ldtx)) %>%
    mutate(mort = ifelse(mort==1 & months_to_event>240 & sim=="Simulated",
                         0, mort)) %>%
    mutate(remove = ifelse(remove==1 & months_to_event>240 & sim=="Simulated",
                           0, remove)) %>%
    mutate(months_to_event = ifelse(months_to_event>240 & sim=="Simulated",
                                    240, months_to_event)) %>%
    mutate(graft_loss = ifelse(graft_loss==1 & months_to_gl>120,
                               0, graft_loss)) %>%
    mutate(gs_death = ifelse(gs_death==1 & months_to_gs_mort>120,
                             0, gs_death)) %>%
    mutate(gl_death = ifelse(gl_death==1 & months_to_gl_mort>120,
                             0, gl_death)) %>%
    mutate(months_to_gl = ifelse(months_to_gl>120,
                                 120, months_to_gl)) %>%
    mutate(months_to_gs_mort = ifelse(months_to_gs_mort>120,
                                      120, months_to_gs_mort)) %>%
    mutate(months_to_gl_mort = ifelse(months_to_gl_mort>120,
                                      120, months_to_gl_mort)) %>%
    mutate(tx_outcome = ifelse(ddtx!=1, NA, tx_outcome)) %>%
    filter(!(ddtx==1 & rec_age_at_tx_c < 0))
}

clean_sim <- function(data = sim_df) {
  tmp_df <<- data %>%
    dplyr::select(-c(optn_reg2, optn_reg3,
                     optn_reg4, optn_reg5, optn_reg6, optn_reg7, optn_reg8,
                     optn_reg9, optn_reg10, optn_reg11)) %>%
    rename(months_to_gl = gl_time) %>%
    rename(months_to_gs_mort = gs_death_time) %>%
    rename(months_to_gl_mort = gl_death_time) %>%
    mutate(can_age_at_listing_c = rec_age_at_tx_c - (months_to_event/12)) %>%
    mutate(baseline_yrs_dial = years_dial - (months_to_event/12)) %>%
    mutate(list_year_c = tx_year_c - (months_to_event/12)) %>%
    mutate(start_time = ifelse(can_age_at_listing_c >=65, 0,
                               abs(can_age_at_listing_c)-0.5)) %>%
    mutate(months_to_event = ifelse(months_to_event < 0, 0, months_to_event)) %>%
    mutate(ddtx = ifelse(ddtx==1 & months_to_event>240,
                         0, ddtx)) %>%
    mutate(ldtx = ifelse(ldtx==1 & months_to_event>240,
                         0, ldtx)) %>%
    mutate(mort = ifelse(mort==1 & months_to_event>240,
                         0, mort)) %>%
    mutate(remove = ifelse(remove==1 & months_to_event>240,
                           0, remove)) %>%
    mutate(months_to_event = ifelse(months_to_event>240,
                                    240, months_to_event)) %>%
    mutate(graft_loss = ifelse(graft_loss==1 & months_to_gl>120,
                               0, graft_loss)) %>%
    mutate(gs_death = ifelse(gs_death==1 & months_to_gs_mort>120,
                             0, gs_death)) %>%
    mutate(gl_death = ifelse(gl_death==1 & months_to_gl_mort>120,
                             0, gl_death)) %>%
    mutate(months_to_gl = ifelse(months_to_gl>120,
                                 120, months_to_gl)) %>%
    mutate(months_to_gs_mort = ifelse(months_to_gs_mort>120,
                                      120, months_to_gs_mort)) %>%
    mutate(months_to_gl_mort = ifelse(months_to_gl_mort>120,
                                      120, months_to_gl_mort))

  return(tmp_df)
}

get_survival_time <- function(data) {
  tmp_df <- data %>%
    mutate(age_at_listing = rec_age_at_tx_c - (months_to_event/12) + 65) %>%
    mutate(age_at_event = age_at_listing + months_to_event) %>%
    mutate(age_cat = ifelse(age_at_event < 70, 0,
                            ifelse(age_at_event >= 70 & age_at_event < 75, 1, 2))) %>%
    mutate(post_tx_surv = ifelse(gs_death==1, months_to_gs_mort, months_to_gl)) %>%
    mutate(post_tx_surv = ifelse(gl_death==1, months_to_gl_mort + months_to_gl, post_tx_surv)) %>%
    mutate(survtime = ifelse(ddtx==1, post_tx_surv + months_to_event, months_to_event)) %>%
    mutate(survtime = ifelse(ldtx==1 & age_cat==0 & male==0, months_to_event + surv.ldtx.65.f,
                             ifelse(ldtx==1 & age_cat==1 & male==0, months_to_event + surv.ldtx.70.f,
                                    ifelse(ldtx==1 & age_cat==2 & male==0, months_to_event + surv.ldtx.75.f, survtime)))) %>%
    mutate(survtime = ifelse(ldtx==1 & age_cat==0 & male==1, months_to_event + surv.ldtx.65.m,
                             ifelse(ldtx==1 & age_cat==1 & male==1, months_to_event + surv.ldtx.70.m,
                                    ifelse(ldtx==1 & age_cat==2 & male==1, months_to_event + surv.ldtx.75.m, survtime)))) %>%
    mutate(survtime = ifelse(remove==1 & age_cat==0 & male==0, months_to_event + surv.remove.65.f,
                             ifelse(remove==1 & age_cat==1 & male==0, months_to_event + surv.remove.70.f,
                                    ifelse(remove==1 & age_cat==2 & male==0, months_to_event + surv.remove.75.f, survtime)))) %>%
    mutate(survtime = ifelse(remove==1 & age_cat==0 & male==1, months_to_event + surv.remove.65.m,
                             ifelse(remove==1 & age_cat==1 & male==1, months_to_event + surv.remove.70.m,
                                    ifelse(remove==1 & age_cat==2 & male==1, months_to_event + surv.remove.75.m, survtime)))) %>%
    mutate(survtime = ifelse(mort==1, months_to_event, survtime)) %>%
    mutate(no_tx = ifelse(mort==1 | remove==1, 1, 0)) %>%
    mutate(kdpi_cat = ifelse(kdpi <= 20, 0,
                             ifelse(kdpi > 20 & kdpi <= 34, 1,
                                    ifelse(kdpi >= 35 & kdpi <= 85, 2, 3)))) %>%
    mutate(grp = ifelse(no_tx == 1, 0,
                        ifelse(kdpi_cat==0, 1,
                               ifelse(kdpi_cat==1, 2,
                                      ifelse(kdpi_cat==2, 3,
                                             ifelse(kdpi_cat==3, 4, NA)))))) %>%
    mutate(grp = factor(grp, levels = c(0, 1, 2, 3, 4),
                        labels = c("No Transplant", "KDPI 0-20", "KDPI 21-34",
                                   "KDPI 35-85", "KDPI 86+"))) %>%
    mutate(race = ifelse(black==0 & hispanic==0 & other_race==0, 0,
                         ifelse(black==1 & hispanic==0 & other_race==0, 1,
                                ifelse(black==0 & hispanic==1 & other_race==0, 2, 3)))) %>%
    mutate(race = factor(race, levels = c(0, 1, 2, 3),
                         labels = c("White", "Black",
                                    "Hispanic", "Other People of Color"))) %>%
    mutate(list_age_cat = ifelse(can_age_at_listing_c < 5, 0,
                                 ifelse(can_age_at_listing_c >= 5 & can_age_at_listing_c < 10, 1, 2))) %>%
    mutate(list_age_cat = factor(list_age_cat, levels = c(0, 1, 2),
                                 labels = c("65-69", "70-74", "75+"))) %>%
    dplyr::select(c(survtime, grp, run, race, list_age_cat))

  return(tmp_df)
}

get_time_to_event <- function(data) {
  tmp_df <- data %>%
    mutate(graft_loss = ifelse(graft_loss==1 & months_to_gl>120,
                               0, graft_loss)) %>%
    mutate(gs_death = ifelse(gs_death==1 & months_to_gs_mort>120,
                             0, gs_death)) %>%
    mutate(diab_stat = factor(diab_stat, levels = c(0, 1),
                              labels = c("No", "Yes"))) %>%
    mutate(graft_loss = ifelse(graft_loss==1 & months_to_gl>60,
                               0, graft_loss)) %>%
    mutate(months_to_gl_mort = ifelse(months_to_gl_mort>60,
                                      60, months_to_gl_mort)) %>%
    dplyr::select(c(months_to_event, ddtx, ldtx, mort, remove, run, graft_loss,
                    gs_death, diab_stat, months_to_gl_mort, gl_death, race))

  return(tmp_df)
}


## Function to get costs and qalys
f_costs_qalys <- function(data, r) {
  # sum(z/((1+r)^(t-1)))
  # z: vector of costs
  # r: discount rate
  # t: vector of months

  # if we use age-dependent costs
  # v_age <- seq(from = data$age_at_listing[i], to = 100, by = 1/12)
  # v_qaly_gf_1_12 <- c(rep(u_posttx_03, 2), rep(u_posttx_48, 5),
  # rep(u_posttx_912, 4))
  cost_health    <- rep(NULL, nrow(data))
  cost_societal  <- rep(NULL, nrow(data))
  qalys          <- rep(NULL, nrow(data))


  for(i in seq_len(nrow(data))) {
    male           <- data$male[i]
    age_list       <- data$can_age_at_listing_c[i] + 65
    age_list_event <- data$can_age_at_listing_c[i] + 65 +
                      (data$months_to_event[i] / 12)
    ## Waitlist Costs & QALYs
    l_cq <- get_costs_qalys_dial(age_dial = age_list,
                                 dialysis_time = data$months_to_event[i] - 1,
                                 sex = male)
    v_costs_health   <- l_cq$v_cost_health
    v_costs_societal <- l_cq$v_cost_soc
    v_qalys          <- l_cq$v_qaly

    ## Deceased Donor Costs & QALYs
    # 30-Day Outcomes
    if(data$ddtx[i] == 1) {
      pos <- get_position_qaly_vector(age_list_event)
      if(data$tx_outcome[i] == 0) {  # If graft success
        v_costs_health   <- c(v_costs_health, c_ddtx + c_oacc)
        v_costs_societal <- c(v_costs_societal, c_ddtx + c_oacc)
        if (male == 1) {
          v_qalys        <- c(v_qalys, v_u_posttx_03_m[pos])
        } else {
          v_qalys        <- c(v_qalys, v_u_posttx_03_f[pos])
        }
      } else if(data$tx_outcome[i] == 1) {  # If graft failure
        v_costs_health   <- c(v_costs_health, c_ddtx + c_oacc + c_graftloss)
        v_costs_societal <- c(v_costs_societal, c_ddtx + c_oacc + c_graftloss)
        if (male == 1) {
          v_qalys        <- c(v_qalys, v_u_dialysis_m[pos])
        } else {
          v_qalys        <- c(v_qalys, v_u_dialysis_f[pos])
        }
      } else if(data$tx_outcome[i] == 2) {  # If delayed graft function
        v_costs_health   <- c(v_costs_health, c_ddtx + c_oacc + c_dgf)
        v_costs_societal <- c(v_costs_societal, c_ddtx + c_oacc + c_dgf)
        if (male == 1) {
          v_qalys        <- c(v_qalys, 0) # updated to be 0 for month of event
        } else {
          v_qalys        <- c(v_qalys, 0) # updated to be 0 for month of event
        }
      } else if(data$tx_outcome[i] == 3) {  # If death
        v_costs_health   <- c(v_costs_health, c_ddtx + c_oacc)
        v_costs_societal <- c(v_costs_societal, c_ddtx + c_oacc)
        v_qalys          <- c(v_qalys, 0)
      }

      # Post-Tx (after 30 days)
      l_post_cq <- get_costs_qalys_posttx(age_tx = age_list_event,
                                          gs_time = max(data$gl_time[i] - 1, 0),
                                          sex = male,
                                          graft_loss = data$graft_loss[i],
                                          gl_death_time = data$gl_death_time[i],
                                          c_graftloss = c_graftloss,
                                          c_giver_time = v_c_giver_posttx,
                                          v_u_posttx_03_m = v_u_posttx_03_m,
                                          v_u_posttx_48_m = v_u_posttx_48_m,
                                          v_u_posttx_912_m = v_u_posttx_912_m,
                                          v_u_posttx_13_m = v_u_posttx_13_m,
                                          v_u_posttx_03_f = v_u_posttx_03_f,
                                          v_u_posttx_48_f = v_u_posttx_48_f,
                                          v_u_posttx_912_f = v_u_posttx_912_f,
                                          v_u_posttx_13_f = v_u_posttx_13_f)
      v_costs_health   <- c(v_costs_health, l_post_cq$v_cost_health)
      v_costs_societal <- c(v_costs_societal, l_post_cq$v_cost_soc)
      v_qalys          <- c(v_qalys, l_post_cq$v_qaly)

    } else if(data$ldtx[i] == 1) {
      pos <- get_position_qaly_vector(age_list_event)
      le_ldtx <- max(round(rgamma(n=1, shape = (120^2)/1296,
                                  scale = (1296/120))), 1)
      l_post_cq <- get_costs_qalys_posttx(age_tx = age_list_event,
                                          gs_time = le_ldtx - 1,
                                          sex = male,
                                          graft_loss = 0,
                                          v_posttx_cost = v_c_posttx,
                                          v_dial_cost = v_c_dialysis,
                                          v_dial_qaly_m = v_u_dialysis_m,
                                          v_dial_qaly_f = v_u_dialysis_f,
                                          c_patient_time = v_c_patient_posttx,
                                          c_giver_time = v_c_giver_posttx,
                                          c_patient_time_gl = c_patient_dialysis,
                                          c_giver_time_gl = c_giver_dialysis,
                                          c_graftloss = c_graftloss,
                                          v_u_posttx_03_m = v_u_posttx_03_m,
                                          v_u_posttx_48_m = v_u_posttx_48_m,
                                          v_u_posttx_912_m = v_u_posttx_912_m,
                                          v_u_posttx_13_m = v_u_posttx_13_m,
                                          v_u_posttx_03_f = v_u_posttx_03_f,
                                          v_u_posttx_48_f = v_u_posttx_48_f,
                                          v_u_posttx_912_f = v_u_posttx_912_f,
                                          v_u_posttx_13_f = v_u_posttx_13_f)
      v_costs_health   <- c(v_costs_health, c_ldtx, l_post_cq$v_cost_health)
      v_costs_societal <- c(v_costs_societal, c_ldtx, l_post_cq$v_cost_soc)
      if (male == 1) {
        v_qalys        <- c(v_qalys, v_u_posttx_03_m[pos], l_post_cq$v_qaly)
      } else {
        v_qalys        <- c(v_qalys, v_u_posttx_03_f[pos], l_post_cq$v_qaly)
      }

    } else if(data$mort[i] == 1) {
      v_costs_health   <- c(v_costs_health, 0)
      v_costs_societal <- c(v_costs_societal, 0)
      v_qalys          <- c(v_qalys, 0)
    } else if(data$remove[i]==1) {
      le_remove <- max(round(rgamma(n=1, shape = (30^2)/900,
                                    scale = (900/30))), 1)
      l_cq <- get_costs_qalys_dial(age_dial = age_list,
                                   dialysis_time = le_remove - 1,
                                   sex = male,
                                   v_dial_cost = v_c_dialysis * rr_removal)

      v_costs_health   <- c(v_costs_health, l_cq$v_cost_health)
      v_costs_societal <- c(v_costs_societal, l_cq$v_cost_soc)
      v_qalys          <- c(v_qalys, l_cq$v_qaly)
    }

    v_costs_health   <- na.omit(v_costs_health)
    v_costs_societal <- na.omit(v_costs_societal)
    v_qalys          <- na.omit(v_qalys)

    v_qalys <- unlist(v_qalys)

    t_cost_health    <- seq_along(v_costs_health)
    t_costs_societal <- seq_along(v_costs_societal)
    t_qaly           <- seq_along(v_qalys)

    cost_health[i]   <- sum(v_costs_health/((1+r)^(t_cost_health-1)))
    cost_societal[i] <- sum(v_costs_societal/((1+r)^(t_costs_societal-1)))
    qalys[i]         <- sum(v_qalys/((1+r)^(t_qaly-1)))
  }


  ## Get total cost and QALYs
  tot_cost_health   <- sum(cost_health)
  tot_cost_societal <- sum(cost_societal)
  tot_qalys         <- sum(qalys)

  ## Get total cost and QALYs by subgroup
  # Diabetes
  cost_health_diab_0   <- cost_health[which(data$diab_stat == 0)]
  cost_societal_diab_0 <- cost_societal[which(data$diab_stat == 0)]
  qalys_diab_0         <- qalys[which(data$diab_stat == 0)]

  cost_health_diab_1   <- cost_health[which(data$diab_stat == 1)]
  cost_societal_diab_1 <- cost_societal[which(data$diab_stat == 1)]
  qalys_diab_1         <- qalys[which(data$diab_stat == 1)]

  tot_cost_health_diab_0    <- sum(cost_health_diab_0)
  tot_cost_societal_diab_0  <- sum(cost_societal_diab_0)
  tot_qalys_diab_0          <- sum(qalys_diab_0)

  tot_cost_health_diab_1    <- sum(cost_health_diab_1)
  tot_cost_societal_diab_1  <- sum(cost_societal_diab_1)
  tot_qalys_diab_1          <- sum(qalys_diab_1)


  # Race/Ethnicity
  cost_health_white   <- cost_health[which(data$race == 0)]
  cost_societal_white <- cost_societal[which(data$race == 0)]
  qalys_white         <- qalys[which(data$race == 0)]

  cost_health_black   <- cost_health[which(data$race == 1)]
  cost_societal_black <- cost_societal[which(data$race == 1)]
  qalys_black         <- qalys[which(data$race == 1)]

  cost_health_hispanic   <- cost_health[which(data$race == 2)]
  cost_societal_hispanic <- cost_societal[which(data$race == 2)]
  qalys_hispanic         <- qalys[which(data$race == 2)]

  cost_health_other   <- cost_health[which(data$race == 3)]
  cost_societal_other <- cost_societal[which(data$race == 3)]
  qalys_other         <- qalys[which(data$race == 3)]

  tot_cost_health_white    <- sum(cost_health_white)
  tot_cost_societal_white  <- sum(cost_societal_white)
  tot_qalys_white          <- sum(qalys_white)

  tot_cost_health_black    <- sum(cost_health_black)
  tot_cost_societal_black  <- sum(cost_societal_black)
  tot_qalys_black          <- sum(qalys_black)

  tot_cost_health_hispanic    <- sum(cost_health_hispanic)
  tot_cost_societal_hispanic  <- sum(cost_societal_hispanic)
  tot_qalys_hispanic          <- sum(qalys_hispanic)

  tot_cost_health_other    <- sum(cost_health_other)
  tot_cost_societal_other  <- sum(cost_societal_other)
  tot_qalys_other          <- sum(qalys_other)

  # Age at Listing
  cost_health_65   <- cost_health[which(data$can_age_at_listing_c < 5)]
  cost_societal_65 <- cost_societal[which(data$can_age_at_listing_c < 5)]
  qalys_65         <- qalys[which(data$can_age_at_listing_c < 5)]

  cost_health_70   <- cost_health[which(data$can_age_at_listing_c >= 5 &
                                        data$can_age_at_listing_c < 10)]
  cost_societal_70 <- cost_societal[which(data$can_age_at_listing_c >= 5 &
                                          data$can_age_at_listing_c < 10)]
  qalys_70         <- qalys[which(data$can_age_at_listing_c >= 5 &
                                  data$can_age_at_listing_c < 10)]

  cost_health_75   <- cost_health[which(data$can_age_at_listing_c >= 10)]
  cost_societal_75 <- cost_societal[which(data$can_age_at_listing_c >= 10)]
  qalys_75         <- qalys[which(data$can_age_at_listing_c >= 10)]

  tot_cost_health_65    <- sum(cost_health_65)
  tot_cost_societal_65  <- sum(cost_societal_65)
  tot_qalys_65          <- sum(qalys_65)

  tot_cost_health_70    <- sum(cost_health_70)
  tot_cost_societal_70  <- sum(cost_societal_70)
  tot_qalys_70          <- sum(qalys_70)

  tot_cost_health_75    <- sum(cost_health_75)
  tot_cost_societal_75  <- sum(cost_societal_75)
  tot_qalys_75          <- sum(qalys_75)

  ## Event counts by subgroup
  # All candidates
  ddtx_count      <- sum(data$ddtx == 1)
  ldtx_count      <- sum(data$ldtx == 1)
  wait_mort_count <- sum(data$mort == 1)
  removal_count   <- sum(data$remove == 1)
  dgf_count       <- sum(data$tx_outcome == 2, na.rm = T)
  gl_count        <- sum(data$graft_loss == 1, na.rm = T)

  error_count <- function(outcome, outcome_val, subgroup, subgroup_val) {
    tryCatch(expr = {
      tmp_count <- sum(outcome == outcome_val & subgroup == subgroup_val, na.rm = T)
    }, error = function(e) {tmp_count <- NaN}
    )
    return(tmp_count)
  }
  # By diabetes status
  ddtx_count_diab0      <- error_count(data$ddtx, 1, data$diab_stat, 0)
  ldtx_count_diab0      <- error_count(data$ldtx, 1, data$diab_stat, 0)
  wait_mort_count_diab0 <- error_count(data$mort, 1, data$diab_stat, 0)
  removal_count_diab0   <- error_count(data$remove, 1, data$diab_stat, 0)
  dgf_count_diab0       <- error_count(data$tx_outcome, 2, data$diab_stat, 0)
  gl_count_diab0        <- error_count(data$graft_loss, 1, data$diab_stat, 0)

  ddtx_count_diab1      <- error_count(data$ddtx, 1, data$diab_stat, 1)
  ldtx_count_diab1      <- error_count(data$ldtx, 1, data$diab_stat, 1)
  wait_mort_count_diab1 <- error_count(data$mort, 1, data$diab_stat, 1)
  removal_count_diab1   <- error_count(data$remove, 1, data$diab_stat, 1)
  dgf_count_diab1       <- error_count(data$tx_outcome, 2, data$diab_stat, 1)
  gl_count_diab1        <- error_count(data$graft_loss, 1, data$diab_stat, 1)

  # ddtx_count_diab0      <- sum(data$ddtx == 1 & data$diab_stat == 0)
  # ldtx_count_diab0      <- sum(data$ldtx == 1 & data$diab_stat == 0)
  # wait_mort_count_diab0 <- sum(data$mort == 1 & data$diab_stat == 0)
  # removal_count_diab0   <- sum(data$remove == 1 & data$diab_stat == 0)
  # dgf_count_diab0       <- sum(data$tx_outcome == 2 & data$diab_stat == 0)
  # gl_count_diab0        <- sum(data$graft_loss == 1 & data$diab_stat == 0)

  # ddtx_count_diab1      <- sum(data$ddtx == 1 & data$diab_stat == 1)
  # ldtx_count_diab1      <- sum(data$ldtx == 1 & data$diab_stat == 1)
  # wait_mort_count_diab1 <- sum(data$mort == 1 & data$diab_stat == 1)
  # removal_count_diab1   <- sum(data$remove == 1 & data$diab_stat == 1)
  # dgf_count_diab1       <- sum(data$tx_outcome == 2 & data$diab_stat == 1)
  # gl_count_diab1        <- sum(data$graft_loss == 1 & data$diab_stat == 1)

  # By race/ethnicity
  ddtx_count_white      <- error_count(data$ddtx, 1, data$race, 0)
  ldtx_count_white      <- error_count(data$ldtx, 1, data$race, 0)
  wait_mort_count_white <- error_count(data$mort, 1, data$race, 0)
  removal_count_white   <- error_count(data$remove, 1, data$race, 0)
  dgf_count_white       <- error_count(data$tx_outcome, 2, data$race, 0)
  gl_count_white        <- error_count(data$graft_loss, 1, data$race, 0)

  ddtx_count_black      <- error_count(data$ddtx, 1, data$race, 1)
  ldtx_count_black      <- error_count(data$ldtx, 1, data$race, 1)
  wait_mort_count_black <- error_count(data$mort, 1, data$race, 1)
  removal_count_black   <- error_count(data$remove, 1, data$race, 1)
  dgf_count_black       <- error_count(data$tx_outcome, 2, data$race, 1)
  gl_count_black        <- error_count(data$graft_loss, 1, data$race, 1)

  ddtx_count_hisp      <- error_count(data$ddtx, 1, data$race, 2)
  ldtx_count_hisp      <- error_count(data$ldtx, 1, data$race, 2)
  wait_mort_count_hisp <- error_count(data$mort, 1, data$race, 2)
  removal_count_hisp   <- error_count(data$remove, 1, data$race, 2)
  dgf_count_hisp       <- error_count(data$tx_outcome, 2, data$race, 2)
  gl_count_hisp        <- error_count(data$graft_loss, 1, data$race, 2)

  ddtx_count_other      <- error_count(data$ddtx, 1, data$race, 3)
  ldtx_count_other      <- error_count(data$ldtx, 1, data$race, 3)
  wait_mort_count_other <- error_count(data$mort, 1, data$race, 3)
  removal_count_other   <- error_count(data$remove, 1, data$race, 3)
  dgf_count_other       <- error_count(data$tx_outcome, 2, data$race, 3)
  gl_count_other        <- error_count(data$graft_loss, 1, data$race, 3)

  # ddtx_count_white      <- sum(data$ddtx == 1 & data$race == 0)
  # ldtx_count_white      <- sum(data$ldtx == 1 & data$race == 0)
  # wait_mort_count_white <- sum(data$mort == 1 & data$race == 0)
  # removal_count_white   <- sum(data$remove == 1 & data$race == 0)
  # dgf_count_white       <- sum(data$tx_outcome == 2 & data$race == 0)
  # gl_count_white        <- sum(data$graft_loss == 1 & data$race == 0)

  # ddtx_count_black      <- sum(data$ddtx == 1 & data$race == 1)
  # ldtx_count_black      <- sum(data$ldtx == 1 & data$race == 1)
  # wait_mort_count_black <- sum(data$mort == 1 & data$race == 1)
  # removal_count_black   <- sum(data$remove == 1 & data$race == 1)
  # dgf_count_black       <- sum(data$tx_outcome == 2 & data$race == 1)
  # gl_count_black        <- sum(data$graft_loss == 1 & data$race == 1)

  # ddtx_count_hisp      <- sum(data$ddtx == 1 & data$race == 2)
  # ldtx_count_hisp      <- sum(data$ldtx == 1 & data$race == 2)
  # wait_mort_count_hisp <- sum(data$mort == 1 & data$race == 2)
  # removal_count_hisp   <- sum(data$remove == 1 & data$race == 2)
  # dgf_count_hisp       <- sum(data$tx_outcome == 2 & data$race == 2)
  # gl_count_hisp        <- sum(data$graft_loss == 1 & data$race == 2)

  # ddtx_count_other      <- sum(data$ddtx == 1 & data$race == 3)
  # ldtx_count_other      <- sum(data$ldtx == 1 & data$race == 3)
  # wait_mort_count_other <- sum(data$mort == 1 & data$race == 3)
  # removal_count_other   <- sum(data$remove == 1 & data$race == 3)
  # dgf_count_other       <- sum(data$tx_outcome == 2 & data$race == 3)
  # gl_count_other        <- sum(data$graft_loss == 1 & data$race == 3)

  # By age at listing
  data$age_grp <- ifelse(data$can_age_at_listing_c < 5, 1,
                          ifelse(data$can_age_at_listing_c >= 5 & data$can_age_at_listing_c < 10, 2, 3))

  ddtx_count_65      <- error_count(data$ddtx, 1, data$age_grp, 1)
  ldtx_count_65      <- error_count(data$ldtx, 1, data$age_grp, 1)
  wait_mort_count_65 <- error_count(data$mort, 1, data$age_grp, 1)
  removal_count_65   <- error_count(data$remove, 1, data$age_grp, 1)
  dgf_count_65       <- error_count(data$tx_outcome, 2, data$age_grp, 1)
  gl_count_65        <- error_count(data$graft_loss, 1, data$age_grp, 1)

  ddtx_count_70      <- error_count(data$ddtx, 1, data$age_grp, 2)
  ldtx_count_70      <- error_count(data$ldtx, 1, data$age_grp, 2)
  wait_mort_count_70 <- error_count(data$mort, 1, data$age_grp, 2)
  removal_count_70   <- error_count(data$remove, 1, data$age_grp, 2)
  dgf_count_70       <- error_count(data$tx_outcome, 2, data$age_grp, 2)
  gl_count_70        <- error_count(data$graft_loss, 1, data$age_grp, 2)

  ddtx_count_75      <- error_count(data$ddtx, 1, data$age_grp, 3)
  ldtx_count_75      <- error_count(data$ldtx, 1, data$age_grp, 3)
  wait_mort_count_75 <- error_count(data$mort, 1, data$age_grp, 3)
  removal_count_75   <- error_count(data$remove, 1, data$age_grp, 3)
  dgf_count_75       <- error_count(data$tx_outcome, 2, data$age_grp, 3)
  gl_count_75        <- error_count(data$graft_loss, 1, data$age_grp, 3)

  # ddtx_count_65      <- sum(data$ddtx == 1 &
  #                                 data$can_age_at_listing_c < 5)
  # ldtx_count_65      <- sum(data$ldtx == 1 &
  #                                 data$can_age_at_listing_c < 5)
  # wait_mort_count_65 <- sum(data$mort == 1 &
  #                                 data$can_age_at_listing_c < 5)
  # removal_count_65   <- sum(data$remove == 1 &
  #                                 data$can_age_at_listing_c < 5)
  # dgf_count_65       <- sum(data$tx_outcome == 2 &
  #                                 data$can_age_at_listing_c < 5)
  # gl_count_65        <- sum(data$graft_loss == 1 &
  #                                 data$can_age_at_listing_c < 5)

  # ddtx_count_70      <- sum(data$ddtx == 1 &
  #                                 (data$can_age_at_listing_c >= 5 &
  #                                 data$can_age_at_listing_c < 10))
  # ldtx_count_70      <- sum(data$ldtx == 1 &
  #                                 (data$can_age_at_listing_c >= 5 &
  #                                 data$can_age_at_listing_c < 10))
  # wait_mort_count_70 <- sum(data$mort == 1 &
  #                                 (data$can_age_at_listing_c >= 5 &
  #                                 data$can_age_at_listing_c < 10))
  # removal_count_70   <- sum(data$remove == 1 &
  #                                 (data$can_age_at_listing_c >= 5 &
  #                                 data$can_age_at_listing_c < 10))
  # dgf_count_70       <- sum(data$tx_outcome == 2 &
  #                                 (data$can_age_at_listing_c >= 5 &
  #                                 data$can_age_at_listing_c < 10))
  # gl_count_70        <- sum(data$graft_loss == 1 &
  #                                 (data$can_age_at_listing_c >= 5 &
  #                                 data$can_age_at_listing_c < 10))

  # ddtx_count_75      <- sum(data$ddtx == 1 &
  #                                 data$can_age_at_listing_c >= 10)
  # ldtx_count_75      <- sum(data$ldtx == 1 &
  #                                 data$can_age_at_listing_c >= 10)
  # wait_mort_count_75 <- sum(data$mort == 1 &
  #                                 data$can_age_at_listing_c >= 10)
  # removal_count_75   <- sum(data$remove == 1 &
  #                                 data$can_age_at_listing_c >= 10)
  # dgf_count_75       <- sum(data$tx_outcome == 2 &
  #                                 data$can_age_at_listing_c >= 10)
  # gl_count_75        <- sum(data$graft_loss == 1 &
  #                                 data$can_age_at_listing_c >= 10)

  event_counts <- c(ddtx_count, ldtx_count, wait_mort_count,
                    removal_count, dgf_count, gl_count,
                    ddtx_count_diab0, ldtx_count_diab0, wait_mort_count_diab0,
                    removal_count_diab0, dgf_count_diab0, gl_count_diab0,
                    ddtx_count_diab1, ldtx_count_diab1, wait_mort_count_diab1,
                    removal_count_diab1, dgf_count_diab1, gl_count_diab1,
                    ddtx_count_white, ldtx_count_white, wait_mort_count_white,
                    removal_count_white, dgf_count_white, gl_count_white,
                    ddtx_count_black, ldtx_count_black, wait_mort_count_black,
                    removal_count_black, dgf_count_black, gl_count_black,
                    ddtx_count_hisp, ldtx_count_hisp, wait_mort_count_hisp,
                    removal_count_hisp, dgf_count_hisp, gl_count_hisp,
                    ddtx_count_other, ldtx_count_other, wait_mort_count_other,
                    removal_count_other, dgf_count_other, gl_count_other,
                    ddtx_count_65, ldtx_count_65, wait_mort_count_65,
                    removal_count_65, dgf_count_65, gl_count_65,
                    ddtx_count_70, ldtx_count_70, wait_mort_count_70,
                    removal_count_70, dgf_count_70, gl_count_70,
                    ddtx_count_75, ldtx_count_75, wait_mort_count_75,
                    removal_count_75, dgf_count_75, gl_count_75) # dgf_count_75, gl_count_75

  cost_qalys <- c(tot_cost_health, tot_cost_societal, tot_qalys, tot_cost_health_diab_0,
                  tot_cost_societal_diab_0, tot_qalys_diab_0, tot_cost_health_diab_1,
                  tot_cost_societal_diab_1, tot_qalys_diab_1, tot_cost_health_white,
                  tot_cost_societal_white, tot_qalys_white, tot_cost_health_black,
                  tot_cost_societal_black, tot_qalys_black, tot_cost_health_hispanic,
                  tot_cost_societal_hispanic, tot_qalys_hispanic, tot_cost_health_other,
                  tot_cost_societal_other, tot_qalys_other, tot_cost_health_65,
                  tot_cost_societal_65, tot_qalys_65, tot_cost_health_70,
                  tot_cost_societal_70, tot_qalys_70, tot_cost_health_75,
                  tot_cost_societal_75, tot_qalys_75, event_counts)
  return(cost_qalys)
}

get_costs_qalys_dial <- function(age_dial, dialysis_time,
                                 v_dial_cost = v_c_dialysis,
                                 v_dial_qaly_m = v_u_dialysis_m,
                                 v_dial_qaly_f = v_u_dialysis_f,
                                 sex,
                                 c_patient_time = c_patient_dialysis,
                                 c_giver_time = c_giver_dialysis) {

  if (dialysis_time == 0) {
    v_pos_cost_wait <- get_position_cost_vector(age_dial)
    v_pos_qaly_wait <- get_position_qaly_vector(age_dial)
    v_cost_health   <- get_cq_vector(v_pos_cost_wait, v_dial_cost)
    v_cost_soc      <- v_cost_health + c_patient_time + c_giver_time
    v_qaly          <- c()
  } else {
    max_age_list    <- age_dial + ((dialysis_time - 1)/12)
    v_age_wait      <- seq(from = age_dial, to = max_age_list, by = (1/12))
    v_pos_cost_wait <- get_position_cost_vector(v_age_wait)
    v_pos_qaly_wait <- get_position_qaly_vector(v_age_wait)
    v_cost_health   <- get_cq_vector(v_pos_cost_wait, v_dial_cost)
    v_cost_soc      <- v_cost_health + c_patient_time + c_giver_time
    if(sex==1){
      v_qaly <- get_cq_vector(v_pos_qaly_wait, v_dial_qaly_m)
    } else {
      v_qaly <- get_cq_vector(v_pos_qaly_wait, v_dial_qaly_f)
    }
  }

  return(list(v_cost_soc = v_cost_soc,
              v_cost_health = v_cost_health,
              v_qaly = v_qaly))
}

get_costs_qalys_posttx <- function(age_tx, gs_time,
                                   graft_loss, gl_death_time = NULL,
                                   v_posttx_cost = v_c_posttx,
                                   sex,
                                   v_dial_cost = v_c_dialysis,
                                   v_dial_qaly_m = v_u_dialysis_m,
                                   v_dial_qaly_f = v_u_dialysis_f,
                                   c_patient_time = v_c_patient_posttx,
                                   c_giver_time = v_c_giver_posttx,
                                   c_patient_time_gl = c_patient_dialysis,
                                   c_giver_time_gl = c_giver_dialysis,
                                   c_graftloss,
                                   v_u_posttx_03_m,
                                   v_u_posttx_48_m,
                                   v_u_posttx_912_m,
                                   v_u_posttx_13_m,
                                   v_u_posttx_03_f,
                                   v_u_posttx_48_f,
                                   v_u_posttx_912_f,
                                   v_u_posttx_13_f) {

  # Functioning Graft
  max_age         <- age_tx + (gs_time/12)
  v_age           <- seq(from = age_tx, to = max(age_tx, max_age - (1/12)),
                         by = (1/12))
  v_time_gs       <- seq(from = 1, to = max(1, gs_time), by = 1)
  v_pos_cost      <- get_position_cost_vector(v_age)
  v_pos_qaly      <- get_position_qaly_vector_posttx(v_age, v_time_gs)
  v_pos_patient   <- get_position_patient_time_vector(v_time_gs)
  v_cost_health   <- get_cq_vector(v_pos_cost, v_posttx_cost)
  v_cost_patient  <- get_cq_vector(v_pos_patient, c_patient_time)
  v_cost_giver    <- get_cq_vector(v_pos_patient, c_giver_time)
  v_cost_soc      <- v_cost_health + v_cost_patient + v_cost_giver

  if (sex == 1) {
    v_qaly <- get_q_vector(vec_age = v_pos_qaly$vec_age,
                           vec_gs = v_pos_qaly$vec_gs,
                           v_q_03 = v_u_posttx_03_m,
                           v_q_48 = v_u_posttx_48_m,
                           v_q_912 = v_u_posttx_912_m,
                           v_q_13 = v_u_posttx_13_m)
  } else {
    v_qaly <- get_q_vector(vec_age = v_pos_qaly$vec_age,
                           vec_gs = v_pos_qaly$vec_gs,
                           v_q_03 = v_u_posttx_03_f,
                           v_q_48 = v_u_posttx_48_f,
                           v_q_912 = v_u_posttx_912_f,
                           v_q_13 = v_u_posttx_13_f)
  }

  age_gl <- max_age


  # Graft Loss
  if (graft_loss == 0) {
    v_cost_health <- c(v_cost_health, 0)
    v_cost_soc    <- c(v_cost_soc, 0)
    v_qaly        <- c(v_qaly, 0)
  } else if (graft_loss == 1 && gl_death_time == 0) {
    v_cost_health <- c(v_cost_health, c_graftloss + 0)
    v_cost_soc    <- c(v_cost_soc, c_graftloss + 0)
    v_qaly        <- c(v_qaly, 0)
  } else if (graft_loss == 1 && gl_death_time > 0) {
    max_age_gl          <- age_gl + (gl_death_time / 12)
    v_age_gl            <- seq(from = age_gl, to = max_age_gl, by = (1/12))
    v_pos_cost_gl       <- get_position_cost_vector(v_age_gl)
    v_pos_qaly_gl       <- get_position_qaly_vector(v_age_gl)
    v_cost_gl_health    <- c()
    v_cost_gl_health    <- get_cq_vector(v_pos_cost_gl, v_dial_cost)
    v_cost_gl_health[1] <- c_graftloss
    v_cost_gl_soc       <- v_cost_gl_health + c_patient_time_gl + c_giver_time_gl
    if (sex == 1) {
      v_gl_qaly <- get_cq_vector(v_pos_qaly_gl, v_dial_qaly_m)
    } else {
      v_gl_qaly <- get_cq_vector(v_pos_qaly_gl, v_dial_qaly_f)
    }
    v_gl_qaly[1] <- 0  # updated to be 0 at month of event

    v_cost_health <- c(v_cost_health, v_cost_gl_health)
    v_cost_soc    <- c(v_cost_soc, v_cost_gl_soc)
    v_qaly        <- c(v_qaly, v_gl_qaly)

  }

  return(list(v_cost_soc = v_cost_soc,
              v_cost_health = v_cost_health,
              v_qaly = v_qaly))
}

get_position_cost_vector <- function(vec) {
  vec <- replace(vec, vec < 70, 1)
  vec <- replace(vec, vec < 75 & vec >= 70, 2)
  vec <- replace(vec, vec < 80 & vec >= 70, 3)
  vec <- replace(vec, vec < 85 & vec >= 70, 4)
  vec <- replace(vec, vec >= 85, 5)
  return(vec)
}

get_position_patient_time_vector <- function(vec) {
  vec <- replace(vec, vec <= 1, 1)
  vec <- replace(vec, vec >= 2 & vec <= 3, 2)
  vec <- replace(vec, vec >= 4, 3)
  return(vec)
}

get_position_qaly_vector <- function(vec) {
  vec <- replace(vec, vec < 70, 1)
  vec <- replace(vec, vec < 80 & vec >= 70, 2)
  vec <- replace(vec, vec >= 80, 3)
  return(vec)
}

get_position_qaly_vector_posttx <- function(vec_age, vec_gs) {
  vec_age <- replace(vec_age, vec_age<70, 1)
  vec_age <- replace(vec_age, vec_age<80 & vec_age>=70, 2)
  vec_age <- replace(vec_age, vec_age>=80, 3)

  vec_gs  <- replace(vec_gs, vec_gs <= 3, 1)
  vec_gs  <- replace(vec_gs, vec_gs >= 4 & vec_gs <= 8, 2)
  vec_gs  <- replace(vec_gs, vec_gs >= 9 & vec_gs <= 12, 3)
  vec_gs  <- replace(vec_gs, vec_gs >= 13, 4)
  return(list(vec_age = vec_age,
              vec_gs = vec_gs))
}


get_cq_vector <- function(vec, v_ref) {
  for(i in min(vec):max(vec)){
    vec <- replace(vec, vec==i, v_ref[i])
  }
  return(vec)
}

# get_q_vector <- function(vec_age, vec_gs, v_q_03, v_q_48, v_q_912, v_q_13) {
#   # print(length(vec_age) == length(vec_gs))
#   for(i in min(vec_age):max(vec_age)){
#     if(vec_gs[i]==1){
#       vec_age <- replace(vec_age, vec_age==i, v_q_03[i])
#     } else if(vec_gs[i]==2) {
#       vec_age <- replace(vec_age, vec_age==i, v_q_48[i])
#     } else if(vec_gs[i]==3) {
#       vec_age <- replace(vec_age, vec_age==i, v_q_912[i])
#     } else if(vec_gs[i]==4) {
#       vec_age <- replace(vec_age, vec_age==i, v_q_13[i])
#     }
#   }
#   return(vec_age)
# }

get_q_vector <- function(vec_age, vec_gs, v_q_03, v_q_48, v_q_912, v_q_13) {
  # print(length(vec_age) == length(vec_gs))
  vec_q <- c()
  for (i in seq_len(length(vec_gs))) {
    if (vec_age[i] == 1) {
      if (vec_gs[i] == 1) {
        vec_q <- c(vec_q, v_q_03[1])
      } else if (vec_gs[i] == 2) {
        vec_q <- c(vec_q, v_q_48[1])
      } else if (vec_gs[i] == 3) {
        vec_q <- c(vec_q, v_q_912[1])
      } else if (vec_gs[i] == 4) {
        vec_q <- c(vec_q, v_q_13[1])
      }
    } else if (vec_age[i] == 2) {
      if (vec_gs[i] == 1) {
        vec_q <- c(vec_q, v_q_03[2])
      } else if (vec_gs[i] == 2) {
        vec_q <- c(vec_q, v_q_48[2])
      } else if (vec_gs[i] == 3) {
        vec_q <- c(vec_q, v_q_912[2])
      } else if (vec_gs[i] == 4) {
        vec_q <- c(vec_q, v_q_13[2])
      }
    } else if (vec_age[i] == 3) {
      if (vec_gs[i] == 1) {
        vec_q <- c(vec_q, v_q_03[3])
      } else if (vec_gs[i] == 2) {
        vec_q <- c(vec_q, v_q_48[3])
      } else if (vec_gs[i] == 3) {
        vec_q <- c(vec_q, v_q_912[3])
      } else if (vec_gs[i] == 4) {
        vec_q <- c(vec_q, v_q_13[3])
      }
    }
  }
  return(vec_age)
}

CorrelateUtils <- function(U, Q, epsilon, delta){
  n <- nrow(U) #number of PSA samples
  s <- ncol(U) #number of states
  R <- matrix(rnorm(n*s,0,1), n, s) #the reference matrix.
  C <- matrix(0,s,s) #a place holder for the correlation matrix
  Viol <- matrix(0,s,s) #violations matrix.
  for (j in 2:s){ # {j,k} is the selected pair of state comparisons
    for (k in 1:(j-1)){
      rho <- 1 # #bivariate correlations
      X = U[,c(j,k)] #selected columns of U
      Y = R[,c(j,k)] #selected columns of R
      viol <- 0
      while(viol<epsilon & rho>=0){ #if these conditions are met, continue
        rho <- rho - delta #reduce correltion
        Xstar = induceRankCorrelation(X, Y, rho) #correlated utilities.
        viol = mean((Q[j,k] * Xstar[,1]) < (Q[j,k] * Xstar[,2])) #compute %violations between the col. vectors.
      }
      #Viol[j,k] <- viol
      C[j,k] <- rho + delta #record the desired correlation.
    }
    # print(j) #just to show the column indices.
    # print(k)
  }
  #Fill in the other elements of C.
  C = C + t(C)
  for (j in 1:s){
    C[j,j] <- 1 # % the diagonal of ones.
  }
  ## Eigenvectors and Eigenvalues correction of C
  eigenResults <- eigen(C)
  B <- eigenResults$values
  V <- eigenResults$vectors
  B[B<=0] <- 0.0001 #to make sure C is positive definite, set eigenvalues<=0 to a very small positive number
  Cstar <- V %*% diag(B) %*% solve(V) #reconstruct C
  Ustar <- induceRankCorrelation(U, R, Cstar) #similar to above, induce the correlation.
  return(Ustar)
}
## To induce Rank correlation: inputs X: QoL vectors, Y is the reference vectors, and
## Sigma is the correlation matrix.

induceRankCorrelation <- function(X, Y, Sigma){
  if (length(Sigma)==1){ #if Sigma is a single value, convert it to a 2x2 matrix.
    Sigma <- matrix(c(1, Sigma,
                      Sigma, 1), 2, 2)
  }
  n <- nrow(X)
  s <- ncol(X)

  #Initialize matrices.
  Xsorted <- matrix(0, n, s)
  Yrank <- matrix(0, n, s)
  Xstar <- matrix(0, n, s)
  P <- chol(Sigma) #compute the upper triangular matrix
  Ystar <- Y %*% P #Sort the values in the reference vectors by multiplying by P
  cor(Ystar)
  for (j in 1:s){
    Xsorted[,j] <- sort(X[,j]) #Sort each variable
    Yrank[order(Ystar[,j]),j] <- seq(1:n) #Reverse sort
    Xstar[,j]=Xsorted[Yrank[,j],j] #sort Xsorted to have the same ranks as Ystar.
  }
  return(Xstar) #return the sorted vectors.
}

# Function to generate ordinal QALYs
ordered_qalys <- function(n.runs, alphas, betas) {
  n <- n.runs
  s <- 3

  m_U <- matrix(0, n, s)
  m_U[,1] <- rbeta(n, alphas[1], betas[1])
  m_U[,2] <- rbeta(n, alphas[2], betas[2])
  m_U[,3] <- rbeta(n, alphas[3], betas[3])

  m_Q <- matrix(c(0,0,0,
                  -1,0,0,
                  -1,-1,0), s, s, byrow = T)

  epsilon <- 0.05
  delta <- 0.01

  m_Ustar <- CorrelateUtils(m_U, m_Q, epsilon, delta) #Induce correlation, and return Ustar.

  return(m_Ustar)
}

# Function to generate ordinal costs
ordered_costs <- function(n.runs, means) {
  n <- n.runs
  s <- 3

  m_U <- matrix(0, n, s)
  m_U[,1] <- runif(n, min = 0.8 * means[1], max = 1.2 * means[1])
  m_U[,2] <- runif(n, min = 0.8 * means[2], max = 1.2 * means[2])
  m_U[,3] <- runif(n, min = 0.8 * means[3], max = 1.2 * means[3])

  m_Q <- matrix(c(0,0,0,
                  -1,0,0,
                  -1,-1,0), s, s, byrow = T)

  epsilon <- 0.05
  delta <- 0.01

  m_Ustar <- CorrelateUtils(m_U, m_Q, epsilon, delta) #Induce correlation, and return Ustar.

  return(m_Ustar)
}

get_psa_draws <- function(l_qaly_params, l_cost_params, n_draws) {
  # Load input parameters
  v_u_dialysis_m_alpha   = l_qaly_params$v_u_dialysis_m_alpha
  v_u_dialysis_m_beta    = l_qaly_params$v_u_dialysis_m_beta
  v_u_dialysis_f_alpha   = l_qaly_params$v_u_dialysis_f_alpha
  v_u_dialysis_f_beta    = l_qaly_params$v_u_dialysis_f_beta
  v_u_posttx_03_m_alpha  = l_qaly_params$v_u_posttx_03_m_alpha
  v_u_posttx_03_m_beta   = l_qaly_params$v_u_posttx_03_m_beta
  v_u_posttx_03_f_alpha  = l_qaly_params$v_u_posttx_03_f_alpha
  v_u_posttx_03_f_beta   = l_qaly_params$v_u_posttx_03_f_beta
  v_u_posttx_48_m_alpha  = l_qaly_params$v_u_posttx_48_m_alpha
  v_u_posttx_48_m_beta   = l_qaly_params$v_u_posttx_48_m_beta
  v_u_posttx_48_f_alpha  = l_qaly_params$v_u_posttx_48_f_alpha
  v_u_posttx_48_f_beta   = l_qaly_params$v_u_posttx_48_f_beta
  v_u_posttx_912_m_alpha = l_qaly_params$v_u_posttx_912_m_alpha
  v_u_posttx_912_m_beta  = l_qaly_params$v_u_posttx_912_m_beta
  v_u_posttx_912_f_alpha = l_qaly_params$v_u_posttx_912_f_alpha
  v_u_posttx_912_f_beta  = l_qaly_params$v_u_posttx_912_f_beta
  v_u_posttx_13_m_alpha  = l_qaly_params$v_u_posttx_13_m_alpha
  v_u_posttx_13_m_beta   = l_qaly_params$v_u_posttx_13_m_beta
  v_u_posttx_13_f_alpha  = l_qaly_params$v_u_posttx_13_f_alpha
  v_u_posttx_13_f_beta   = l_qaly_params$v_u_posttx_13_f_beta

  c_dial_65_params        = l_cost_params$c_dial_65_params
  c_dial_70_params        = l_cost_params$c_dial_70_params
  c_dial_75_params        = l_cost_params$c_dial_75_params
  c_dial_80_params        = l_cost_params$c_dial_80_params
  c_dial_85_params        = l_cost_params$c_dial_85_params
  c_posttx_65_params      = l_cost_params$c_posttx_65_params
  c_posttx_70_params      = l_cost_params$c_posttx_70_params
  c_posttx_75_params      = l_cost_params$c_posttx_75_params
  c_posttx_80_params      = l_cost_params$c_posttx_80_params
  c_posttx_85_params      = l_cost_params$c_posttx_85_params
  c_ddtx_params           = l_cost_params$c_ddtx_params
  c_oacc_params           = l_cost_params$c_oacc_params
  c_ldtx_params           = l_cost_params$c_ldtx_params
  c_dgf_params            = l_cost_params$c_dgf_params
  C_graftloss_params      = l_cost_params$c_graftloss_params
  c_giver_dialysis_params = l_cost_params$c_giver_dialysis_params
  c_giver_posttx_params   = l_cost_params$c_giver_posttx_params

  rr_removal             = l_cost_params$rr_removal
  c_patient_dialysis     = l_cost_params$c_patient_dialysis
  v_c_patient_posttx     = l_cost_params$v_c_patient_posttx
  v_c_giver_posttx       = l_cost_params$v_c_giver_posttx

  # QALY draws
  m_u_dialysis_m <- ordered_qalys(n.runs = n_draws,
                                  alphas = v_u_dialysis_m_alpha,
                                  betas = v_u_dialysis_m_beta)
  m_u_dialysis_f <- ordered_qalys(n.runs = n_draws,
                                  alphas = v_u_dialysis_f_alpha,
                                  betas = v_u_dialysis_f_beta)
  m_u_posttx_03_m <- ordered_qalys(n.runs = n_draws,
                                   alphas = v_u_posttx_03_m_alpha,
                                   betas = v_u_posttx_03_m_beta)
  m_u_posttx_03_f <- ordered_qalys(n.runs = n_draws,
                                   alphas = v_u_posttx_03_f_alpha,
                                   betas = v_u_posttx_03_f_beta)
  m_u_posttx_48_m <- ordered_qalys(n.runs = n_draws,
                                   alphas = v_u_posttx_48_m_alpha,
                                   betas = v_u_posttx_48_m_beta)
  m_u_posttx_48_f <- ordered_qalys(n.runs = n_draws,
                                   alphas = v_u_posttx_48_f_alpha,
                                   betas = v_u_posttx_48_f_beta)
  m_u_posttx_912_m <- ordered_qalys(n.runs = n_draws,
                                    alphas = v_u_posttx_912_m_alpha,
                                    betas = v_u_posttx_912_m_beta)
  m_u_posttx_912_f <- ordered_qalys(n.runs = n_draws,
                                    alphas = v_u_posttx_912_f_alpha,
                                    betas = v_u_posttx_912_f_beta)
  m_u_posttx_13_m <- ordered_qalys(n.runs = n_draws,
                                   alphas = v_u_posttx_13_m_alpha,
                                   betas = v_u_posttx_13_m_beta)
  m_u_posttx_13_f <- ordered_qalys(n.runs = n_draws,
                                   alphas = v_u_posttx_13_f_alpha,
                                   betas = v_u_posttx_13_f_beta)

  # Cost draws
  v_c_dialysis <- cbind(rgamma(n_draws,
                               shape = c_dial_65_params[1],
                               rate = c_dial_65_params[2]),
                        rgamma(n_draws,
                               shape = c_dial_70_params[1],
                               rate = c_dial_70_params[2]),
                        rgamma(n_draws,
                               shape = c_dial_75_params[1],
                               rate = c_dial_75_params[2]),
                        rgamma(n_draws,
                               shape = c_dial_80_params[1],
                               rate = c_dial_80_params[2]),
                        rgamma(n_draws,
                               shape = c_dial_85_params[1],
                               rate = c_dial_85_params[2]))
  v_c_posttx <- cbind(rgamma(n_draws,
                             shape = c_posttx_65_params[1],
                             rate = c_posttx_65_params[2]),
                      rgamma(n_draws,
                             shape = c_posttx_70_params[1],
                             rate = c_posttx_70_params[2]),
                      rgamma(n_draws,
                             shape = c_posttx_75_params[1],
                             rate = c_posttx_75_params[2]),
                      rgamma(n_draws,
                             shape = c_posttx_80_params[1],
                             rate = c_posttx_80_params[2]),
                      rgamma(n_draws,
                             shape = c_posttx_85_params[1],
                             rate = c_posttx_85_params[2]))
    c_ddtx <- rgamma(n_draws, shape = c_ddtx_params[1],
                     rate = c_ddtx_params[2])
    c_oacc <- rgamma(n_draws, shape = c_oacc_params[1],
                     rate = c_oacc_params[2])
    c_ldtx <- rgamma(n_draws, shape = c_ldtx_params[1],
                     rate = c_ldtx_params[2])
    c_dgf <- rgamma(n_draws, shape = c_dgf_params[1],
                    rate = c_dgf_params[2])
    c_graftloss <- rgamma(n_draws, shape = C_graftloss_params[1],
                          rate = C_graftloss_params[2])
    c_giver_dialysis <- rgamma(n_draws,
                               shape = c_giver_dialysis_params[1],
                               rate = c_giver_dialysis_params[2])
    c_giver_posttx <- rgamma(n_draws,
                             shape = c_giver_posttx_params[1],
                             rate = c_giver_posttx_params[2])
    c_patient_dialysis <- runif(n_draws, min = 0.8 * c_patient_dialysis,
                                max = 1.2 * c_patient_dialysis)

    v_c_patient_posttx <- ordered_costs(n.runs = n_draws,
                                        means = v_c_patient_posttx)
    v_c_giver_posttx <- ordered_costs(n.runs = n_draws,
                                      means = v_c_giver_posttx)
    rr_removal <- runif(n_draws,
                        min = rr_removal - 0.2,
                        max = rr_removal + 0.2)

    return(list(m_u_dialysis_m = m_u_dialysis_m,
                m_u_dialysis_f = m_u_dialysis_f,
                m_u_posttx_03_m = m_u_posttx_03_m,
                m_u_posttx_03_f = m_u_posttx_03_f,
                m_u_posttx_48_m = m_u_posttx_48_m,
                m_u_posttx_48_f = m_u_posttx_48_f,
                m_u_posttx_912_m = m_u_posttx_912_m,
                m_u_posttx_912_f = m_u_posttx_912_f,
                m_u_posttx_13_m = m_u_posttx_13_m,
                m_u_posttx_13_f = m_u_posttx_13_f,
                v_c_dialysis = v_c_dialysis,
                v_c_posttx = v_c_posttx,
                c_ddtx = c_ddtx,
                c_oacc = c_oacc,
                c_ldtx = c_ldtx,
                c_dgf = c_dgf,
                c_graftloss = c_graftloss,
                c_giver_dialysis = c_giver_dialysis,
                c_giver_posttx = c_giver_posttx,
                c_patient_dialysis = c_patient_dialysis,
                v_c_patient_posttx = v_c_patient_posttx,
                v_c_giver_posttx = v_c_giver_posttx,
                rr_removal = rr_removal))
}

f_costs_qalys_psa <- function(data, r, params, x) {
  # sum(z/((1+r)^(t-1)))
  # z: vector of costs
  # r: discount rate
  # t: vector of months

  v_u_dialysis_m <- params$m_u_dialysis_m[x,]
  v_u_dialysis_f <- params$m_u_dialysis_f[x,]
  v_u_posttx_03_m <- params$m_u_posttx_03_m[x,]
  v_u_posttx_03_f <- params$m_u_posttx_03_f[x,]
  v_u_posttx_48_m <- params$m_u_posttx_48_m[x,]
  v_u_posttx_48_f <- params$m_u_posttx_48_f[x,]
  v_u_posttx_912_m <- params$m_u_posttx_912_m[x,]
  v_u_posttx_912_f <- params$m_u_posttx_912_f[x,]
  v_u_posttx_13_m <- params$m_u_posttx_13_m[x,]
  v_u_posttx_13_f <- params$m_u_posttx_13_f[x,]
  v_c_dialysis <- params$v_c_dialysis[x]
  v_c_posttx <- params$v_c_posttx[x]
  c_ddtx <- params$c_ddtx[x]
  c_oacc <- params$c_oacc[x]
  c_ldtx <- params$c_ldtx[x]
  c_dgf <- params$c_dgf[x]
  c_graftloss <- params$c_graftloss[x]
  c_giver_dialysis <- params$c_giver_dialysis[x]
  c_giver_posttx <- params$c_giver_posttx[x]
  c_patient_dialysis <- params$c_patient_dialysis[x]
  v_c_patient_posttx <- params$v_c_patient_posttx[x,]
  rr_removal <- params$rr_removal[x]

  # if we use age-dependent costs
  # v_age <- seq(from = data$age_at_listing[i], to = 100, by = 1/12)
  # v_qaly_gf_1_12 <- c(rep(u_posttx_03, 2), rep(u_posttx_48, 5),
  # rep(u_posttx_912, 4))
  cost_health    <- rep(NULL, nrow(data))
  cost_societal  <- rep(NULL, nrow(data))
  qalys          <- rep(NULL, nrow(data))


  for(i in seq_len(nrow(data))) {
    male           <- data$male[i]
    age_list       <- data$can_age_at_listing_c[i] + 65
    age_list_event <- data$can_age_at_listing_c[i] + 65 +
                      (data$months_to_event[i] / 12)
    ## Waitlist Costs & QALYs
    l_cq <- get_costs_qalys_dial(age_dial = age_list,
                                 dialysis_time = data$months_to_event[i] - 1,
                                 sex = male,
                                 v_dial_cost = v_c_dialysis,
                                 v_dial_qaly_m = v_u_dialysis_m,
                                 v_dial_qaly_f = v_u_dialysis_f,
                                 c_patient_time = c_patient_dialysis,
                                 c_giver_time = c_giver_dialysis)
    v_costs_health   <- l_cq$v_cost_health
    v_costs_societal <- l_cq$v_cost_soc
    v_qalys          <- l_cq$v_qaly

    ## Deceased Donor Costs & QALYs
    # 30-Day Outcomes
    if(data$ddtx[i] == 1) {
      pos <- get_position_qaly_vector(age_list_event)
      if(data$tx_outcome[i] == 0) {  # If graft success
        v_costs_health   <- c(v_costs_health, c_ddtx + c_oacc)
        v_costs_societal <- c(v_costs_societal, c_ddtx + c_oacc)
        if (male == 1) {
          v_qalys        <- c(v_qalys, v_u_posttx_03_m[pos])
        } else {
          v_qalys        <- c(v_qalys, v_u_posttx_03_f[pos])
        }
      } else if(data$tx_outcome[i] == 1) {  # If graft failure
        v_costs_health   <- c(v_costs_health, c_ddtx + c_oacc + c_graftloss)
        v_costs_societal <- c(v_costs_societal, c_ddtx + c_oacc + c_graftloss)
        if (male == 1) {
          v_qalys        <- c(v_qalys, v_u_dialysis_m[pos])
        } else {
          v_qalys        <- c(v_qalys, v_u_dialysis_f[pos])
        }
      } else if(data$tx_outcome[i] == 2) {  # If delayed graft function
        v_costs_health   <- c(v_costs_health, c_ddtx + c_oacc + c_dgf)
        v_costs_societal <- c(v_costs_societal, c_ddtx + c_oacc + c_dgf)
        if (male == 1) {
          v_qalys        <- c(v_qalys, 0) # updated to be 0 for month of event
        } else {
          v_qalys        <- c(v_qalys, 0) # updated to be 0 for month of event
        }
      } else if(data$tx_outcome[i] == 3) {  # If death
        v_costs_health   <- c(v_costs_health, c_ddtx + c_oacc)
        v_costs_societal <- c(v_costs_societal, c_ddtx + c_oacc)
        v_qalys          <- c(v_qalys, 0)
      }

      # Post-Tx (after 30 days)
      l_post_cq <- get_costs_qalys_posttx(age_tx = age_list_event,
                                          gs_time = max(data$gl_time[i] - 1, 0),
                                          sex = male,
                                          graft_loss = data$graft_loss[i],
                                          gl_death_time = data$gl_death_time[i],
                                          v_posttx_cost = v_c_posttx,
                                          v_dial_cost = v_c_dialysis,
                                          v_dial_qaly_m = v_u_dialysis_m,
                                          v_dial_qaly_f = v_u_dialysis_f,
                                          c_patient_time = v_c_patient_posttx,
                                          c_giver_time = c_giver_posttx,
                                          c_patient_time_gl = c_patient_dialysis,
                                          c_giver_time_gl = c_giver_dialysis,
                                          c_graftloss = c_graftloss,
                                          v_u_posttx_03_m = v_u_posttx_03_m,
                                          v_u_posttx_48_m = v_u_posttx_48_m,
                                          v_u_posttx_912_m = v_u_posttx_912_m,
                                          v_u_posttx_13_m = v_u_posttx_13_m,
                                          v_u_posttx_03_f = v_u_posttx_03_f,
                                          v_u_posttx_48_f = v_u_posttx_48_f,
                                          v_u_posttx_912_f = v_u_posttx_912_f,
                                          v_u_posttx_13_f = v_u_posttx_13_f)
      v_costs_health   <- c(v_costs_health, l_post_cq$v_cost_health)
      v_costs_societal <- c(v_costs_societal, l_post_cq$v_cost_soc)
      v_qalys          <- c(v_qalys, l_post_cq$v_qaly)

    } else if(data$ldtx[i] == 1) {
      pos <- get_position_qaly_vector(age_list_event)
      le_ldtx <- max(round(rgamma(n=1, shape = (120^2)/1296,
                                  scale = (1296/120))), 1)
      l_post_cq <- get_costs_qalys_posttx(age_tx = age_list_event,
                                          gs_time = le_ldtx - 1,
                                          sex = male,
                                          graft_loss = 0,
                                          v_posttx_cost = v_c_posttx,
                                          v_dial_cost = v_c_dialysis,
                                          v_dial_qaly_m = v_u_dialysis_m,
                                          v_dial_qaly_f = v_u_dialysis_f,
                                          c_patient_time = v_c_patient_posttx,
                                          c_giver_time = c_giver_posttx,
                                          c_patient_time_gl = c_patient_dialysis,
                                          c_giver_time_gl = c_giver_dialysis,
                                          c_graftloss = c_graftloss,
                                          v_u_posttx_03_m = v_u_posttx_03_m,
                                          v_u_posttx_48_m = v_u_posttx_48_m,
                                          v_u_posttx_912_m = v_u_posttx_912_m,
                                          v_u_posttx_13_m = v_u_posttx_13_m,
                                          v_u_posttx_03_f = v_u_posttx_03_f,
                                          v_u_posttx_48_f = v_u_posttx_48_f,
                                          v_u_posttx_912_f = v_u_posttx_912_f,
                                          v_u_posttx_13_f = v_u_posttx_13_f)
      v_costs_health   <- c(v_costs_health, c_ldtx, l_post_cq$v_cost_health)
      v_costs_societal <- c(v_costs_societal, c_ldtx, l_post_cq$v_cost_soc)
      if (male == 1) {
        v_qalys        <- c(v_qalys, v_u_posttx_03_m[pos], l_post_cq$v_qaly)
      } else {
        v_qalys        <- c(v_qalys, v_u_posttx_03_f[pos], l_post_cq$v_qaly)
      }

    } else if(data$mort[i] == 1) {
      v_costs_health   <- c(v_costs_health, 0)
      v_costs_societal <- c(v_costs_societal, 0)
      v_qalys          <- c(v_qalys, 0)
    } else if(data$remove[i]==1) {
      le_remove <- max(round(rgamma(n=1, shape = (30^2)/900,
                                    scale = (900/30))), 1)
      l_cq <- get_costs_qalys_dial(age_dial = age_list,
                                   dialysis_time = le_remove - 1,
                                   sex = male,
                                   v_dial_cost = v_c_dialysis * rr_removal,
                                   v_dial_qaly_m = v_u_dialysis_m,
                                   v_dial_qaly_f = v_u_dialysis_f,
                                   c_patient_time = c_patient_dialysis,
                                   c_giver_time = c_giver_dialysis)

      v_costs_health   <- c(v_costs_health, l_cq$v_cost_health)
      v_costs_societal <- c(v_costs_societal, l_cq$v_cost_soc)
      v_qalys          <- c(v_qalys, l_cq$v_qaly)
    }

    v_costs_health   <- na.omit(v_costs_health)
    v_costs_societal <- na.omit(v_costs_societal)
    v_qalys          <- na.omit(v_qalys)

    v_qalys <- unlist(v_qalys)

    t_cost_health    <- seq_along(v_costs_health)
    t_costs_societal <- seq_along(v_costs_societal)
    t_qaly           <- seq_along(v_qalys)

    cost_health[i]   <- sum(v_costs_health/((1+r)^(t_cost_health-1)))
    cost_societal[i] <- sum(v_costs_societal/((1+r)^(t_costs_societal-1)))
    qalys[i]         <- sum(v_qalys/((1+r)^(t_qaly-1)))
  }


  ## Get total cost and QALYs
  tot_cost_health   <- sum(cost_health)
  tot_cost_societal <- sum(cost_societal)
  tot_qalys         <- sum(qalys)

  ## Get total cost and QALYs by subgroup
  # Diabetes
  cost_health_diab_0   <- cost_health[which(data$diab_stat == 0)]
  cost_societal_diab_0 <- cost_societal[which(data$diab_stat == 0)]
  qalys_diab_0         <- qalys[which(data$diab_stat == 0)]

  cost_health_diab_1   <- cost_health[which(data$diab_stat == 1)]
  cost_societal_diab_1 <- cost_societal[which(data$diab_stat == 1)]
  qalys_diab_1         <- qalys[which(data$diab_stat == 1)]

  tot_cost_health_diab_0    <- sum(cost_health_diab_0)
  tot_cost_societal_diab_0  <- sum(cost_societal_diab_0)
  tot_qalys_diab_0          <- sum(qalys_diab_0)

  tot_cost_health_diab_1    <- sum(cost_health_diab_1)
  tot_cost_societal_diab_1  <- sum(cost_societal_diab_1)
  tot_qalys_diab_1          <- sum(qalys_diab_1)


  # Race/Ethnicity
  cost_health_white   <- cost_health[which(data$race == 0)]
  cost_societal_white <- cost_societal[which(data$race == 0)]
  qalys_white         <- qalys[which(data$race == 0)]

  cost_health_black   <- cost_health[which(data$race == 1)]
  cost_societal_black <- cost_societal[which(data$race == 1)]
  qalys_black         <- qalys[which(data$race == 1)]

  cost_health_hispanic   <- cost_health[which(data$race == 2)]
  cost_societal_hispanic <- cost_societal[which(data$race == 2)]
  qalys_hispanic         <- qalys[which(data$race == 2)]

  cost_health_other   <- cost_health[which(data$race == 3)]
  cost_societal_other <- cost_societal[which(data$race == 3)]
  qalys_other         <- qalys[which(data$race == 3)]

  tot_cost_health_white    <- sum(cost_health_white)
  tot_cost_societal_white  <- sum(cost_societal_white)
  tot_qalys_white          <- sum(qalys_white)

  tot_cost_health_black    <- sum(cost_health_black)
  tot_cost_societal_black  <- sum(cost_societal_black)
  tot_qalys_black          <- sum(qalys_black)

  tot_cost_health_hispanic    <- sum(cost_health_hispanic)
  tot_cost_societal_hispanic  <- sum(cost_societal_hispanic)
  tot_qalys_hispanic          <- sum(qalys_hispanic)

  tot_cost_health_other    <- sum(cost_health_other)
  tot_cost_societal_other  <- sum(cost_societal_other)
  tot_qalys_other          <- sum(qalys_other)

  # Age at Listing
  cost_health_65   <- cost_health[which(data$can_age_at_listing_c < 5)]
  cost_societal_65 <- cost_societal[which(data$can_age_at_listing_c < 5)]
  qalys_65         <- qalys[which(data$can_age_at_listing_c < 5)]

  cost_health_70   <- cost_health[which(data$can_age_at_listing_c >= 5 &
                                        data$can_age_at_listing_c < 10)]
  cost_societal_70 <- cost_societal[which(data$can_age_at_listing_c >= 5 &
                                          data$can_age_at_listing_c < 10)]
  qalys_70         <- qalys[which(data$can_age_at_listing_c >= 5 &
                                  data$can_age_at_listing_c < 10)]

  cost_health_75   <- cost_health[which(data$can_age_at_listing_c >= 10)]
  cost_societal_75 <- cost_societal[which(data$can_age_at_listing_c >= 10)]
  qalys_75         <- qalys[which(data$can_age_at_listing_c >= 10)]

  tot_cost_health_65    <- sum(cost_health_65)
  tot_cost_societal_65  <- sum(cost_societal_65)
  tot_qalys_65          <- sum(qalys_65)

  tot_cost_health_70    <- sum(cost_health_70)
  tot_cost_societal_70  <- sum(cost_societal_70)
  tot_qalys_70          <- sum(qalys_70)

  tot_cost_health_75    <- sum(cost_health_75)
  tot_cost_societal_75  <- sum(cost_societal_75)
  tot_qalys_75          <- sum(qalys_75)

  ## Event counts by subgroup
  # All candidates
  ddtx_count      <- sum(data$ddtx == 1)
  ldtx_count      <- sum(data$ldtx == 1)
  wait_mort_count <- sum(data$mort == 1)
  removal_count   <- sum(data$remove == 1)
  dgf_count       <- sum(data$tx_outcome == 2, na.rm = T)
  gl_count        <- sum(data$graft_loss == 1, na.rm = T)

  error_count <- function(outcome, outcome_val, subgroup, subgroup_val) {
    tryCatch(expr = {
      tmp_count <- sum(outcome == outcome_val & subgroup == subgroup_val, na.rm = T)
    }, error = function(e) {tmp_count <- NaN}
    )
    return(tmp_count)
  }
  # By diabetes status
  ddtx_count_diab0      <- error_count(data$ddtx, 1, data$diab_stat, 0)
  ldtx_count_diab0      <- error_count(data$ldtx, 1, data$diab_stat, 0)
  wait_mort_count_diab0 <- error_count(data$mort, 1, data$diab_stat, 0)
  removal_count_diab0   <- error_count(data$remove, 1, data$diab_stat, 0)
  dgf_count_diab0       <- error_count(data$tx_outcome, 2, data$diab_stat, 0)
  gl_count_diab0        <- error_count(data$graft_loss, 1, data$diab_stat, 0)

  ddtx_count_diab1      <- error_count(data$ddtx, 1, data$diab_stat, 1)
  ldtx_count_diab1      <- error_count(data$ldtx, 1, data$diab_stat, 1)
  wait_mort_count_diab1 <- error_count(data$mort, 1, data$diab_stat, 1)
  removal_count_diab1   <- error_count(data$remove, 1, data$diab_stat, 1)
  dgf_count_diab1       <- error_count(data$tx_outcome, 2, data$diab_stat, 1)
  gl_count_diab1        <- error_count(data$graft_loss, 1, data$diab_stat, 1)

  # By race/ethnicity
  ddtx_count_white      <- error_count(data$ddtx, 1, data$race, 0)
  ldtx_count_white      <- error_count(data$ldtx, 1, data$race, 0)
  wait_mort_count_white <- error_count(data$mort, 1, data$race, 0)
  removal_count_white   <- error_count(data$remove, 1, data$race, 0)
  dgf_count_white       <- error_count(data$tx_outcome, 2, data$race, 0)
  gl_count_white        <- error_count(data$graft_loss, 1, data$race, 0)

  ddtx_count_black      <- error_count(data$ddtx, 1, data$race, 1)
  ldtx_count_black      <- error_count(data$ldtx, 1, data$race, 1)
  wait_mort_count_black <- error_count(data$mort, 1, data$race, 1)
  removal_count_black   <- error_count(data$remove, 1, data$race, 1)
  dgf_count_black       <- error_count(data$tx_outcome, 2, data$race, 1)
  gl_count_black        <- error_count(data$graft_loss, 1, data$race, 1)

  ddtx_count_hisp      <- error_count(data$ddtx, 1, data$race, 2)
  ldtx_count_hisp      <- error_count(data$ldtx, 1, data$race, 2)
  wait_mort_count_hisp <- error_count(data$mort, 1, data$race, 2)
  removal_count_hisp   <- error_count(data$remove, 1, data$race, 2)
  dgf_count_hisp       <- error_count(data$tx_outcome, 2, data$race, 2)
  gl_count_hisp        <- error_count(data$graft_loss, 1, data$race, 2)

  ddtx_count_other      <- error_count(data$ddtx, 1, data$race, 3)
  ldtx_count_other      <- error_count(data$ldtx, 1, data$race, 3)
  wait_mort_count_other <- error_count(data$mort, 1, data$race, 3)
  removal_count_other   <- error_count(data$remove, 1, data$race, 3)
  dgf_count_other       <- error_count(data$tx_outcome, 2, data$race, 3)
  gl_count_other        <- error_count(data$graft_loss, 1, data$race, 3)

  # By age at listing
  data$age_grp <- ifelse(data$can_age_at_listing_c < 5, 1,
                          ifelse(data$can_age_at_listing_c >= 5 & data$can_age_at_listing_c < 10, 2, 3))

  ddtx_count_65      <- error_count(data$ddtx, 1, data$age_grp, 1)
  ldtx_count_65      <- error_count(data$ldtx, 1, data$age_grp, 1)
  wait_mort_count_65 <- error_count(data$mort, 1, data$age_grp, 1)
  removal_count_65   <- error_count(data$remove, 1, data$age_grp, 1)
  dgf_count_65       <- error_count(data$tx_outcome, 2, data$age_grp, 1)
  gl_count_65        <- error_count(data$graft_loss, 1, data$age_grp, 1)

  ddtx_count_70      <- error_count(data$ddtx, 1, data$age_grp, 2)
  ldtx_count_70      <- error_count(data$ldtx, 1, data$age_grp, 2)
  wait_mort_count_70 <- error_count(data$mort, 1, data$age_grp, 2)
  removal_count_70   <- error_count(data$remove, 1, data$age_grp, 2)
  dgf_count_70       <- error_count(data$tx_outcome, 2, data$age_grp, 2)
  gl_count_70        <- error_count(data$graft_loss, 1, data$age_grp, 2)

  ddtx_count_75      <- error_count(data$ddtx, 1, data$age_grp, 3)
  ldtx_count_75      <- error_count(data$ldtx, 1, data$age_grp, 3)
  wait_mort_count_75 <- error_count(data$mort, 1, data$age_grp, 3)
  removal_count_75   <- error_count(data$remove, 1, data$age_grp, 3)
  dgf_count_75       <- error_count(data$tx_outcome, 2, data$age_grp, 3)
  gl_count_75        <- error_count(data$graft_loss, 1, data$age_grp, 3)

  event_counts <- c(ddtx_count, ldtx_count, wait_mort_count,
                    removal_count, dgf_count, gl_count,
                    ddtx_count_diab0, ldtx_count_diab0, wait_mort_count_diab0,
                    removal_count_diab0, dgf_count_diab0, gl_count_diab0,
                    ddtx_count_diab1, ldtx_count_diab1, wait_mort_count_diab1,
                    removal_count_diab1, dgf_count_diab1, gl_count_diab1,
                    ddtx_count_white, ldtx_count_white, wait_mort_count_white,
                    removal_count_white, dgf_count_white, gl_count_white,
                    ddtx_count_black, ldtx_count_black, wait_mort_count_black,
                    removal_count_black, dgf_count_black, gl_count_black,
                    ddtx_count_hisp, ldtx_count_hisp, wait_mort_count_hisp,
                    removal_count_hisp, dgf_count_hisp, gl_count_hisp,
                    ddtx_count_other, ldtx_count_other, wait_mort_count_other,
                    removal_count_other, dgf_count_other, gl_count_other,
                    ddtx_count_65, ldtx_count_65, wait_mort_count_65,
                    removal_count_65, dgf_count_65, gl_count_65,
                    ddtx_count_70, ldtx_count_70, wait_mort_count_70,
                    removal_count_70, dgf_count_70, gl_count_70,
                    ddtx_count_75, ldtx_count_75, wait_mort_count_75,
                    removal_count_75, dgf_count_75, gl_count_75) # dgf_count_75, gl_count_75

  cost_qalys <- c(tot_cost_health, tot_cost_societal, tot_qalys, tot_cost_health_diab_0,
                  tot_cost_societal_diab_0, tot_qalys_diab_0, tot_cost_health_diab_1,
                  tot_cost_societal_diab_1, tot_qalys_diab_1, tot_cost_health_white,
                  tot_cost_societal_white, tot_qalys_white, tot_cost_health_black,
                  tot_cost_societal_black, tot_qalys_black, tot_cost_health_hispanic,
                  tot_cost_societal_hispanic, tot_qalys_hispanic, tot_cost_health_other,
                  tot_cost_societal_other, tot_qalys_other, tot_cost_health_65,
                  tot_cost_societal_65, tot_qalys_65, tot_cost_health_70,
                  tot_cost_societal_70, tot_qalys_70, tot_cost_health_75,
                  tot_cost_societal_75, tot_qalys_75, event_counts)
  return(cost_qalys)
}

generate_sse <- function(beta, mu_target, sd_target) {
  alpha <- mu_target * beta
  v_samples <- rgamma(n = 100000, shape = alpha , rate = beta)
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
                        lower = 0,
                        upper = 1)
  beta  <- optim_result$par
  alpha <- mu * beta

  return(list(alpha = alpha, beta = beta))
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
