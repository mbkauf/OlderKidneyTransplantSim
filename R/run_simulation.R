#' Run the simulation
#'
#' @param seed random seed
#' @param data data of simulated individuals
#' @param ddtx_coef deceased donor transplant coefficients
#' @param ddtx_p deceased donor transplant shape parameter
#' @param ldtx_coef living donor transplant coefficients
#' @param ldtx_gamma living donor transplant shape parameter
#' @param mort_coef waitlist mortality coefficients
#' @param mort_shape waitlist mortality shape parameter
#' @param remove_coef waitlist removal coefficients
#' @param remove_p waitlist removal shape parameter
#' @param ddtx_rate increased rate of transplantation?
#' @param france use KDPI distribution from France?
#' @param gl30_coef 30-Day Multinomial Logit: graft loss coefficients
#' @param dgf_coef 30-Day Multinomial Logit: delayed graft function coefficients
#' @param mort30_coef 30-Day Multinomial Logit: mortality coefficients
#' @param gl_coef graft loss coefficients
#' @param gl_shape graft loss shape parameter
#' @param gs_mort_coef death with function coefficients
#' @param gs_mort_shape death with function shape parameter
#' @param dial_mort_coef death after graft loss coefficients
#' @param dial_mort_shape death after graft loss shape parameter
#' @param gl_mort_coef death same day as graft loss coefficients
#'
#' @return data frame with results for each individual
#' @export
#'
#' @examples
#' results <- run.simulation(seed = 11111)
run_simulation  <- function(seed, data = sample_df,
                            ddtx_coef, ddtx_p,
                            ldtx_coef, ldtx_gamma,
                            mort_coef, mort_shape,
                            remove_coef, remove_p,
                            ddtx_rate = 0, france = FALSE,
                            gl30_coef, dgf_coef, mort30_coef,
                            gl_coef, gl_shape,
                            gs_mort_coef, gs_mort_shape,
                            dial_mort_coef, dial_mort_shape,
                            gl_mort_coef) {
  ### Function Inputs
  # seed: used to replicate results
  # data: dataset with simulated patient cohort
  # ddtx_coef & ddtx_p: coefficient and shape parameters for deceased donor tx
  # ldtx_coef & ldtx_gamma: coefficient and shape parameters for living donor tx
  # mort_coef & mort_shape: coefficient and shape parameters for waitlist death
  # remove_coef & remove_p: coefficient and shape parameters for waitlist removal
  # ddtx_rate: how much to increase rate of transplantation
  # gl30_coef: coefficients for 30-day graft loss
  # dgf_coef: coefficients for delayed graft function
  # mort30_coef: coefficients for 30-day mortality
  # gl_coef & gl_shape: coefficient and shape parameters for graft loss
  # gs_mort_coef & gs_mort_shape: coefficient and shape parameters for death w/ function
  # dial_mort_coef & dial_mort_shape: coefficient and shape parameters for death after graft loss
  # gl_mort_coef: coefficients for death at graft loss
  # france: use kdpi distribution of french allocation system?


  # Set up simulation
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

  ### Run Simulation
  i <- 0
  while (all(event_results == 1) == FALSE) {
    i = i + 1
    ### Deceased Donor Tx?
    # Generate vector of ddtx probabilities
    ddtx_haz  <- (exp(ddtx_p) * ddtx_lambda * i^(exp(ddtx_p) - 1))*(1+ddtx_rate)
    ddtx_prob <- 1 - exp(-ddtx_haz)

    # Determine who has event by using random uniform draws and determine who
    # had a new event this period
    rand_ddtx_surv <- runif(num_patients, min = 0, max = 1)  # Vector of random probabilities from uniform
    y              <- as.integer(rand_ddtx_surv < ddtx_prob)  # vector of those who had event in period i
    new_ddtx       <- as.integer(y==1 & pmax(y, ddtx_results) > ddtx_results)  # vector of new events

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
                                 pmax(ddtx_results, mort_results, remove_results)==1)
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
    new_mort       <- as.integer(x==1 & pmax(x, mort_results) > mort_results)

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

    # Create matrix that tracks cumulative events, create patient trace
    # remove_trace[,i] <- remove_results

    ### Time to any event
    event_results <- pmax(ddtx_results, ldtx_results,
                          mort_results, remove_results)
    m2event       <- m2event + pmax(new_ddtx, new_ldtx,
                                    new_mort, new_remove)*i
  }

  results_sim <- as.data.frame(cbind(data, event_results, ddtx_results,
                                     ldtx_results, mort_results, remove_results,
                                     m2event))
  results_sim <- results_sim %>%
    rename(
      event = event_results,
      ddtx = ddtx_results,
      ldtx = ldtx_results,
      mort = mort_results,
      remove = remove_results,
      months_to_event = m2event
    )

  # Prepare dataset for post-transplant simulation
  no_ddtx <- results_sim %>%
    dplyr::filter(ddtx==0)

  posttx_df <- results_sim %>%
    dplyr::filter(ddtx==1)

  colvars = names(posttx_df)
  colvars = colnames(posttx_df)
  start_loc = match("optn_reg2",colvars)
  end_loc = match("optn_reg11",colvars)
  optn <- posttx_df[,start_loc:end_loc]


  if(france==FALSE){
    Z <- invisible(rand_kdpi_cold(n = nrow(posttx_df)))
  } else {
    Z <- invisible(rand_kdpi_cold(n = nrow(posttx_df), france=TRUE))
  }

  posttx_df <- cbind(posttx_df, Z)
  names(posttx_df)[ncol(posttx_df)] <- "kdpi"
  names(posttx_df)[ncol(posttx_df)-1] <- "rec_cold_isch_tm"

  # data <- as.data.frame(data)

  ddtx_df <- posttx_df

  posttx_df <- posttx_df %>%
    dplyr::mutate(rec_age_at_tx_c = can_age_at_listing_c + (months_to_event/12)) %>%
    dplyr::mutate(years_dial = baseline_yrs_dial + (months_to_event/12)) %>%
    dplyr::mutate(tx_year_c = list_year_c + (months_to_event/12)) %>%
    dplyr::mutate(tx_year_c = ifelse(tx_year_c > 10, 10, tx_year_c)) %>%
    dplyr::mutate(constant = 1) %>%
    dplyr::select(rec_age_at_tx_c,	male,	black, hispanic, other_race,
                  years_dial, canhx_cpra, kdpi, can_blood_ab, can_blood_b,
                  can_blood_o, rec_cold_isch_tm, tx_year_c, diab_stat,
                  copd, pvd,	ang_cad, constant)

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


  data_gl_mort <- data_gl %>%
    dplyr::mutate(m2gl = rep(0, nrow(data_gl))) %>%
    dplyr::select(rec_age_at_tx_c,	male,	black, hispanic, other_race,
                  years_dial, canhx_cpra, kdpi, can_blood_ab, can_blood_b,
                  can_blood_o, rec_cold_isch_tm, tx_year_c, diab_stat, copd,
                  pvd, ang_cad, dgf_results, m2gl, constant)
  data_gs_mort <- data_gl %>%
    dplyr::select(rec_age_at_tx_c,	male,	black, hispanic, other_race,
                  years_dial, canhx_cpra, kdpi, can_blood_ab, can_blood_b,
                  can_blood_o, rec_cold_isch_tm, diab_stat, copd,
                  pvd, ang_cad, dgf_results, constant) # tx_year_c
  data_gl <- data_gl %>%
    dplyr::select(rec_age_at_tx_c,	male,	black, hispanic, other_race,
                  years_dial, canhx_cpra, kdpi, can_blood_ab, can_blood_b,
                  can_blood_o, rec_cold_isch_tm, tx_year_c, diab_stat, copd,
                  pvd, ang_cad, dgf_results, constant)

  data_gl_mort <- as.matrix(data_gl_mort)
  data_gl      <- as.matrix(data_gl)
  data_gs_mort <- as.matrix(data_gs_mort)

  # build.data(n=num_patients)
  ### Graft loss (PH Gompertz)
  gl_lambda <- exp(data_gl %*% gl_coef)

  ### Death w/ functioning graft (PH Gompertz)
  gs_mort_lambda <- exp(data_gs_mort %*% gs_mort_coef)

  ### Initialize vectors and matrices
  gs_mort_results <- rep(0, num_patients)
  # gl_results      <- max(rep(0, num_patients), gl30_results)
  gl_mort_results <- rep(0, num_patients)

  m2gl      <- rep(0, num_patients)
  m2gs_mort <- rep(0, num_patients)
  m2gl_mort <- rep(0, num_patients)

  death    <- rep(0, num_patients)
  print(num_patients)
  print(length(mort30_results))
  m2age100 <- (rep(35, num_patients) - data_gs_mort[,1]) * 12

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

    ### Next we determine who had a graft loss event
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

  list_vars_tx = c("gl_results", "gs_mort_results", "gl_mort_results", "m2gl",
                   "m2gs_mort", "m2gl_mort", "dgf_results", "tx_outcome",
                   "kdpi", "rec_cold_isch_tm")

  no_ddtx <- cbind(no_ddtx, setNames(lapply(list_vars_tx,
                                            function(x) x=NA), list_vars_tx))

  sim_df <- as.data.frame(cbind(ddtx_df, gl_results, gs_mort_results,
                                gl_mort_results, m2gl, m2gs_mort, m2gl_mort,
                                dgf_results, tx_outcome))
  all_df <- rbind(sim_df, no_ddtx)

  return(all_df)
}
