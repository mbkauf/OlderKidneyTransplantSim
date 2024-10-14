library(flexsurv)
library(survival)
library(dplyr)
library(mvtnorm)
library(haven)
library(MCMCpack)
library(LaplacesDemon)


get_likelihood_targets <- function(data_list, data_tx) {

  ## Set up data for analysis
  ### Test some code
  data <- data_list %>%
    mutate(days_to_65 = round((65 - can_age_at_listing - 0.5) * 365)) %>%
    mutate(approx_65_dt = can_listing_dt + days_to_65) %>%
    mutate(start_dt = can_listing_dt) %>%
    mutate(start_dt = if_else(can_age_at_listing < 65, approx_65_dt, start_dt)) %>%
    mutate(t_start = as.numeric(targetsdifftime(as.Date(start_dt), as.Date(can_listing_dt), units = "days")) / 30.4375) %>%
    mutate(censor_time = as.numeric(targetsdifftime(as.Date(can_rem_dt), as.Date(start_dt), units = "days")) / 30.4375) %>%
    mutate(censor_time = ifelse(is.na(censor_time), as.numeric(targetsdifftime(as.Date(can_endwlfu), as.Date(start_dt), units = "days")) / 30.4375, censor_time)) %>%
    mutate(t = censor_time + t_start) %>%
    mutate(t0 = t_start) %>%
    filter(test==1)

  ddtx_surv <- matrix(nrow = 2, ncol = length(seq(from=3, to=180, by=3)))
  ldtx_surv <- matrix(nrow = 2, ncol = length(seq(from=3, to=180, by=3)))
  mort_surv <- matrix(nrow = 2, ncol = length(seq(from=3, to=180, by=3)))
  remove_surv <- matrix(nrow = 2, ncol = length(seq(from=3, to=180, by=3)))
  ddtx_surv_diab0 <- matrix(nrow = 2, ncol = length(seq(from=3, to=180, by=3)))
  ldtx_surv_diab0 <- matrix(nrow = 2, ncol = length(seq(from=3, to=180, by=3)))
  mort_surv_diab0 <- matrix(nrow = 2, ncol = length(seq(from=3, to=180, by=3)))
  remove_surv_diab0 <- matrix(nrow = 2, ncol = length(seq(from=3, to=180, by=3)))
  ddtx_surv_diab1 <- matrix(nrow = 2, ncol = length(seq(from=3, to=180, by=3)))
  ldtx_surv_diab1 <- matrix(nrow = 2, ncol = length(seq(from=3, to=180, by=3)))
  mort_surv_diab1 <- matrix(nrow = 2, ncol = length(seq(from=3, to=180, by=3)))
  remove_surv_diab1 <- matrix(nrow = 2, ncol = length(seq(from=3, to=180, by=3)))
  ddtx_surv_white <- matrix(nrow = 2, ncol = length(seq(from=3, to=180, by=3)))
  ldtx_surv_white <- matrix(nrow = 2, ncol = length(seq(from=3, to=180, by=3)))
  mort_surv_white <- matrix(nrow = 2, ncol = length(seq(from=3, to=180, by=3)))
  remove_surv_white <- matrix(nrow = 2, ncol = length(seq(from=3, to=180, by=3)))
  ddtx_surv_black <- matrix(nrow = 2, ncol = length(seq(from=3, to=180, by=3)))
  ldtx_surv_black <- matrix(nrow = 2, ncol = length(seq(from=3, to=180, by=3)))
  mort_surv_black <- matrix(nrow = 2, ncol = length(seq(from=3, to=180, by=3)))
  remove_surv_black <- matrix(nrow = 2, ncol = length(seq(from=3, to=180, by=3)))
  ddtx_surv_hispanic <- matrix(nrow = 2, ncol = length(seq(from=3, to=180, by=3)))
  ldtx_surv_hispanic <- matrix(nrow = 2, ncol = length(seq(from=3, to=180, by=3)))
  mort_surv_hispanic <- matrix(nrow = 2, ncol = length(seq(from=3, to=180, by=3)))
  remove_surv_hispanic <- matrix(nrow = 2, ncol = length(seq(from=3, to=180, by=3)))
  ddtx_surv_other <- matrix(nrow = 2, ncol = length(seq(from=3, to=180, by=3)))
  ldtx_surv_other <- matrix(nrow = 2, ncol = length(seq(from=3, to=180, by=3)))
  mort_surv_other <- matrix(nrow = 2, ncol = length(seq(from=3, to=180, by=3)))
  remove_surv_other <- matrix(nrow = 2, ncol = length(seq(from=3, to=180, by=3)))

  ### All Candidates
  # Kaplan-Meier
  ddtx_km_fit    <- survival::survfit(Surv(time = t0, time2 = t, event_cd==0) ~ 1, data = data)
  ldtx_km_fit    <- survival::survfit(Surv(time = t0, time2 = t, event_cd==1) ~ 1, data = data)
  mort_km_fit    <- survival::survfit(Surv(time = t0, time2 = t, event_cd==2) ~ 1, data = data)
  remove_km_fit  <- survival::survfit(Surv(time = t0, time2 = t, event_cd==3) ~ 1, data = data)

  ddtx_surv[1,]   <- t(summary(ddtx_km_fit, times = seq(from=3, to=180, by=3), extend=TRUE)$surv)
  ldtx_surv[1,]   <- t(summary(ldtx_km_fit, times = seq(from=3, to=180, by=3), extend=TRUE)$surv)
  mort_surv[1,]   <- t(summary(mort_km_fit, times = seq(from=3, to=180, by=3), extend=TRUE)$surv)
  remove_surv[1,] <- t(summary(remove_km_fit, times = seq(from=3, to=180, by=3), extend=TRUE)$surv)

  ddtx_surv[2,]   <- t(summary(ddtx_km_fit, times = seq(from=3, to=180, by=3), extend=TRUE)$std.err)
  ldtx_surv[2,]   <- t(summary(ldtx_km_fit, times = seq(from=3, to=180, by=3), extend=TRUE)$std.err)
  mort_surv[2,]   <- t(summary(mort_km_fit, times = seq(from=3, to=180, by=3), extend=TRUE)$std.err)
  remove_surv[2,] <- t(summary(remove_km_fit, times = seq(from=3, to=180, by=3), extend=TRUE)$std.err)

  # Parametric Fits
  ddtx_weibull_fit = flexsurv::flexsurvreg(Surv(time = t0, time2 = t, event_cd==0)~1,
                                 data=data, dist="weibull")
  ldtx_gompertz_fit = flexsurv::flexsurvreg(Surv(time = t0, time2 = t, event_cd==1)~1,
                                  data=data, dist="gompertz")
  mort_weibull_fit = flexsurv::flexsurvreg(Surv(time = t0, time2 = t, event_cd==2)~1,
                                 data=data, dist="weibull")
  remove_weibull_fit = flexsurv::flexsurvreg(Surv(time = t0, time2 = t, event_cd==3)~1,
                                   data=data, dist="weibull")

  ddtx_shape_all <- ddtx_weibull_fit$coefficients
  ldtx_shape_all <- ldtx_gompertz_fit$coefficients
  mort_shape_all <- mort_weibull_fit$coefficients
  remove_shape_all <- remove_weibull_fit$coefficients
  ddtx_shape_all_cov <- ddtx_weibull_fit$cov
  ldtx_shape_all_cov <- ldtx_gompertz_fit$cov
  mort_shape_all_cov <- mort_weibull_fit$cov
  remove_shape_all_cov <- remove_weibull_fit$cov

  ### By Diabetes Status
  df_diab <- split(data, data$diab_stat)

  ## No Diabetes
  # Kaplan-Meier
  ddtx_km_fit_diab0    <- survival::survfit(Surv(time = t0, time2 = t, event_cd==0) ~ 1, data = df_diab[[1]])
  ldtx_km_fit_diab0    <- survival::survfit(Surv(time = t0, time2 = t, event_cd==1) ~ 1, data = df_diab[[1]])
  mort_km_fit_diab0    <- survival::survfit(Surv(time = t0, time2 = t, event_cd==2) ~ 1, data = df_diab[[1]])
  remove_km_fit_diab0  <- survival::survfit(Surv(time = t0, time2 = t, event_cd==3) ~ 1, data = df_diab[[1]])

  ddtx_surv_diab0[1,]   <- t(summary(ddtx_km_fit_diab0, times = seq(from=3, to=180, by=3), extend=TRUE)$surv)
  ldtx_surv_diab0[1,]   <- t(summary(ldtx_km_fit_diab0, times = seq(from=3, to=180, by=3), extend=TRUE)$surv)
  mort_surv_diab0[1,]   <- t(summary(mort_km_fit_diab0, times = seq(from=3, to=180, by=3), extend=TRUE)$surv)
  remove_surv_diab0[1,] <- t(summary(remove_km_fit_diab0, times = seq(from=3, to=180, by=3), extend=TRUE)$surv)

  ddtx_surv_diab0[2,]   <- t(summary(ddtx_km_fit_diab0, times = seq(from=3, to=180, by=3), extend=TRUE)$std.err)
  ldtx_surv_diab0[2,]   <- t(summary(ldtx_km_fit_diab0, times = seq(from=3, to=180, by=3), extend=TRUE)$std.err)
  mort_surv_diab0[2,]   <- t(summary(mort_km_fit_diab0, times = seq(from=3, to=180, by=3), extend=TRUE)$std.err)
  remove_surv_diab0[2,] <- t(summary(remove_km_fit_diab0, times = seq(from=3, to=180, by=3), extend=TRUE)$std.err)

  # Parametric Fits
  ddtx_weibull_fit_diab0 = flexsurv::flexsurvreg(Surv(time = t0, time2 = t, event_cd==0)~1,
                                       data=df_diab[[1]], dist="weibull")
  ldtx_gompertz_fit_diab0 = flexsurv::flexsurvreg(Surv(time = t0, time2 = t, event_cd==1)~1,
                                        data=df_diab[[1]], dist="gompertz")
  mort_weibull_fit_diab0 = flexsurv::flexsurvreg(Surv(time = t0, time2 = t, event_cd==2)~1,
                                       data=df_diab[[1]], dist="weibull")
  remove_weibull_fit_diab0 = flexsurv::flexsurvreg(Surv(time = t0, time2 = t, event_cd==3)~1,
                                         data=df_diab[[1]], dist="weibull")

  ddtx_shape_diab0 <- ddtx_weibull_fit_diab0$coefficients
  ldtx_shape_diab0 <- ldtx_gompertz_fit_diab0$coefficients
  mort_shape_diab0 <- mort_weibull_fit_diab0$coefficients
  remove_shape_diab0 <- remove_weibull_fit_diab0$coefficients
  ddtx_shape_diab0_cov <- ddtx_weibull_fit_diab0$cov
  ldtx_shape_diab0_cov <- ldtx_gompertz_fit_diab0$cov
  mort_shape_diab0_cov <- mort_weibull_fit_diab0$cov
  remove_shape_diab0_cov <- remove_weibull_fit_diab0$cov

  ## History of Diabetes
  # Kaplan-Meier
  ddtx_km_fit_diab1    <- survival::survfit(Surv(time = t0, time2 = t, event_cd==0) ~ 1, data = df_diab[[2]])
  ldtx_km_fit_diab1    <- survival::survfit(Surv(time = t0, time2 = t, event_cd==1) ~ 1, data = df_diab[[2]])
  mort_km_fit_diab1    <- survival::survfit(Surv(time = t0, time2 = t, event_cd==2) ~ 1, data = df_diab[[2]])
  remove_km_fit_diab1  <- survival::survfit(Surv(time = t0, time2 = t, event_cd==3) ~ 1, data = df_diab[[2]])

  ddtx_surv_diab1[1,]   <- t(summary(ddtx_km_fit_diab1, times = seq(from=3, to=180, by=3), extend=TRUE)$surv)
  ldtx_surv_diab1[1,]   <- t(summary(ldtx_km_fit_diab1, times = seq(from=3, to=180, by=3), extend=TRUE)$surv)
  mort_surv_diab1[1,]   <- t(summary(mort_km_fit_diab1, times = seq(from=3, to=180, by=3), extend=TRUE)$surv)
  remove_surv_diab1[1,] <- t(summary(remove_km_fit_diab1, times = seq(from=3, to=180, by=3), extend=TRUE)$surv)

  ddtx_surv_diab1[2,]   <- t(summary(ddtx_km_fit_diab1, times = seq(from=3, to=180, by=3), extend=TRUE)$std.err)
  ldtx_surv_diab1[2,]   <- t(summary(ldtx_km_fit_diab1, times = seq(from=3, to=180, by=3), extend=TRUE)$std.err)
  mort_surv_diab1[2,]   <- t(summary(mort_km_fit_diab1, times = seq(from=3, to=180, by=3), extend=TRUE)$std.err)
  remove_surv_diab1[2,] <- t(summary(remove_km_fit_diab1, times = seq(from=3, to=180, by=3), extend=TRUE)$std.err)

  # Parametric Fits
  ddtx_weibull_fit_diab1 = flexsurv::flexsurvreg(Surv(time = t0, time2 = t, event_cd==0)~1,
                                       data=df_diab[[2]], dist="weibull")
  ldtx_gompertz_fit_diab1 = flexsurv::flexsurvreg(Surv(time = t0, time2 = t, event_cd==1)~1,
                                        data=df_diab[[2]], dist="gompertz")
  mort_weibull_fit_diab1 = flexsurv::flexsurvreg(Surv(time = t0, time2 = t, event_cd==2)~1,
                                       data=df_diab[[2]], dist="weibull")
  remove_weibull_fit_diab1 = flexsurv::flexsurvreg(Surv(time = t0, time2 = t, event_cd==3)~1,
                                         data=df_diab[[2]], dist="weibull")

  ddtx_shape_diab1 <- ddtx_weibull_fit_diab1$coefficients
  ldtx_shape_diab1 <- ldtx_gompertz_fit_diab1$coefficients
  mort_shape_diab1 <- mort_weibull_fit_diab1$coefficients
  remove_shape_diab1 <- remove_weibull_fit_diab1$coefficients
  ddtx_shape_diab1_cov <- ddtx_weibull_fit_diab1$cov
  ldtx_shape_diab1_cov <- ldtx_gompertz_fit_diab1$cov
  mort_shape_diab1_cov <- mort_weibull_fit_diab1$cov
  remove_shape_diab1_cov <- remove_weibull_fit_diab1$cov

  ### By Race/Ethnicity
  df_race <- split(data, data$race_ethnic)

  ## White Candidates
  # Kaplan-Meier
  ddtx_km_fit_white    <- survival::survfit(Surv(time = t0, time2 = t, event_cd==0) ~ 1, data = df_race[[1]])
  ldtx_km_fit_white    <- survival::survfit(Surv(time = t0, time2 = t, event_cd==1) ~ 1, data = df_race[[1]])
  mort_km_fit_white    <- survival::survfit(Surv(time = t0, time2 = t, event_cd==2) ~ 1, data = df_race[[1]])
  remove_km_fit_white  <- survival::survfit(Surv(time = t0, time2 = t, event_cd==3) ~ 1, data = df_race[[1]])

  ddtx_surv_white[1,]   <- t(summary(ddtx_km_fit_white, times = seq(from=3, to=180, by=3), extend=TRUE)$surv)
  ldtx_surv_white[1,]   <- t(summary(ldtx_km_fit_white, times = seq(from=3, to=180, by=3), extend=TRUE)$surv)
  mort_surv_white[1,]   <- t(summary(mort_km_fit_white, times = seq(from=3, to=180, by=3), extend=TRUE)$surv)
  remove_surv_white[1,] <- t(summary(remove_km_fit_white, times = seq(from=3, to=180, by=3), extend=TRUE)$surv)

  ddtx_surv_white[2,]   <- t(summary(ddtx_km_fit_white, times = seq(from=3, to=180, by=3), extend=TRUE)$std.err)
  ldtx_surv_white[2,]   <- t(summary(ldtx_km_fit_white, times = seq(from=3, to=180, by=3), extend=TRUE)$std.err)
  mort_surv_white[2,]   <- t(summary(mort_km_fit_white, times = seq(from=3, to=180, by=3), extend=TRUE)$std.err)
  remove_surv_white[2,] <- t(summary(remove_km_fit_white, times = seq(from=3, to=180, by=3), extend=TRUE)$std.err)

  # Parametric Fits
  ddtx_weibull_fit_white = flexsurv::flexsurvreg(Surv(time = t0, time2 = t, event_cd==0)~1,
                                       data=df_race[[1]], dist="weibull")
  ldtx_gompertz_fit_white = flexsurv::flexsurvreg(Surv(time = t0, time2 = t, event_cd==1)~1,
                                        data=df_race[[1]], dist="gompertz")
  mort_weibull_fit_white = flexsurv::flexsurvreg(Surv(time = t0, time2 = t, event_cd==2)~1,
                                       data=df_race[[1]], dist="weibull")
  remove_weibull_fit_white = flexsurv::flexsurvreg(Surv(time = t0, time2 = t, event_cd==3)~1,
                                         data=df_race[[1]], dist="weibull")

  ddtx_shape_white <- ddtx_weibull_fit_white$coefficients
  ldtx_shape_white <- ldtx_gompertz_fit_white$coefficients
  mort_shape_white <- mort_weibull_fit_white$coefficients
  remove_shape_white <- remove_weibull_fit_white$coefficients
  ddtx_shape_white_cov <- ddtx_weibull_fit_white$cov
  ldtx_shape_white_cov <- ldtx_gompertz_fit_white$cov
  mort_shape_white_cov <- mort_weibull_fit_white$cov
  remove_shape_white_cov <- remove_weibull_fit_white$cov

  ## Black Candidates
  # Kaplan-Meier
  ddtx_km_fit_black    <- survival::survfit(Surv(time = t0, time2 = t, event_cd==0) ~ 1, data = df_race[[2]])
  ldtx_km_fit_black    <- survival::survfit(Surv(time = t0, time2 = t, event_cd==1) ~ 1, data = df_race[[2]])
  mort_km_fit_black    <- survival::survfit(Surv(time = t0, time2 = t, event_cd==2) ~ 1, data = df_race[[2]])
  remove_km_fit_black  <- survival::survfit(Surv(time = t0, time2 = t, event_cd==3) ~ 1, data = df_race[[2]])

  ddtx_surv_black[1,]   <- t(summary(ddtx_km_fit_black, times = seq(from=3, to=180, by=3), extend=TRUE)$surv)
  ldtx_surv_black[1,]   <- t(summary(ldtx_km_fit_black, times = seq(from=3, to=180, by=3), extend=TRUE)$surv)
  mort_surv_black[1,]   <- t(summary(mort_km_fit_black, times = seq(from=3, to=180, by=3), extend=TRUE)$surv)
  remove_surv_black[1,] <- t(summary(remove_km_fit_black, times = seq(from=3, to=180, by=3), extend=TRUE)$surv)

  ddtx_surv_black[2,]   <- t(summary(ddtx_km_fit_black, times = seq(from=3, to=180, by=3), extend=TRUE)$std.err)
  ldtx_surv_black[2,]   <- t(summary(ldtx_km_fit_black, times = seq(from=3, to=180, by=3), extend=TRUE)$std.err)
  mort_surv_black[2,]   <- t(summary(mort_km_fit_black, times = seq(from=3, to=180, by=3), extend=TRUE)$std.err)
  remove_surv_black[2,] <- t(summary(remove_km_fit_black, times = seq(from=3, to=180, by=3), extend=TRUE)$std.err)

  # Parametric Fits
  ddtx_weibull_fit_black = flexsurv::flexsurvreg(Surv(time = t0, time2 = t, event_cd==0)~1,
                                       data=df_race[[2]], dist="weibull")
  ldtx_gompertz_fit_black = flexsurv::flexsurvreg(Surv(time = t0, time2 = t, event_cd==1)~1,
                                        data=df_race[[2]], dist="gompertz")
  mort_weibull_fit_black = flexsurv::flexsurvreg(Surv(time = t0, time2 = t, event_cd==2)~1,
                                       data=df_race[[2]], dist="weibull")
  remove_weibull_fit_black = flexsurv::flexsurvreg(Surv(time = t0, time2 = t, event_cd==3)~1,
                                         data=df_race[[2]], dist="weibull")

  ddtx_shape_black <- ddtx_weibull_fit_black$coefficients
  ldtx_shape_black <- ldtx_gompertz_fit_black$coefficients
  mort_shape_black <- mort_weibull_fit_black$coefficients
  remove_shape_black <- remove_weibull_fit_black$coefficients
  ddtx_shape_black_cov <- ddtx_weibull_fit_black$cov
  ldtx_shape_black_cov <- ldtx_gompertz_fit_black$cov
  mort_shape_black_cov <- mort_weibull_fit_black$cov
  remove_shape_black_cov <- remove_weibull_fit_black$cov

  ## Hispanic Candidates
  # Kaplan-Meier
  ddtx_km_fit_hispanic    <- survival::survfit(Surv(time = t0, time2 = t, event_cd==0) ~ 1, data = df_race[[3]])
  ldtx_km_fit_hispanic    <- survival::survfit(Surv(time = t0, time2 = t, event_cd==1) ~ 1, data = df_race[[3]])
  mort_km_fit_hispanic    <- survival::survfit(Surv(time = t0, time2 = t, event_cd==2) ~ 1, data = df_race[[3]])
  remove_km_fit_hispanic  <- survival::survfit(Surv(time = t0, time2 = t, event_cd==3) ~ 1, data = df_race[[3]])

  ddtx_surv_hispanic[1,]   <- t(summary(ddtx_km_fit_hispanic, times = seq(from=3, to=180, by=3), extend=TRUE)$surv)
  ldtx_surv_hispanic[1,]   <- t(summary(ldtx_km_fit_hispanic, times = seq(from=3, to=180, by=3), extend=TRUE)$surv)
  mort_surv_hispanic[1,]   <- t(summary(mort_km_fit_hispanic, times = seq(from=3, to=180, by=3), extend=TRUE)$surv)
  remove_surv_hispanic[1,] <- t(summary(remove_km_fit_hispanic, times = seq(from=3, to=180, by=3), extend=TRUE)$surv)

  ddtx_surv_hispanic[2,]   <- t(summary(ddtx_km_fit_hispanic, times = seq(from=3, to=180, by=3), extend=TRUE)$std.err)
  ldtx_surv_hispanic[2,]   <- t(summary(ldtx_km_fit_hispanic, times = seq(from=3, to=180, by=3), extend=TRUE)$std.err)
  mort_surv_hispanic[2,]   <- t(summary(mort_km_fit_hispanic, times = seq(from=3, to=180, by=3), extend=TRUE)$std.err)
  remove_surv_hispanic[2,] <- t(summary(remove_km_fit_hispanic, times = seq(from=3, to=180, by=3), extend=TRUE)$std.err)

  # Parametric Fits
  ddtx_weibull_fit_hispanic = flexsurv::flexsurvreg(Surv(time = t0, time2 = t, event_cd==0)~1,
                                          data=df_race[[3]], dist="weibull")
  ldtx_gompertz_fit_hispanic = flexsurv::flexsurvreg(Surv(time = t0, time2 = t, event_cd==1)~1,
                                           data=df_race[[3]], dist="gompertz")
  mort_weibull_fit_hispanic = flexsurv::flexsurvreg(Surv(time = t0, time2 = t, event_cd==2)~1,
                                          data=df_race[[3]], dist="weibull")
  remove_weibull_fit_hispanic = flexsurv::flexsurvreg(Surv(time = t0, time2 = t, event_cd==3)~1,
                                            data=df_race[[3]], dist="weibull")

  ddtx_shape_hispanic <- ddtx_weibull_fit_hispanic$coefficients
  ldtx_shape_hispanic <- ldtx_gompertz_fit_hispanic$coefficients
  mort_shape_hispanic <- mort_weibull_fit_hispanic$coefficients
  remove_shape_hispanic <- remove_weibull_fit_hispanic$coefficients
  ddtx_shape_hispanic_cov <- ddtx_weibull_fit_hispanic$cov
  ldtx_shape_hispanic_cov <- ldtx_gompertz_fit_hispanic$cov
  mort_shape_hispanic_cov <- mort_weibull_fit_hispanic$cov
  remove_shape_hispanic_cov <- remove_weibull_fit_hispanic$cov

  ## Other Candidates
  # Kaplan-Meier
  ddtx_km_fit_other    <- survival::survfit(Surv(time = t0, time2 = t, event_cd==0) ~ 1, data = df_race[[4]])
  ldtx_km_fit_other    <- survival::survfit(Surv(time = t0, time2 = t, event_cd==1) ~ 1, data = df_race[[4]])
  mort_km_fit_other    <- survival::survfit(Surv(time = t0, time2 = t, event_cd==2) ~ 1, data = df_race[[4]])
  remove_km_fit_other  <- survival::survfit(Surv(time = t0, time2 = t, event_cd==3) ~ 1, data = df_race[[4]])

  ddtx_surv_other[1,]   <- t(summary(ddtx_km_fit_other, times = seq(from=3, to=180, by=3), extend=TRUE)$surv)
  ldtx_surv_other[1,]   <- t(summary(ldtx_km_fit_other, times = seq(from=3, to=180, by=3), extend=TRUE)$surv)
  mort_surv_other[1,]   <- t(summary(mort_km_fit_other, times = seq(from=3, to=180, by=3), extend=TRUE)$surv)
  remove_surv_other[1,] <- t(summary(remove_km_fit_other, times = seq(from=3, to=180, by=3), extend=TRUE)$surv)

  ddtx_surv_other[2,]   <- t(summary(ddtx_km_fit_other, times = seq(from=3, to=180, by=3), extend=TRUE)$std.err)
  ldtx_surv_other[2,]   <- t(summary(ldtx_km_fit_other, times = seq(from=3, to=180, by=3), extend=TRUE)$std.err)
  mort_surv_other[2,]   <- t(summary(mort_km_fit_other, times = seq(from=3, to=180, by=3), extend=TRUE)$std.err)
  remove_surv_other[2,] <- t(summary(remove_km_fit_other, times = seq(from=3, to=180, by=3), extend=TRUE)$std.err)

  # Parametric Fits
  ddtx_weibull_fit_other = flexsurv::flexsurvreg(Surv(time = t0, time2 = t, event_cd==0)~1,
                                       data=df_race[[4]], dist="weibull")
  ldtx_gompertz_fit_other = flexsurv::flexsurvreg(Surv(time = t0, time2 = t, event_cd==1)~1,
                                        data=df_race[[4]], dist="gompertz")
  mort_weibull_fit_other = flexsurv::flexsurvreg(Surv(time = t0, time2 = t, event_cd==2)~1,
                                       data=df_race[[4]], dist="weibull")
  remove_weibull_fit_other = flexsurv::flexsurvreg(Surv(time = t0, time2 = t, event_cd==3)~1,
                                         data=df_race[[4]], dist="weibull")

  ddtx_shape_other <- ddtx_weibull_fit_other$coefficients
  ldtx_shape_other <- ldtx_gompertz_fit_other$coefficients
  mort_shape_other <- mort_weibull_fit_other$coefficients
  remove_shape_other <- remove_weibull_fit_other$coefficients
  ddtx_shape_other_cov <- ddtx_weibull_fit_other$cov
  ldtx_shape_other_cov <- ldtx_gompertz_fit_other$cov
  mort_shape_other_cov <- mort_weibull_fit_other$cov
  remove_shape_other_cov <- remove_weibull_fit_other$cov

  ### Table 1
  data_tx <- data_tx %>%
    dplyr::mutate(male = ifelse(gender==1, 0, 1)) %>%
    dplyr::mutate(white = ifelse(race_ethnic==1, 1, 0)) %>%
    dplyr::mutate(black = ifelse(race_ethnic==2, 1, 0)) %>%
    dplyr::mutate(hispanic = ifelse(race_ethnic==3, 1, 0)) %>%
    dplyr::mutate(other_race = ifelse(race_ethnic==4, 1, 0)) %>%
    dplyr::mutate(can_blood_a = ifelse(can_abo_cat==1, 1, 0)) %>%
    dplyr::mutate(can_blood_ab = ifelse(can_abo_cat==2, 1, 0)) %>%
    dplyr::mutate(can_blood_b = ifelse(can_abo_cat==3, 1, 0)) %>%
    dplyr::mutate(can_blood_o = ifelse(can_abo_cat==4, 1, 0)) %>%
    dplyr::mutate(optn_reg1 = ifelse(region==1, 1, 0)) %>%
    dplyr::mutate(optn_reg2 = ifelse(region==2, 1, 0)) %>%
    dplyr::mutate(optn_reg3 = ifelse(region==3, 1, 0)) %>%
    dplyr::mutate(optn_reg4 = ifelse(region==4, 1, 0)) %>%
    dplyr::mutate(optn_reg5 = ifelse(region==5, 1, 0)) %>%
    dplyr::mutate(optn_reg6 = ifelse(region==6, 1, 0)) %>%
    dplyr::mutate(optn_reg7 = ifelse(region==7, 1, 0)) %>%
    dplyr::mutate(optn_reg8 = ifelse(region==8, 1, 0)) %>%
    dplyr::mutate(optn_reg9 = ifelse(region==9, 1, 0)) %>%
    dplyr::mutate(optn_reg10 = ifelse(region==10, 1, 0)) %>%
    dplyr::mutate(optn_reg11 = ifelse(region==11, 1, 0))



  sim_table1 <- matrix(data = NA, nrow = 2, ncol = 27)
  sim_table1[1,1] <- mean(data_tx$rec_age_at_tx)
  sim_table1[1,2] <- mean(data_tx$male)
  sim_table1[1,3] <- mean(data_tx$white)
  sim_table1[1,4] <- mean(data_tx$black)
  sim_table1[1,5] <- mean(data_tx$hispanic)
  sim_table1[1,6] <- mean(data_tx$other_race)
  sim_table1[1,7] <- mean(data_tx$can_blood_a)
  sim_table1[1,8] <- mean(data_tx$can_blood_ab)
  sim_table1[1,9] <- mean(data_tx$can_blood_b)
  sim_table1[1,10] <- mean(data_tx$can_blood_o)
  sim_table1[1,11] <- mean(data_tx$years_dial)
  sim_table1[1,12] <- mean(data_tx$diab_stat, na.rm = TRUE)
  sim_table1[1,13] <- mean(data_tx$copd, na.rm = TRUE)
  sim_table1[1,14] <- mean(data_tx$pvd, na.rm = TRUE)
  sim_table1[1,15] <- mean(data_tx$ang_cad, na.rm = TRUE)
  sim_table1[1,16] <- mean(data_tx$canhx_cpra)
  sim_table1[1,17] <- mean(data_tx$optn_reg1, na.rm = TRUE)
  sim_table1[1,18] <- mean(data_tx$optn_reg2, na.rm = TRUE)
  sim_table1[1,19] <- mean(data_tx$optn_reg3, na.rm = TRUE)
  sim_table1[1,20] <- mean(data_tx$optn_reg4, na.rm = TRUE)
  sim_table1[1,21] <- mean(data_tx$optn_reg5, na.rm = TRUE)
  sim_table1[1,22] <- mean(data_tx$optn_reg6, na.rm = TRUE)
  sim_table1[1,23] <- mean(data_tx$optn_reg7, na.rm = TRUE)
  sim_table1[1,24] <- mean(data_tx$optn_reg8, na.rm = TRUE)
  sim_table1[1,25] <- mean(data_tx$optn_reg9, na.rm = TRUE)
  sim_table1[1,26] <- mean(data_tx$optn_reg10, na.rm = TRUE)
  sim_table1[1,27] <- mean(data_tx$optn_reg11, na.rm = TRUE)
  sim_table1[2,1] <- sd(data_tx$rec_age_at_tx)
  sim_table1[2,2] <- sd(data_tx$male)
  sim_table1[2,3] <- sd(data_tx$white)
  sim_table1[2,4] <- sd(data_tx$black)
  sim_table1[2,5] <- sd(data_tx$hispanic)
  sim_table1[2,6] <- sd(data_tx$other_race)
  sim_table1[2,7] <- sd(data_tx$can_blood_a)
  sim_table1[2,8] <- sd(data_tx$can_blood_ab)
  sim_table1[2,9] <- sd(data_tx$can_blood_b)
  sim_table1[2,10] <- sd(data_tx$can_blood_o)
  sim_table1[2,11] <- sd(data_tx$years_dial)
  sim_table1[2,12] <- sd(data_tx$diab_stat, na.rm = TRUE)
  sim_table1[2,13] <- sd(data_tx$copd, na.rm = TRUE)
  sim_table1[2,14] <- sd(data_tx$pvd, na.rm = TRUE)
  sim_table1[2,15] <- sd(data_tx$ang_cad, na.rm = TRUE)
  sim_table1[2,16] <- sd(data_tx$canhx_cpra)
  sim_table1[2,17] <- sd(data_tx$optn_reg1, na.rm = TRUE)
  sim_table1[2,18] <- sd(data_tx$optn_reg2, na.rm = TRUE)
  sim_table1[2,19] <- sd(data_tx$optn_reg3, na.rm = TRUE)
  sim_table1[2,20] <- sd(data_tx$optn_reg4, na.rm = TRUE)
  sim_table1[2,21] <- sd(data_tx$optn_reg5, na.rm = TRUE)
  sim_table1[2,22] <- sd(data_tx$optn_reg6, na.rm = TRUE)
  sim_table1[2,23] <- sd(data_tx$optn_reg7, na.rm = TRUE)
  sim_table1[2,24] <- sd(data_tx$optn_reg8, na.rm = TRUE)
  sim_table1[2,25] <- sd(data_tx$optn_reg9, na.rm = TRUE)
  sim_table1[2,26] <- sd(data_tx$optn_reg10, na.rm = TRUE)
  sim_table1[2,27] <- sd(data_tx$optn_reg11, na.rm = TRUE)

  ### Prepare for data export
  # Combine shape and scale parameters
  m_ddtx_shape <- cbind(t(ddtx_shape_all), t(ddtx_shape_diab0), t(ddtx_shape_diab0),
                        t(ddtx_shape_white), t(ddtx_shape_black),
                        t(ddtx_shape_hispanic), t(ddtx_shape_other))
  m_ldtx_shape <- cbind(t(ldtx_shape_all), t(ldtx_shape_diab0), t(ldtx_shape_diab0),
                        t(ldtx_shape_white), t(ldtx_shape_black),
                        t(ldtx_shape_hispanic), t(ldtx_shape_other))
  m_mort_shape <- cbind(t(mort_shape_all), t(mort_shape_diab0), t(mort_shape_diab0),
                        t(mort_shape_white), t(mort_shape_black),
                        t(mort_shape_hispanic), t(mort_shape_other))
  m_remove_shape <- cbind(t(remove_shape_all), t(remove_shape_diab0), t(remove_shape_diab0),
                          t(remove_shape_white), t(remove_shape_black),
                          t(remove_shape_hispanic), t(remove_shape_other))

  m_ddtx_shape_cov <- cbind(t(ddtx_shape_all_cov), t(ddtx_shape_diab0_cov), t(ddtx_shape_diab0_cov),
                        t(ddtx_shape_white_cov), t(ddtx_shape_black_cov),
                        t(ddtx_shape_hispanic_cov), t(ddtx_shape_other_cov))
  m_ldtx_shape_cov <- cbind(t(ldtx_shape_all_cov), t(ldtx_shape_diab0_cov), t(ldtx_shape_diab0_cov),
                        t(ldtx_shape_white_cov), t(ldtx_shape_black_cov),
                        t(ldtx_shape_hispanic_cov), t(ldtx_shape_other_cov))
  m_mort_shape_cov <- cbind(t(mort_shape_all_cov), t(mort_shape_diab0_cov), t(mort_shape_diab0_cov),
                        t(mort_shape_white_cov), t(mort_shape_black_cov),
                        t(mort_shape_hispanic_cov), t(mort_shape_other_cov))
  m_remove_shape_cov <- cbind(t(remove_shape_all_cov), t(remove_shape_diab0_cov), t(remove_shape_diab0_cov),
                          t(remove_shape_white_cov), t(remove_shape_black_cov),
                          t(remove_shape_hispanic_cov), t(remove_shape_other_cov))

  m_ddtx_shape <- rbind(m_ddtx_shape, m_ddtx_shape_cov)
  m_ldtx_shape <- rbind(m_ldtx_shape, m_ldtx_shape_cov)
  m_mort_shape <- rbind(m_mort_shape, m_mort_shape_cov)
  m_remove_shape <- rbind(m_remove_shape, m_remove_shape_cov)

  return(list(ddtx_surv,
              ldtx_surv,
              mort_surv,
              remove_surv,
              ddtx_surv_diab0,
              ldtx_surv_diab0,
              mort_surv_diab0,
              remove_surv_diab0,
              ddtx_surv_diab1,
              ldtx_surv_diab1,
              mort_surv_diab1,
              remove_surv_diab1,
              ddtx_surv_white,
              ldtx_surv_white,
              mort_surv_white,
              remove_surv_white,
              ddtx_surv_black,
              ldtx_surv_black,
              mort_surv_black,
              remove_surv_black,
              ddtx_surv_hispanic,
              ldtx_surv_hispanic,
              mort_surv_hispanic,
              remove_surv_hispanic,
              ddtx_surv_other,
              ldtx_surv_other,
              mort_surv_other,
              remove_surv_other,
              m_ddtx_shape,
              m_ldtx_shape,
              m_mort_shape,
              m_remove_shape,
              sim_table1))

}


### Load data
df_list_events <- read_dta("/Volumes/Nephrology/Projects/UNOS/KaufmannM/data/list_events.dta")
data_tx <- read_dta("/Volumes/Nephrology/Projects/UNOS/KaufmannM/data/graft_loss.dta")
data_gl_mort <- read_dta("/Volumes/Nephrology/Projects/UNOS/KaufmannM/data/graft_loss_mort.dta")

### Get targets
waitlist_targets <- get_likelihood_targets(data_list = df_list_events,
                                           data_tx = df_tx)
out_path <- "~/Documents/Simulation Model/Model Outputs/"
for(i in 1:length(waitlist_targets)){
  write.csv(waitlist_targets[i], paste0(out_path, "Recalibration/waitlist/target_", i, ".csv"))
}


get_gof_likelihood <- function(targets, sim_results, iter) {

  # All Patients - Survival
  gof_ddtx_all   <- log(dnorm(sim_results[[1]][iter,], mean = as.numeric(targets[[1]][1,]),
                              sd = as.numeric(targets[[1]][2,])))
  gof_ldtx_all   <- log(dnorm(sim_results[[2]][iter,], mean = as.numeric(targets[[2]][1,]),
                              sd = as.numeric(targets[[2]][2,])))
  gof_mort_all   <- log(dnorm(sim_results[[3]][iter,], mean = as.numeric(targets[[3]][1,]),
                              sd = as.numeric(targets[[3]][2,])))
  gof_remove_all <- log(dnorm(sim_results[[4]][iter,], mean = as.numeric(targets[[4]][1,]),
                              sd = as.numeric(targets[[4]][2,])))
  # No Diabetes - Survival
  gof_ddtx_diab0   <- log(dnorm(sim_results[[5]][iter,], mean = as.numeric(targets[[5]][1,]),
                              sd = as.numeric(targets[[5]][2,])))
  gof_ldtx_diab0   <- log(dnorm(sim_results[[6]][iter,], mean = as.numeric(targets[[6]][1,]),
                              sd = as.numeric(targets[[6]][2,])))
  gof_mort_diab0   <- log(dnorm(sim_results[[7]][iter,], mean = as.numeric(targets[[7]][1,]),
                              sd = as.numeric(targets[[7]][2,])))
  gof_remove_diab0 <- log(dnorm(sim_results[[8]][iter,], mean = as.numeric(targets[[8]][1,]),
                              sd = as.numeric(targets[[8]][2,])))
  # Diabetes - Survival
  gof_ddtx_diab1   <- log(dnorm(sim_results[[9]][iter,], mean = as.numeric(targets[[9]][1,]),
                              sd = as.numeric(targets[[9]][2,])))
  gof_ldtx_diab1   <- log(dnorm(sim_results[[10]][iter,], mean = as.numeric(targets[[10]][1,]),
                              sd = as.numeric(targets[[10]][2,])))
  gof_mort_diab1   <- log(dnorm(sim_results[[11]][iter,], mean = as.numeric(targets[[11]][1,]),
                              sd = as.numeric(targets[[11]][2,])))
  gof_remove_diab1 <- log(dnorm(sim_results[[12]][iter,], mean = as.numeric(targets[[12]][1,]),
                              sd = as.numeric(targets[[12]][2,])))
  # White Patients - Survival
  gof_ddtx_white   <- log(dnorm(sim_results[[13]][iter,], mean = as.numeric(targets[[13]][1,]),
                              sd = as.numeric(targets[[13]][2,])))
  gof_ldtx_white   <- log(dnorm(sim_results[[14]][iter,], mean = as.numeric(targets[[14]][1,]),
                              sd = as.numeric(targets[[14]][2,])))
  gof_mort_white   <- log(dnorm(sim_results[[15]][iter,], mean = as.numeric(targets[[15]][1,]),
                              sd = as.numeric(targets[[15]][2,])))
  gof_remove_white <- log(dnorm(sim_results[[16]][iter,], mean = as.numeric(targets[[16]][1,]),
                              sd = as.numeric(targets[[16]][2,])))
  # Black Patients - Survival
  gof_ddtx_black   <- log(dnorm(sim_results[[17]][iter,], mean = as.numeric(targets[[17]][1,]),
                              sd = as.numeric(targets[[17]][2,])))
  gof_ldtx_black   <- log(dnorm(sim_results[[18]][iter,], mean = as.numeric(targets[[18]][1,]),
                              sd = as.numeric(targets[[18]][2,])))
  gof_mort_black   <- log(dnorm(sim_results[[19]][iter,], mean = as.numeric(targets[[19]][1,]),
                              sd = as.numeric(targets[[19]][2,])))
  gof_remove_black <- log(dnorm(sim_results[[20]][iter,], mean = as.numeric(targets[[20]][1,]),
                              sd = as.numeric(targets[[20]][2,])))
  # Hispanic Patients - Survival
  gof_ddtx_hispanic   <- log(dnorm(sim_results[[21]][iter,], mean = as.numeric(targets[[21]][1,]),
                              sd = as.numeric(targets[[21]][2,])))
  gof_ldtx_hispanic   <- log(dnorm(sim_results[[22]][iter,], mean = as.numeric(targets[[22]][1,]),
                              sd = as.numeric(targets[[22]][2,])))
  gof_mort_hispanic   <- log(dnorm(sim_results[[23]][iter,], mean = as.numeric(targets[[23]][1,]),
                              sd = as.numeric(targets[[23]][2,])))
  gof_remove_hispanic <- log(dnorm(sim_results[[24]][iter,], mean = as.numeric(targets[[24]][1,]),
                              sd = as.numeric(targets[[24]][2,])))
  # Other Patients - Survival
  gof_ddtx_other   <- log(dnorm(sim_results[[25]][iter,], mean = as.numeric(targets[[25]][1,]),
                              sd = as.numeric(targets[[25]][2,])))
  gof_ldtx_other   <- log(dnorm(sim_results[[26]][iter,], mean = as.numeric(targets[[26]][1,]),
                              sd = as.numeric(targets[[26]][2,])))
  gof_mort_other   <- log(dnorm(sim_results[[27]][iter,], mean = as.numeric(targets[[27]][1,]),
                              sd = as.numeric(targets[[27]][2,])))
  gof_remove_other <- log(dnorm(sim_results[[28]][iter,], mean = as.numeric(targets[[28]][1,]),
                              sd = as.numeric(targets[[28]][2,])))

  ### Shape parameters
  # All Patients
  gof_ddtx_shape_all   <- log(dmvnorm(sim_results[[29]][iter,1:2],
                                      mean = as.numeric(targets[[29]][1,1:2]),
                                      sigma = as.matrix(targets[[29]][2:3,1:2])))
  gof_ldtx_shape_all   <- log(dmvnorm(sim_results[[30]][iter,1:2],
                                      mean = as.numeric(targets[[30]][1,1:2]),
                                      sigma = as.matrix(targets[[30]][2:3,1:2])))
  gof_mort_shape_all   <- log(dmvnorm(sim_results[[31]][iter,1:2],
                                      mean = as.numeric(targets[[31]][1,1:2]),
                                      sigma = as.matrix(targets[[31]][2:3,1:2])))
  gof_remove_shape_all <- log(dmvnorm(sim_results[[32]][iter,1:2],
                                      mean = as.numeric(targets[[32]][1,1:2]),
                                      sigma = as.matrix(targets[[32]][2:3,1:2])))
  # No Diabetes
  gof_ddtx_shape_diab0   <- log(dmvnorm(sim_results[[29]][iter,3:4],
                                        mean = as.numeric(targets[[29]][1,3:4]),
                                        sigma = as.matrix(targets[[29]][2:3,3:4])))
  gof_ldtx_shape_diab0   <- log(dmvnorm(sim_results[[30]][iter,3:4],
                                        mean = as.numeric(targets[[30]][1,3:4]),
                                        sigma = as.matrix(targets[[30]][2:3,3:4])))
  gof_mort_shape_diab0   <- log(dmvnorm(sim_results[[31]][iter,3:4],
                                        mean = as.numeric(targets[[31]][1,3:4]),
                                        sigma = as.matrix(targets[[31]][2:3,3:4])))
  gof_remove_shape_diab0 <- log(dmvnorm(sim_results[[32]][iter,3:4],
                                        mean = as.numeric(targets[[32]][1,3:4]),
                                        sigma = as.matrix(targets[[32]][2:3,3:4])))
  # Diabetes
  gof_ddtx_shape_diab1   <- log(dmvnorm(sim_results[[29]][iter,5:6],
                                        mean = as.numeric(targets[[29]][1,5:6]),
                                        sigma = as.matrix(targets[[29]][2:3,5:6])))
  gof_ldtx_shape_diab1   <- log(dmvnorm(sim_results[[30]][iter,5:6],
                                        mean = as.numeric(targets[[30]][1,5:6]),
                                        sigma = as.matrix(targets[[30]][2:3,5:6])))
  gof_mort_shape_diab1   <- log(dmvnorm(sim_results[[31]][iter,5:6],
                                        mean = as.numeric(targets[[31]][1,5:6]),
                                        sigma = as.matrix(targets[[31]][2:3,5:6])))
  gof_remove_shape_diab1 <- log(dmvnorm(sim_results[[32]][iter,5:6],
                                        mean = as.numeric(targets[[32]][1,5:6]),
                                        sigma = as.matrix(targets[[32]][2:3,5:6])))
  # White Patients
  gof_ddtx_shape_white   <- log(dmvnorm(sim_results[[29]][iter,7:8],
                                        mean = as.numeric(targets[[29]][1,7:8]),
                                        sigma = as.matrix(targets[[29]][2:3,7:8])))
  gof_ldtx_shape_white   <- log(dmvnorm(sim_results[[30]][iter,7:8],
                                        mean = as.numeric(targets[[30]][1,7:8]),
                                        sigma = as.matrix(targets[[30]][2:3,7:8])))
  gof_mort_shape_white   <- log(dmvnorm(sim_results[[31]][iter,7:8],
                                        mean = as.numeric(targets[[31]][1,7:8]),
                                        sigma = as.matrix(targets[[31]][2:3,7:8])))
  gof_remove_shape_white <- log(dmvnorm(sim_results[[32]][iter,7:8],
                                        mean = as.numeric(targets[[32]][1,7:8]),
                                        sigma = as.matrix(targets[[32]][2:3,7:8])))
  # Black Patients
  gof_ddtx_shape_black   <- log(dmvnorm(sim_results[[29]][iter,9:10],
                                        mean = as.numeric(targets[[29]][1,9:10]),
                                        sigma = as.matrix(targets[[29]][2:3,9:10])))
  gof_ldtx_shape_black   <- log(dmvnorm(sim_results[[30]][iter,9:10],
                                        mean = as.numeric(targets[[30]][1,9:10]),
                                        sigma = as.matrix(targets[[30]][2:3,9:10])))
  gof_mort_shape_black   <- log(dmvnorm(sim_results[[31]][iter,9:10],
                                        mean = as.numeric(targets[[31]][1,9:10]),
                                        sigma = as.matrix(targets[[31]][2:3,9:10])))
  gof_remove_shape_black <- log(dmvnorm(sim_results[[32]][iter,9:10],
                                        mean = as.numeric(targets[[32]][1,9:10]),
                                        sigma = as.matrix(targets[[32]][2:3,9:10])))
  # Hispanic Patients
  gof_ddtx_shape_hispanic   <- log(dmvnorm(sim_results[[29]][iter,11:12],
                                           mean = as.numeric(targets[[29]][1,11:12]),
                                           sigma = as.matrix(targets[[29]][2:3,11:12])))
  gof_ldtx_shape_hispanic   <- log(dmvnorm(sim_results[[30]][iter,11:12],
                                           mean = as.numeric(targets[[30]][1,11:12]),
                                           sigma = as.matrix(targets[[30]][2:3,11:12])))
  gof_mort_shape_hispanic   <- log(dmvnorm(sim_results[[31]][iter,11:12],
                                           mean = as.numeric(targets[[31]][1,11:12]),
                                           sigma = as.matrix(targets[[31]][2:3,11:12])))
  gof_remove_shape_hispanic <- log(dmvnorm(sim_results[[32]][iter,11:12],
                                           mean = as.numeric(targets[[32]][1,11:12]),
                                           sigma = as.matrix(targets[[32]][2:3,11:12])))
  # Other Patients
  gof_ddtx_shape_other   <- log(dmvnorm(sim_results[[29]][iter,13:14],
                                        mean = as.numeric(targets[[29]][1,13:14]),
                                        sigma = as.matrix(targets[[29]][2:3,13:14])))
  gof_ldtx_shape_other   <- log(dmvnorm(sim_results[[30]][iter,13:14],
                                        mean = as.numeric(targets[[30]][1,13:14]),
                                        sigma = as.matrix(targets[[30]][2:3,13:14])))
  gof_mort_shape_other   <- log(dmvnorm(sim_results[[31]][iter,13:14],
                                        mean = as.numeric(targets[[31]][1,13:14]),
                                        sigma = as.matrix(targets[[31]][2:3,13:14])))
  gof_remove_shape_other <- log(dmvnorm(sim_results[[32]][iter,13:14],
                                        mean = as.numeric(targets[[32]][1,13:14]),
                                        sigma = as.matrix(targets[[32]][2:3,13:14])))

  ### Table 1
  gof_table1 <- log(dnorm(sim_results[[33]][iter,], mean = as.numeric(targets[[33]][1,]),
                          sd = as.numeric(targets[[33]][2,])))

  ### Sum GOF
  gof_surv  <- -1 * sum(gof_ddtx_all, gof_ldtx_all, gof_mort_all, gof_remove_all,
                        gof_ddtx_diab0, gof_ldtx_diab0, gof_mort_diab0, gof_remove_diab0,
                        gof_ddtx_diab1, gof_ldtx_diab1, gof_mort_diab1, gof_remove_diab1,
                        gof_ddtx_white, gof_ldtx_white, gof_mort_white, gof_remove_white,
                        gof_ddtx_black, gof_ldtx_black, gof_mort_black, gof_remove_black,
                        gof_ddtx_hispanic, gof_ldtx_hispanic, gof_mort_hispanic, gof_remove_hispanic,
                        gof_ddtx_other, gof_ldtx_other, gof_mort_other, gof_remove_other, na.rm = TRUE)

  gof_shape <- -1 * sum(gof_ddtx_shape_all, gof_ldtx_shape_all, gof_mort_shape_all, gof_remove_shape_all,
                        gof_ddtx_shape_diab0, gof_ldtx_shape_diab0, gof_mort_shape_diab0, gof_remove_shape_diab0,
                        gof_ddtx_shape_diab1, gof_ldtx_shape_diab1, gof_mort_shape_diab1, gof_remove_shape_diab1,
                        gof_ddtx_shape_white, gof_ldtx_shape_white, gof_mort_shape_white, gof_remove_shape_white,
                        gof_ddtx_shape_black, gof_ldtx_shape_black, gof_mort_shape_black, gof_remove_shape_black,
                        gof_ddtx_shape_hispanic, gof_ldtx_shape_hispanic, gof_mort_shape_hispanic, gof_remove_shape_hispanic,
                        gof_ddtx_shape_other, gof_ldtx_shape_other, gof_mort_shape_other, gof_remove_shape_other)

  gof_table1 <- -1 * sum(gof_table1)

  return(list(gof_surv, gof_shape, gof_table1))
}

get_gof_all_likelihood <- function(targets, sim_results, iter) {

  # All Patients - Survival
  gof_ddtx_all   <- log(dnorm(as.numeric(sim_results[[1]][iter,]), 
                              mean = as.numeric(targets[[1]][1,]),
                              sd = as.numeric(targets[[1]][2,])))
  gof_ldtx_all   <- log(dnorm(as.numeric(sim_results[[2]][iter,]), 
                              mean = as.numeric(targets[[2]][1,]),
                              sd = as.numeric(targets[[2]][2,])))
  gof_mort_all   <- log(dnorm(as.numeric(sim_results[[3]][iter,]), 
                              mean = as.numeric(targets[[3]][1,]),
                              sd = as.numeric(targets[[3]][2,])))
  gof_remove_all <- log(dnorm(as.numeric(sim_results[[4]][iter,]), 
                              mean = as.numeric(targets[[4]][1,]),
                              sd = as.numeric(targets[[4]][2,])))
  # No Diabetes - Survival
  gof_ddtx_diab0   <- log(dnorm(as.numeric(sim_results[[5]][iter,]), 
                                mean = as.numeric(targets[[5]][1,]),
                                sd = as.numeric(targets[[5]][2,])))
  gof_ldtx_diab0   <- log(dnorm(as.numeric(sim_results[[6]][iter,]), 
                                mean = as.numeric(targets[[6]][1,]),
                                sd = as.numeric(targets[[6]][2,])))
  gof_mort_diab0   <- log(dnorm(as.numeric(sim_results[[7]][iter,]), 
                                mean = as.numeric(targets[[7]][1,]),
                                sd = as.numeric(targets[[7]][2,])))
  gof_remove_diab0 <- log(dnorm(as.numeric(sim_results[[8]][iter,]), 
                                mean = as.numeric(targets[[8]][1,]),
                                sd = as.numeric(targets[[8]][2,])))
  # Diabetes - Survival
  gof_ddtx_diab1   <- log(dnorm(as.numeric(sim_results[[9]][iter,]), 
                                mean = as.numeric(targets[[9]][1,]),
                                sd = as.numeric(targets[[9]][2,])))
  gof_ldtx_diab1   <- log(dnorm(as.numeric(sim_results[[10]][iter,]), 
                                mean = as.numeric(targets[[10]][1,]),
                                sd = as.numeric(targets[[10]][2,])))
  gof_mort_diab1   <- log(dnorm(as.numeric(sim_results[[11]][iter,]), 
                                mean = as.numeric(targets[[11]][1,]),
                                sd = as.numeric(targets[[11]][2,])))
  gof_remove_diab1 <- log(dnorm(as.numeric(sim_results[[12]][iter,]), 
                                mean = as.numeric(targets[[12]][1,]),
                                sd = as.numeric(targets[[12]][2,])))
  # White Patients - Survival
  gof_ddtx_white   <- log(dnorm(as.numeric(sim_results[[13]][iter,]), 
                                mean = as.numeric(targets[[13]][1,]),
                                sd = as.numeric(targets[[13]][2,])))
  gof_ldtx_white   <- log(dnorm(as.numeric(sim_results[[14]][iter,]), 
                                mean = as.numeric(targets[[14]][1,]),
                                sd = as.numeric(targets[[14]][2,])))
  gof_mort_white   <- log(dnorm(as.numeric(sim_results[[15]][iter,]), 
                                mean = as.numeric(targets[[15]][1,]),
                                sd = as.numeric(targets[[15]][2,])))
  gof_remove_white <- log(dnorm(as.numeric(sim_results[[16]][iter,]), 
                                mean = as.numeric(targets[[16]][1,]),
                                sd = as.numeric(targets[[16]][2,])))
  # Black Patients - Survival
  gof_ddtx_black   <- log(dnorm(as.numeric(sim_results[[17]][iter,]), 
                                mean = as.numeric(targets[[17]][1,]),
                                sd = as.numeric(targets[[17]][2,])))
  gof_ldtx_black   <- log(dnorm(as.numeric(sim_results[[18]][iter,]), 
                                mean = as.numeric(targets[[18]][1,]),
                                sd = as.numeric(targets[[18]][2,])))
  gof_mort_black   <- log(dnorm(as.numeric(sim_results[[19]][iter,]), 
                                mean = as.numeric(targets[[19]][1,]),
                                sd = as.numeric(targets[[19]][2,])))
  gof_remove_black <- log(dnorm(as.numeric(sim_results[[20]][iter,]), 
                                mean = as.numeric(targets[[20]][1,]),
                                sd = as.numeric(targets[[20]][2,])))
  # Hispanic Patients - Survival
  gof_ddtx_hispanic   <- log(dnorm(as.numeric(sim_results[[21]][iter,]), 
                                   mean = as.numeric(targets[[21]][1,]),
                                   sd = as.numeric(targets[[21]][2,])))
  gof_ldtx_hispanic   <- log(dnorm(as.numeric(sim_results[[22]][iter,]), 
                                   mean = as.numeric(targets[[22]][1,]),
                                   sd = as.numeric(targets[[22]][2,])))
  gof_mort_hispanic   <- log(dnorm(as.numeric(sim_results[[23]][iter,]), 
                                   mean = as.numeric(targets[[23]][1,]),
                                   sd = as.numeric(targets[[23]][2,])))
  gof_remove_hispanic <- log(dnorm(as.numeric(sim_results[[24]][iter,]), 
                                   mean = as.numeric(targets[[24]][1,]),
                                   sd = as.numeric(targets[[24]][2,])))
  # Other Patients - Survival
  gof_ddtx_other   <- log(dnorm(as.numeric(sim_results[[25]][iter,]), 
                                mean = as.numeric(targets[[25]][1,]),
                                sd = as.numeric(targets[[25]][2,])))
  gof_ldtx_other   <- log(dnorm(as.numeric(sim_results[[26]][iter,]), 
                                mean = as.numeric(targets[[26]][1,]),
                                sd = as.numeric(targets[[26]][2,])))
  gof_mort_other   <- log(dnorm(as.numeric(sim_results[[27]][iter,]), 
                                mean = as.numeric(targets[[27]][1,]),
                                sd = as.numeric(targets[[27]][2,])))
  gof_remove_other <- log(dnorm(as.numeric(sim_results[[28]][iter,]), 
                                mean = as.numeric(targets[[28]][1,]),
                                sd = as.numeric(targets[[28]][2,])))

  ### Shape parameters
  # All Patients
  gof_ddtx_shape_all   <- log(dmvnorm(as.numeric(sim_results[[29]][iter,1:2]),
                                      mean = as.numeric(targets[[29]][1,1:2]),
                                      sigma = as.matrix(targets[[29]][2:3,1:2])))
  gof_ldtx_shape_all   <- log(dmvnorm(as.numeric(sim_results[[30]][iter,1:2]),
                                      mean = as.numeric(targets[[30]][1,1:2]),
                                      sigma = as.matrix(targets[[30]][2:3,1:2])))
  gof_mort_shape_all   <- log(dmvnorm(as.numeric(sim_results[[31]][iter,1:2]),
                                      mean = as.numeric(targets[[31]][1,1:2]),
                                      sigma = as.matrix(targets[[31]][2:3,1:2])))
  gof_remove_shape_all <- log(dmvnorm(as.numeric(sim_results[[32]][iter,1:2]),
                                      mean = as.numeric(targets[[32]][1,1:2]),
                                      sigma = as.matrix(targets[[32]][2:3,1:2])))
  # No Diabetes
  gof_ddtx_shape_diab0   <- log(dmvnorm(as.numeric(sim_results[[29]][iter,3:4]),
                                        mean = as.numeric(targets[[29]][1,3:4]),
                                        sigma = as.matrix(targets[[29]][2:3,3:4])))
  gof_ldtx_shape_diab0   <- log(dmvnorm(as.numeric(sim_results[[30]][iter,3:4]),
                                        mean = as.numeric(targets[[30]][1,3:4]),
                                        sigma = as.matrix(targets[[30]][2:3,3:4])))
  gof_mort_shape_diab0   <- log(dmvnorm(as.numeric(sim_results[[31]][iter,3:4]),
                                        mean = as.numeric(targets[[31]][1,3:4]),
                                        sigma = as.matrix(targets[[31]][2:3,3:4])))
  gof_remove_shape_diab0 <- log(dmvnorm(as.numeric(sim_results[[32]][iter,3:4]),
                                        mean = as.numeric(targets[[32]][1,3:4]),
                                        sigma = as.matrix(targets[[32]][2:3,3:4])))
  # Diabetes
  gof_ddtx_shape_diab1   <- log(dmvnorm(as.numeric(sim_results[[29]][iter,5:6]),
                                        mean = as.numeric(targets[[29]][1,5:6]),
                                        sigma = as.matrix(targets[[29]][2:3,5:6])))
  gof_ldtx_shape_diab1   <- log(dmvnorm(as.numeric(sim_results[[30]][iter,5:6]),
                                        mean = as.numeric(targets[[30]][1,5:6]),
                                        sigma = as.matrix(targets[[30]][2:3,5:6])))
  gof_mort_shape_diab1   <- log(dmvnorm(as.numeric(sim_results[[31]][iter,5:6]),
                                        mean = as.numeric(targets[[31]][1,5:6]),
                                        sigma = as.matrix(targets[[31]][2:3,5:6])))
  gof_remove_shape_diab1 <- log(dmvnorm(as.numeric(sim_results[[32]][iter,5:6]),
                                        mean = as.numeric(targets[[32]][1,5:6]),
                                        sigma = as.matrix(targets[[32]][2:3,5:6])))
  # White Patients
  gof_ddtx_shape_white   <- log(dmvnorm(as.numeric(sim_results[[29]][iter,7:8]),
                                        mean = as.numeric(targets[[29]][1,7:8]),
                                        sigma = as.matrix(targets[[29]][2:3,7:8])))
  gof_ldtx_shape_white   <- log(dmvnorm(as.numeric(sim_results[[30]][iter,7:8]),
                                        mean = as.numeric(targets[[30]][1,7:8]),
                                        sigma = as.matrix(targets[[30]][2:3,7:8])))
  gof_mort_shape_white   <- log(dmvnorm(as.numeric(sim_results[[31]][iter,7:8]),
                                        mean = as.numeric(targets[[31]][1,7:8]),
                                        sigma = as.matrix(targets[[31]][2:3,7:8])))
  gof_remove_shape_white <- log(dmvnorm(as.numeric(sim_results[[32]][iter,7:8]),
                                        mean = as.numeric(targets[[32]][1,7:8]),
                                        sigma = as.matrix(targets[[32]][2:3,7:8])))
  # Black Patients
  gof_ddtx_shape_black   <- log(dmvnorm(as.numeric(sim_results[[29]][iter,9:10]),
                                        mean = as.numeric(targets[[29]][1,9:10]),
                                        sigma = as.matrix(targets[[29]][2:3,9:10])))
  gof_ldtx_shape_black   <- log(dmvnorm(as.numeric(sim_results[[30]][iter,9:10]),
                                        mean = as.numeric(targets[[30]][1,9:10]),
                                        sigma = as.matrix(targets[[30]][2:3,9:10])))
  gof_mort_shape_black   <- log(dmvnorm(as.numeric(sim_results[[31]][iter,9:10]),
                                        mean = as.numeric(targets[[31]][1,9:10]),
                                        sigma = as.matrix(targets[[31]][2:3,9:10])))
  gof_remove_shape_black <- log(dmvnorm(as.numeric(sim_results[[32]][iter,9:10]),
                                        mean = as.numeric(targets[[32]][1,9:10]),
                                        sigma = as.matrix(targets[[32]][2:3,9:10])))
  # Hispanic Patients
  gof_ddtx_shape_hispanic   <- log(dmvnorm(as.numeric(sim_results[[29]][iter,11:12]),
                                           mean = as.numeric(targets[[29]][1,11:12]),
                                           sigma = as.matrix(targets[[29]][2:3,11:12])))
  gof_ldtx_shape_hispanic   <- log(dmvnorm(as.numeric(sim_results[[30]][iter,11:12]),
                                           mean = as.numeric(targets[[30]][1,11:12]),
                                           sigma = as.matrix(targets[[30]][2:3,11:12])))
  gof_mort_shape_hispanic   <- log(dmvnorm(as.numeric(sim_results[[31]][iter,11:12]),
                                           mean = as.numeric(targets[[31]][1,11:12]),
                                           sigma = as.matrix(targets[[31]][2:3,11:12])))
  gof_remove_shape_hispanic <- log(dmvnorm(as.numeric(sim_results[[32]][iter,11:12]),
                                           mean = as.numeric(targets[[32]][1,11:12]),
                                           sigma = as.matrix(targets[[32]][2:3,11:12])))
  # Other Patients
  gof_ddtx_shape_other   <- log(dmvnorm(as.numeric(sim_results[[29]][iter,13:14]),
                                        mean = as.numeric(targets[[29]][1,13:14]),
                                        sigma = as.matrix(targets[[29]][2:3,13:14])))
  gof_ldtx_shape_other   <- log(dmvnorm(as.numeric(sim_results[[30]][iter,13:14]),
                                        mean = as.numeric(targets[[30]][1,13:14]),
                                        sigma = as.matrix(targets[[30]][2:3,13:14])))
  gof_mort_shape_other   <- log(dmvnorm(as.numeric(sim_results[[31]][iter,13:14]),
                                        mean = as.numeric(targets[[31]][1,13:14]),
                                        sigma = as.matrix(targets[[31]][2:3,13:14])))
  gof_remove_shape_other <- log(dmvnorm(as.numeric(sim_results[[32]][iter,13:14]),
                                        mean = as.numeric(targets[[32]][1,13:14]),
                                        sigma = as.matrix(targets[[32]][2:3,13:14])))

  ### Table 1
  gof_table1 <- log(dnorm(as.numeric(sim_results[[33]][iter,]), 
                          mean = as.numeric(targets[[33]][1,]),
                          sd = as.numeric(targets[[33]][2,])))

  ### 30-Day Outcomes
  # gof_30_all      <- LaplacesDemon::ddirichlet(sim_results[[34]][iter,1:4], as.numeric(targets[[34]][1,]), log = TRUE)
  # gof_30_diab0    <- LaplacesDemon::ddirichlet(sim_results[[34]][iter,5:8], as.numeric(targets[[34]][2,]), log = TRUE)
  # gof_30_diab1    <- LaplacesDemon::ddirichlet(sim_results[[34]][iter,9:12], as.numeric(targets[[34]][3,]), log = TRUE)
  # gof_30_white    <- LaplacesDemon::ddirichlet(sim_results[[34]][iter,13:16], as.numeric(targets[[34]][4,]), log = TRUE)
  # gof_30_black    <- LaplacesDemon::ddirichlet(sim_results[[34]][iter,17:20], as.numeric(targets[[34]][5,]), log = TRUE)
  # gof_30_hispanic <- LaplacesDemon::ddirichlet(sim_results[[34]][iter,21:24], as.numeric(targets[[34]][6,]), log = TRUE)
  # gof_30_other    <- LaplacesDemon::ddirichlet(sim_results[[34]][iter,25:28], as.numeric(targets[[34]][7,]), log = TRUE)

  scaled_30_all      <- round(as.numeric(sim_results[[34]][iter,1:4]) * sum(targets[[34]][1,]))
  scaled_30_diab0    <- round(as.numeric(sim_results[[34]][iter,5:8]) * sum(targets[[34]][2,]))
  scaled_30_diab1    <- round(as.numeric(sim_results[[34]][iter,9:12]) * sum(targets[[34]][3,]))
  scaled_30_white    <- round(as.numeric(sim_results[[34]][iter,13:16]) * sum(targets[[34]][4,]))
  scaled_30_black    <- round(as.numeric(sim_results[[34]][iter,17:20]) * sum(targets[[34]][5,]))
  scaled_30_hispanic <- round(as.numeric(sim_results[[34]][iter,21:24]) * sum(targets[[34]][6,]))
  scaled_30_other    <- round(as.numeric(sim_results[[34]][iter,25:28]) * sum(targets[[34]][7,]))

  diff_30_all      <- sum(scaled_30_all) - sum(targets[[34]][1,])
  diff_30_diab0    <- sum(scaled_30_diab0) - sum(targets[[34]][2,])
  diff_30_diab1    <- sum(scaled_30_diab1) - sum(targets[[34]][3,])
  diff_30_white    <- sum(scaled_30_white) - sum(targets[[34]][4,])
  diff_30_black    <- sum(scaled_30_black) - sum(targets[[34]][5,])
  diff_30_hispanic <- sum(scaled_30_hispanic) - sum(targets[[34]][6,])
  diff_30_other    <- sum(scaled_30_other) - sum(targets[[34]][7,])

  n_30_all_zero      <- sum(scaled_30_all==0)
  n_30_diab0_zero    <- sum(scaled_30_diab0==0)
  n_30_diab1_zero    <- sum(scaled_30_diab1==0)
  n_30_white_zero    <- sum(scaled_30_white==0)
  n_30_black_zero    <- sum(scaled_30_black==0)
  n_30_hispanic_zero <- sum(scaled_30_hispanic==0)
  n_30_other_zero    <- sum(scaled_30_other==0)

  scaled_30_all[scaled_30_all==0] <- 1
  scaled_30_diab0[scaled_30_diab0==0] <- 1
  scaled_30_diab1[scaled_30_diab1==0] <- 1
  scaled_30_white[scaled_30_white==0] <- 1
  scaled_30_black[scaled_30_black==0] <- 1
  scaled_30_hispanic[scaled_30_hispanic==0] <- 1
  scaled_30_other[scaled_30_other==0] <- 1

  prop_scaled_30_all      <- scaled_30_all / sum(scaled_30_all)
  prop_scaled_30_diab0    <- scaled_30_diab0 / sum(scaled_30_diab0)
  prop_scaled_30_diab1    <- scaled_30_diab1 / sum(scaled_30_diab1)
  prop_scaled_30_white    <- scaled_30_white / sum(scaled_30_white)
  prop_scaled_30_black    <- scaled_30_black / sum(scaled_30_black)
  prop_scaled_30_hispanic <- scaled_30_hispanic / sum(scaled_30_hispanic)
  prop_scaled_30_other    <- scaled_30_other / sum(scaled_30_other)

  gof_30_all      <- LaplacesDemon::ddirichlet(prop_scaled_30_all, as.numeric(targets[[34]][1,]), log = TRUE)
  gof_30_diab0    <- LaplacesDemon::ddirichlet(prop_scaled_30_diab0, as.numeric(targets[[34]][2,]), log = TRUE)
  gof_30_diab1    <- LaplacesDemon::ddirichlet(prop_scaled_30_diab1, as.numeric(targets[[34]][3,]), log = TRUE)
  gof_30_white    <- LaplacesDemon::ddirichlet(prop_scaled_30_white, as.numeric(targets[[34]][4,]), log = TRUE)
  gof_30_black    <- LaplacesDemon::ddirichlet(prop_scaled_30_black, as.numeric(targets[[34]][5,]), log = TRUE)
  gof_30_hispanic <- LaplacesDemon::ddirichlet(prop_scaled_30_hispanic, as.numeric(targets[[34]][6,]), log = TRUE)
  gof_30_other    <- LaplacesDemon::ddirichlet(prop_scaled_30_other, as.numeric(targets[[34]][7,]), log = TRUE)

  gof_30 <- -1 * sum(gof_30_all, gof_30_diab0, gof_30_diab1, gof_30_white,
                     gof_30_black, gof_30_hispanic, gof_30_other)
  ### Death with Function
  ## Survival-based GOF
  gof_gs_mort_all <- dnorm(as.numeric(sim_results[[35]][iter,]),
                               mean = as.numeric(targets[[35]][1,]),
                               sd = as.numeric(targets[[35]][2,]), log = T)
  gof_gs_mort_diab0 <- dnorm(as.numeric(sim_results[[36]][iter,]),
                                 mean = as.numeric(targets[[36]][1,]),
                                 sd = as.numeric(targets[[36]][2,]), log = T)
  gof_gs_mort_diab1 <- dnorm(as.numeric(sim_results[[37]][iter,]),
                                 mean = as.numeric(targets[[37]][1,]),
                                 sd = as.numeric(targets[[37]][2,]), log = T)
  gof_gs_mort_white <- dnorm(as.numeric(sim_results[[38]][iter,]),
                                 mean = as.numeric(targets[[38]][1,]),
                                 sd = as.numeric(targets[[38]][2,]), log = T)
  gof_gs_mort_black <- dnorm(as.numeric(sim_results[[39]][iter,]),
                                 mean = as.numeric(targets[[39]][1,]),
                                 sd = as.numeric(targets[[39]][2,]), log = T)
  gof_gs_mort_hispanic <- dnorm(as.numeric(sim_results[[40]][iter,]),
                                    mean = as.numeric(targets[[40]][1,]),
                                    sd = as.numeric(targets[[40]][2,]), log = T)
  gof_gs_mort_other <- dnorm(as.numeric(sim_results[[41]][iter,]),
                                 mean = as.numeric(targets[[41]][1,]),
                                 sd = as.numeric(targets[[41]][2,]), log = T)

  ## Shape-based GOF
  gof_gs_mort_shape_all <- dmvnorm(as.numeric(sim_results[[56]][iter,1:2]),
                                       mean = as.numeric(targets[[42]][1,1:2]),
                                       sigma = as.matrix(targets[[42]][2:3,1:2]), log = T)
  gof_gs_mort_shape_diab0 <- dmvnorm(as.numeric(sim_results[[56]][iter,3:4]),
                                       mean = as.numeric(targets[[42]][1,3:4]),
                                       sigma = as.matrix(targets[[42]][2:3,3:4]), log = T)
  gof_gs_mort_shape_diab1 <- dmvnorm(as.numeric(sim_results[[56]][iter,5:6]),
                                       mean = as.numeric(targets[[42]][1,5:6]),
                                       sigma = as.matrix(targets[[42]][2:3,5:6]), log = T)
  gof_gs_mort_shape_white <- dmvnorm(as.numeric(sim_results[[56]][iter,7:8]),
                                         mean = as.numeric(targets[[42]][1,7:8]),
                                         sigma = as.matrix(targets[[42]][2:3,7:8]), log = T)
  gof_gs_mort_shape_black <- dmvnorm(as.numeric(sim_results[[56]][iter,9:10]),
                                         mean = as.numeric(targets[[42]][1,9:10]),
                                         sigma = as.matrix(targets[[42]][2:3,9:10]), log = T)
  gof_gs_mort_shape_hispanic <- dmvnorm(as.numeric(sim_results[[56]][iter,11:12]),
                                         mean = as.numeric(targets[[42]][1,11:12]),
                                         sigma = as.matrix(targets[[42]][2:3,11:12]), log = T)
  gof_gs_mort_shape_other <- dmvnorm(as.numeric(sim_results[[56]][iter,13:14]),
                                         mean = as.numeric(targets[[42]][1,13:14]),
                                         sigma = as.matrix(targets[[42]][2:3,13:14]), log = T)

  ### Graft Loss
  ## Survival-based GOF
  gof_gl_all <- dnorm(as.numeric(sim_results[[42]][iter,]),
                          mean = as.numeric(targets[[43]][1,]),
                          sd = as.numeric(targets[[43]][2,]), log = T)
  gof_gl_diab0 <- dnorm(as.numeric(sim_results[[43]][iter,]),
                            mean = as.numeric(targets[[44]][1,]),
                            sd = as.numeric(targets[[44]][2,]), log = T)
  gof_gl_diab1 <- dnorm(as.numeric(sim_results[[44]][iter,]),
                            mean = as.numeric(targets[[45]][1,]),
                            sd = as.numeric(targets[[45]][2,]), log = T)
  gof_gl_white <- dnorm(as.numeric(sim_results[[45]][iter,]),
                            mean = as.numeric(targets[[46]][1,]),
                            sd = as.numeric(targets[[46]][2,]), log = T)
  gof_gl_black <- dnorm(as.numeric(sim_results[[46]][iter,]),
                            mean = as.numeric(targets[[47]][1,]),
                            sd = as.numeric(targets[[47]][2,]), log = T)
  gof_gl_hispanic <- dnorm(as.numeric(sim_results[[47]][iter,]),
                               mean = as.numeric(targets[[48]][1,]),
                               sd = as.numeric(targets[[48]][2,]), log = T)
  gof_gl_other <- dnorm(as.numeric(sim_results[[48]][iter,]),
                            mean = as.numeric(targets[[49]][1,]),
                            sd = as.numeric(targets[[49]][2,]), log = T)

  ## Shape-based GOF
  gof_gl_shape_all <- dmvnorm(as.numeric(sim_results[[57]][iter,1:2]),
                              mean = as.numeric(targets[[50]][1,1:2]),
                              sigma = as.matrix(targets[[50]][2:3,1:2]), log = T)
  gof_gl_shape_diab0 <- dmvnorm(as.numeric(sim_results[[57]][iter,3:4]),
                                         mean = as.numeric(targets[[50]][1,3:4]),
                                         sigma = as.matrix(targets[[50]][2:3,3:4]), log = T)
  gof_gl_shape_diab1 <- dmvnorm(as.numeric(sim_results[[57]][iter,5:6]),
                                         mean = as.numeric(targets[[50]][1,5:6]),
                                         sigma = as.matrix(targets[[50]][2:3,5:6]), log = T)
  gof_gl_shape_white <- dmvnorm(as.numeric(sim_results[[57]][iter,7:8]),
                                         mean = as.numeric(targets[[50]][1,7:8]),
                                         sigma = as.matrix(targets[[50]][2:3,7:8]), log = T)
  gof_gl_shape_black <- dmvnorm(as.numeric(sim_results[[57]][iter,9:10]),
                                         mean = as.numeric(targets[[50]][1,9:10]),
                                         sigma = as.matrix(targets[[50]][2:3,9:10]), log = T)
  gof_gl_shape_hispanic <- dmvnorm(as.numeric(sim_results[[57]][iter,11:12]),
                                            mean = as.numeric(targets[[50]][1,11:12]),
                                            sigma = as.matrix(targets[[50]][2:3,11:12]), log = T)
  gof_gl_shape_other <- dmvnorm(as.numeric(sim_results[[57]][iter,13:14]),
                                         mean = as.numeric(targets[[50]][1,13:14]),
                                         sigma = as.matrix(targets[[50]][2:3,13:14]), log = T)

  ### Death after Graft Loss
  ## Survival-based GOF
  gof_gl_mort_all <- dnorm(as.numeric(sim_results[[49]][iter,]),
                               mean = as.numeric(targets[[51]][1,]),
                               sd = as.numeric(targets[[51]][2,]), log = T)
  gof_gl_mort_diab0 <- dnorm(as.numeric(sim_results[[50]][iter,]),
                                 mean = as.numeric(targets[[52]][1,]),
                                 sd = as.numeric(targets[[52]][2,]), log = T)
  gof_gl_mort_diab1 <- dnorm(as.numeric(sim_results[[51]][iter,]),
                                 mean = as.numeric(targets[[53]][1,]),
                                 sd = as.numeric(targets[[53]][2,]), log = T)
  gof_gl_mort_white <- dnorm(as.numeric(sim_results[[52]][iter,]),
                                 mean = as.numeric(targets[[54]][1,]),
                                 sd = as.numeric(targets[[54]][2,]), log = T)
  gof_gl_mort_black <- dnorm(as.numeric(sim_results[[53]][iter,]),
                                 mean = as.numeric(targets[[55]][1,]),
                                 sd = as.numeric(targets[[55]][2,]), log = T)
  gof_gl_mort_hispanic <- dnorm(as.numeric(sim_results[[54]][iter,]),
                                    mean = as.numeric(targets[[56]][1,]),
                                    sd = as.numeric(targets[[56]][2,]), log = T)
  gof_gl_mort_other <- dnorm(as.numeric(sim_results[[55]][iter,]),
                                 mean = as.numeric(targets[[57]][1,]),
                                 sd = as.numeric(targets[[57]][2,]), log = T)

  ## Shape-based GOF
  gof_gl_mort_shape_all <- dmvnorm(as.numeric(sim_results[[58]][iter,1:2]),
                                       mean = as.numeric(targets[[58]][1,1:2]),
                                       sigma = as.matrix(targets[[58]][2:3,1:2]), log = T)
  gof_gl_mort_shape_diab0 <- dmvnorm(as.numeric(sim_results[[58]][iter,3:4]),
                                         mean = as.numeric(targets[[58]][1,3:4]),
                                         sigma = as.matrix(targets[[58]][2:3,3:4]), log = T)
  gof_gl_mort_shape_diab1 <- dmvnorm(as.numeric(sim_results[[58]][iter,5:6]),
                                         mean = as.numeric(targets[[58]][1,5:6]),
                                         sigma = as.matrix(targets[[58]][2:3,5:6]), log = T)
  gof_gl_mort_shape_white <- dmvnorm(as.numeric(sim_results[[58]][iter,7:8]),
                                         mean = as.numeric(targets[[58]][1,7:8]),
                                         sigma = as.matrix(targets[[58]][2:3,7:8]), log = T)
  gof_gl_mort_shape_black <- dmvnorm(as.numeric(sim_results[[58]][iter,9:10]),
                                         mean = as.numeric(targets[[58]][1,9:10]),
                                         sigma = as.matrix(targets[[58]][2:3,9:10]), log = T)
  gof_gl_mort_shape_hispanic <- dmvnorm(as.numeric(sim_results[[58]][iter,11:12]),
                                            mean = as.numeric(targets[[58]][1,11:12]),
                                            sigma = as.matrix(targets[[58]][2:3,11:12]), log = T)
  gof_gl_mort_shape_other <- dmvnorm(as.numeric(sim_results[[58]][iter,13:14]),
                                         mean = as.numeric(targets[[58]][1,13:14]),
                                         sigma = as.matrix(targets[[58]][2:3,13:14]), log = T)

  ### Death same day as graft loss
  # gof_prop_all      <- dbeta(sim_results[[59]][iter,1], shape1=targets[[59]][1,1], shape2=targets[[59]][1,2], log = T)
  # gof_prop_diab0    <- dbeta(sim_results[[59]][iter,2], shape1=targets[[59]][2,1], shape2=targets[[59]][2,2], log = T)
  # gof_prop_diab1    <- dbeta(sim_results[[59]][iter,3], shape1=targets[[59]][3,1], shape2=targets[[59]][3,2], log = T)
  # gof_prop_white    <- dbeta(sim_results[[59]][iter,4], shape1=targets[[59]][4,1], shape2=targets[[59]][4,2], log = T)
  # gof_prop_black    <- dbeta(sim_results[[59]][iter,5], shape1=targets[[59]][5,1], shape2=targets[[59]][5,2], log = T)
  # gof_prop_hispanic <- dbeta(sim_results[[59]][iter,6], shape1=targets[[59]][6,1], shape2=targets[[59]][6,2], log = T)
  # gof_prop_other    <- dbeta(sim_results[[59]][iter,7], shape1=targets[[59]][7,1], shape2=targets[[59]][7,2], log = T)

  gof_prop_all      <- dbinom(round(as.numeric(sim_results[[59]][iter,1]) * targets[[59]][1,1]),
                              targets[[59]][1,1],
                              targets[[59]][1,1]/(targets[[59]][1,1] + targets[[59]][1,2]),
                              log = T)
  gof_prop_diab0    <- dbinom(round(as.numeric(sim_results[[59]][iter,2]) * targets[[59]][2,1]),
                              targets[[59]][2,1],
                              targets[[59]][2,1]/(targets[[59]][2,1] + targets[[59]][2,2]),
                              log = T)
  gof_prop_diab1    <- dbinom(round(as.numeric(sim_results[[59]][iter,3]) * targets[[59]][3,1]),
                              targets[[59]][3,1],
                              targets[[59]][3,1]/(targets[[59]][3,1] + targets[[59]][3,2]),
                              log = T)
  gof_prop_white    <- dbinom(round(as.numeric(sim_results[[59]][iter,4]) * targets[[59]][4,1]),
                              targets[[59]][4,1],
                              targets[[59]][4,1]/(targets[[59]][4,1] + targets[[59]][4,2]),
                              log = T)
  gof_prop_black    <- dbinom(round(as.numeric(sim_results[[59]][iter,5]) * targets[[59]][5,1]),
                              targets[[59]][5,1],
                              targets[[59]][5,1]/(targets[[59]][5,1] + targets[[59]][5,2]),
                              log = T)
  gof_prop_hispanic <- dbinom(round(as.numeric(sim_results[[59]][iter,6]) * targets[[59]][6,1]),
                              targets[[59]][6,1],
                              targets[[59]][6,1]/(targets[[59]][6,1] + targets[[59]][6,2]),
                              log = T)
  gof_prop_other    <- dbinom(round(as.numeric(sim_results[[59]][iter,7]) * targets[[59]][7,1]),
                              targets[[59]][7,1],
                              targets[[59]][7,1]/(targets[[59]][7,1] + targets[[59]][7,2]),
                              log = T)

  gof_gl_mort_prop <- -1 * sum(gof_prop_all, gof_prop_diab0, gof_prop_diab1,
                          gof_prop_white, gof_prop_black, gof_prop_hispanic, gof_prop_other)

  ### Sum GOF
  # Waitlist
  gof_wait_surv  <- -1 * sum(gof_ddtx_all, gof_ldtx_all, gof_mort_all, gof_remove_all,
                             gof_ddtx_diab0, gof_ldtx_diab0, gof_mort_diab0, gof_remove_diab0,
                             gof_ddtx_diab1, gof_ldtx_diab1, gof_mort_diab1, gof_remove_diab1,
                             gof_ddtx_white, gof_ldtx_white, gof_mort_white, gof_remove_white,
                             gof_ddtx_black, gof_ldtx_black, gof_mort_black, gof_remove_black,
                             gof_ddtx_hispanic, gof_ldtx_hispanic, gof_mort_hispanic, gof_remove_hispanic,
                             gof_ddtx_other, gof_ldtx_other, gof_mort_other, gof_remove_other, na.rm = TRUE)

  gof_wait_shape <- -1 * sum(gof_ddtx_shape_all, gof_ldtx_shape_all, gof_mort_shape_all, gof_remove_shape_all,
                             gof_ddtx_shape_diab0, gof_ldtx_shape_diab0, gof_mort_shape_diab0, gof_remove_shape_diab0,
                             gof_ddtx_shape_diab1, gof_ldtx_shape_diab1, gof_mort_shape_diab1, gof_remove_shape_diab1,
                             gof_ddtx_shape_white, gof_ldtx_shape_white, gof_mort_shape_white, gof_remove_shape_white,
                             gof_ddtx_shape_black, gof_ldtx_shape_black, gof_mort_shape_black, gof_remove_shape_black,
                             gof_ddtx_shape_hispanic, gof_ldtx_shape_hispanic, gof_mort_shape_hispanic, gof_remove_shape_hispanic,
                             gof_ddtx_shape_other, gof_ldtx_shape_other, gof_mort_shape_other, gof_remove_shape_other)

  gof_table1 <- -1 * sum(gof_table1)

  # Post-Transplant
  # Survival
  gof_gs_mort_surv <- -1*sum(gof_gs_mort_all, gof_gs_mort_diab0, gof_gs_mort_diab1,
                             gof_gs_mort_white, gof_gs_mort_black,
                             gof_gs_mort_hispanic, gof_gs_mort_other, na.rm = TRUE)
  gof_gl_surv <- -1*sum(gof_gl_all, gof_gl_diab0, gof_gl_diab1, gof_gl_white,
                        gof_gl_black, gof_gl_hispanic, gof_gl_other, na.rm = TRUE)
  gof_gl_mort_surv <- -1*sum(gof_gl_mort_all, gof_gl_mort_diab0, gof_gl_mort_diab1,
                             gof_gl_mort_white, gof_gl_mort_black,
                             gof_gl_mort_hispanic, gof_gl_mort_other, na.rm = TRUE)

  gof_post_surv  <- -1 * sum(gof_gs_mort_all, gof_gs_mort_diab0, gof_gs_mort_diab1,
                             gof_gs_mort_white, gof_gs_mort_black, gof_gs_mort_hispanic, gof_gs_mort_other,
                             gof_gl_all, gof_gl_diab0, gof_gl_diab1, gof_gl_white, gof_gl_black,
                             gof_gl_hispanic, gof_gl_other, gof_gl_mort_all, gof_gl_mort_diab0, gof_gl_mort_diab1,
                             gof_gl_mort_white, gof_gl_mort_black, gof_gl_mort_hispanic, gof_gl_mort_other, na.rm = TRUE)

  # Shape
  gof_gs_mort_shape <- -1*sum(gof_gs_mort_shape_all, gof_gs_mort_shape_diab0,
                              gof_gs_mort_shape_diab1, gof_gs_mort_shape_white,
                              gof_gs_mort_shape_black, gof_gs_mort_shape_hispanic,
                              gof_gs_mort_shape_other)
  gof_gl_shape <- -1*sum(gof_gl_shape_all, gof_gl_shape_diab0, gof_gl_shape_diab1,
                        gof_gl_shape_white, gof_gl_shape_black,
                        gof_gl_shape_hispanic, gof_gl_shape_other)
  gof_gl_mort_shape <- -1*sum(gof_gl_mort_shape_all, gof_gl_mort_shape_diab0,
                             gof_gl_mort_shape_diab1, gof_gl_mort_shape_white,
                             gof_gl_mort_shape_black, gof_gl_mort_shape_hispanic,
                             gof_gl_mort_shape_other)

  gof_post_shape <- -1 * sum(gof_gs_mort_shape_all, gof_gs_mort_shape_diab0, gof_gs_mort_shape_diab1,
                             gof_gs_mort_shape_white, gof_gs_mort_shape_black, gof_gs_mort_shape_hispanic, gof_gs_mort_shape_other,
                             gof_gl_shape_all, gof_gl_shape_diab0, gof_gl_shape_diab1, gof_gl_shape_white, gof_gl_shape_black,
                             gof_gl_shape_hispanic, gof_gl_shape_other, gof_gl_mort_shape_all, gof_gl_mort_shape_diab0, gof_gl_mort_shape_diab1,
                             gof_gl_mort_shape_white, gof_gl_mort_shape_black, gof_gl_mort_shape_hispanic, gof_gl_mort_shape_other)

  return(list(gof_wait_surv, gof_wait_shape, gof_table1,
              gof_gs_mort_surv, gof_gl_surv, gof_gl_mort_surv,
              gof_gs_mort_shape, gof_gl_shape, gof_gl_mort_shape,
              gof_post_surv, gof_post_shape,
              gof_30, gof_gl_mort_prop))
}


get_likelihood_targets_post <- function(data_gs = data_tx,
                                        data_gl = data_gl_mort) {
  ### Required packages
  require(survival)
  require(flexsurv)

  ### Load Data
  data_30 <- data_gs %>%
    dplyr::filter(test==1)

  data_gs <- data_gs %>%
    dplyr::filter(test==1) %>%
    dplyr::filter(tx_outcome == 0 | tx_outcome == 2) %>%
    mutate(death_at_gl = ifelse(months_to_censor==0 & death==1, 1, 0))

  data_tx_diab <- split(data_gs, data_gs$diab_stat)
  data_tx_race <- split(data_gs, data_gs$race_ethnic)

  data_gl <- data_gl %>%
    dplyr::filter(test==1 & graft_loss==1)  %>%
    mutate(death_at_gl = ifelse(months_to_censor==0 & death==1, 1, 0))

  data_gl_mort_diab <- split(data_gl, data_gl$diab_stat)
  data_gl_mort_race <- split(data_gl, data_gl$race_ethnic)

  data_gl_weibull <- data_gl %>%
    dplyr::filter(test==1 & months_to_censor!=0)

  data_gl_weibull_diab <- split(data_gl_weibull, data_gl_weibull$diab_stat)
  data_gl_weibull_race <- split(data_gl_weibull, data_gl_weibull$race_ethnic)
  ### Multinomial Logit Targets
  alphas_30day       <- matrix(nrow = 7, ncol = 4)
  alphas_30day[1,]   <- as.vector(table(data_30$tx_outcome))
  alphas_30day[2:3,] <- table(data_30$diab_stat, data_30$tx_outcome)
  alphas_30day[4:7,] <- table(data_30$race_ethnic, data_30$tx_outcome)


  ### Initialize matrices to store survival outcome results
  gs_mort_surv          <- matrix(nrow = 2, ncol = length(seq(from=3, to=120, by=3)))
  gs_mort_surv_diab0    <- matrix(nrow = 2, ncol = length(seq(from=3, to=120, by=3)))
  gs_mort_surv_diab1    <- matrix(nrow = 2, ncol = length(seq(from=3, to=120, by=3)))
  gs_mort_surv_white    <- matrix(nrow = 2, ncol = length(seq(from=3, to=120, by=3)))
  gs_mort_surv_black    <- matrix(nrow = 2, ncol = length(seq(from=3, to=120, by=3)))
  gs_mort_surv_hispanic <- matrix(nrow = 2, ncol = length(seq(from=3, to=120, by=3)))
  gs_mort_surv_other    <- matrix(nrow = 2, ncol = length(seq(from=3, to=120, by=3)))
  gs_mort_shape         <- matrix(nrow = 3, ncol = 14)

  gl_surv           <- matrix(nrow = 2, ncol = length(seq(from=3, to=120, by=3)))
  gl_surv_diab0     <- matrix(nrow = 2, ncol = length(seq(from=3, to=120, by=3)))
  gl_surv_diab1     <- matrix(nrow = 2, ncol = length(seq(from=3, to=120, by=3)))
  gl_surv_white     <- matrix(nrow = 2, ncol = length(seq(from=3, to=120, by=3)))
  gl_surv_black     <- matrix(nrow = 2, ncol = length(seq(from=3, to=120, by=3)))
  gl_surv_hispanic  <- matrix(nrow = 2, ncol = length(seq(from=3, to=120, by=3)))
  gl_surv_other     <- matrix(nrow = 2, ncol = length(seq(from=3, to=120, by=3)))
  gl_shape          <- matrix(nrow = 3, ncol = 14)

  gl_mort_surv          <- matrix(nrow = 2, ncol = length(seq(from=3, to=60, by=3)))
  gl_mort_surv_diab0    <- matrix(nrow = 2, ncol = length(seq(from=3, to=60, by=3)))
  gl_mort_surv_diab1    <- matrix(nrow = 2, ncol = length(seq(from=3, to=60, by=3)))
  gl_mort_surv_white    <- matrix(nrow = 2, ncol = length(seq(from=3, to=60, by=3)))
  gl_mort_surv_black    <- matrix(nrow = 2, ncol = length(seq(from=3, to=60, by=3)))
  gl_mort_surv_hispanic <- matrix(nrow = 2, ncol = length(seq(from=3, to=60, by=3)))
  gl_mort_surv_other    <- matrix(nrow = 2, ncol = length(seq(from=3, to=60, by=3)))
  gl_mort_shape         <- matrix(nrow = 3, ncol = 14)

  ### Death with function
  # KM Curves
  gs_mort_fit          <- survival::survfit(Surv(time = months_to_censor, death==1) ~ 1, data = data_gs)
  gs_mort_fit_diab0    <- survival::survfit(Surv(time = months_to_censor, death==1) ~ 1, data = data_tx_diab[[1]])
  gs_mort_fit_diab1    <- survival::survfit(Surv(time = months_to_censor, death==1) ~ 1, data = data_tx_diab[[2]])
  gs_mort_fit_white    <- survival::survfit(Surv(time = months_to_censor, death==1) ~ 1, data = data_tx_race[[1]])
  gs_mort_fit_black    <- survival::survfit(Surv(time = months_to_censor, death==1) ~ 1, data = data_tx_race[[2]])
  gs_mort_fit_hispanic <- survival::survfit(Surv(time = months_to_censor, death==1) ~ 1, data = data_tx_race[[3]])
  gs_mort_fit_other    <- survival::survfit(Surv(time = months_to_censor, death==1) ~ 1, data = data_tx_race[[4]])

  gs_mort_surv[1,]          <- t(summary(gs_mort_fit, times = seq(from=3, to=120, by=3), extend=TRUE)$surv)
  gs_mort_surv_diab0[1,]    <- t(summary(gs_mort_fit_diab0, times = seq(from=3, to=120, by=3), extend=TRUE)$surv)
  gs_mort_surv_diab1[1,]    <- t(summary(gs_mort_fit_diab1, times = seq(from=3, to=120, by=3), extend=TRUE)$surv)
  gs_mort_surv_white[1,]    <- t(summary(gs_mort_fit_white, times = seq(from=3, to=120, by=3), extend=TRUE)$surv)
  gs_mort_surv_black[1,]    <- t(summary(gs_mort_fit_black, times = seq(from=3, to=120, by=3), extend=TRUE)$surv)
  gs_mort_surv_hispanic[1,] <- t(summary(gs_mort_fit_hispanic, times = seq(from=3, to=120, by=3), extend=TRUE)$surv)
  gs_mort_surv_other[1,]    <- t(summary(gs_mort_fit_other, times = seq(from=3, to=120, by=3), extend=TRUE)$surv)

  gs_mort_surv[2,]          <- t(summary(gs_mort_fit, times = seq(from=3, to=120, by=3), extend=TRUE)$std.err)
  gs_mort_surv_diab0[2,]    <- t(summary(gs_mort_fit_diab0, times = seq(from=3, to=120, by=3), extend=TRUE)$std.err)
  gs_mort_surv_diab1[2,]    <- t(summary(gs_mort_fit_diab1, times = seq(from=3, to=120, by=3), extend=TRUE)$std.err)
  gs_mort_surv_white[2,]    <- t(summary(gs_mort_fit_white, times = seq(from=3, to=120, by=3), extend=TRUE)$std.err)
  gs_mort_surv_black[2,]    <- t(summary(gs_mort_fit_black, times = seq(from=3, to=120, by=3), extend=TRUE)$std.err)
  gs_mort_surv_hispanic[2,] <- t(summary(gs_mort_fit_hispanic, times = seq(from=3, to=120, by=3), extend=TRUE)$std.err)
  gs_mort_surv_other[2,]    <- t(summary(gs_mort_fit_other, times = seq(from=3, to=120, by=3), extend=TRUE)$std.err)

  # Parametric Fits
  gs_mort_fit_gompertz          <- flexsurv::flexsurvreg(Surv(time = months_to_censor, death==1) ~ 1, data=data_gs, dist="gompertz")
  gs_mort_fit_gompertz_diab0    <- flexsurv::flexsurvreg(Surv(time = months_to_censor, death==1) ~ 1, data=data_tx_diab[[1]], dist="gompertz")
  gs_mort_fit_gompertz_diab1    <- flexsurv::flexsurvreg(Surv(time = months_to_censor, death==1) ~ 1, data=data_tx_diab[[2]], dist="gompertz")
  gs_mort_fit_gompertz_white    <- flexsurv::flexsurvreg(Surv(time = months_to_censor, death==1) ~ 1, data=data_tx_race[[1]], dist="gompertz")
  gs_mort_fit_gompertz_black    <- flexsurv::flexsurvreg(Surv(time = months_to_censor, death==1) ~ 1, data=data_tx_race[[2]], dist="gompertz")
  gs_mort_fit_gompertz_hispanic <- flexsurv::flexsurvreg(Surv(time = months_to_censor, death==1) ~ 1, data=data_tx_race[[3]], dist="gompertz")
  gs_mort_fit_gompertz_other    <- flexsurv::flexsurvreg(Surv(time = months_to_censor, death==1) ~ 1, data=data_tx_race[[4]], dist="gompertz")

  gs_mort_shape[1,] <- c(gs_mort_fit_gompertz$coefficients,
                         gs_mort_fit_gompertz_diab0$coefficients,
                         gs_mort_fit_gompertz_diab1$coefficients,
                         gs_mort_fit_gompertz_white$coefficients,
                         gs_mort_fit_gompertz_black$coefficients,
                         gs_mort_fit_gompertz_hispanic$coefficients,
                         gs_mort_fit_gompertz_other$coefficients)
  gs_mort_shape[2:3,] <- c(gs_mort_fit_gompertz$cov,
                           gs_mort_fit_gompertz_diab0$cov,
                           gs_mort_fit_gompertz_diab1$cov,
                           gs_mort_fit_gompertz_white$cov,
                           gs_mort_fit_gompertz_black$cov,
                           gs_mort_fit_gompertz_hispanic$cov,
                           gs_mort_fit_gompertz_other$cov)

  ### Graft loss
  gl_fit          <- survival::survfit(Surv(time = months_to_censor, graft_loss==1) ~ 1, data = data_gs)
  gl_fit_diab0    <- survival::survfit(Surv(time = months_to_censor, graft_loss==1) ~ 1, data = data_tx_diab[[1]])
  gl_fit_diab1    <- survival::survfit(Surv(time = months_to_censor, graft_loss==1) ~ 1, data = data_tx_diab[[2]])
  gl_fit_white    <- survival::survfit(Surv(time = months_to_censor, graft_loss==1) ~ 1, data = data_tx_race[[1]])
  gl_fit_black    <- survival::survfit(Surv(time = months_to_censor, graft_loss==1) ~ 1, data = data_tx_race[[2]])
  gl_fit_hispanic <- survival::survfit(Surv(time = months_to_censor, graft_loss==1) ~ 1, data = data_tx_race[[3]])
  gl_fit_other    <- survival::survfit(Surv(time = months_to_censor, graft_loss==1) ~ 1, data = data_tx_race[[4]])

  gl_surv[1,]          <- t(summary(gl_fit, times = seq(from=3, to=120, by=3), extend=TRUE)$surv)
  gl_surv_diab0[1,]    <- t(summary(gl_fit_diab0, times = seq(from=3, to=120, by=3), extend=TRUE)$surv)
  gl_surv_diab1[1,]    <- t(summary(gl_fit_diab1, times = seq(from=3, to=120, by=3), extend=TRUE)$surv)
  gl_surv_white[1,]    <- t(summary(gl_fit_white, times = seq(from=3, to=120, by=3), extend=TRUE)$surv)
  gl_surv_black[1,]    <- t(summary(gl_fit_black, times = seq(from=3, to=120, by=3), extend=TRUE)$surv)
  gl_surv_hispanic[1,] <- t(summary(gl_fit_hispanic, times = seq(from=3, to=120, by=3), extend=TRUE)$surv)
  gl_surv_other[1,]    <- t(summary(gl_fit_other, times = seq(from=3, to=120, by=3), extend=TRUE)$surv)

  gl_surv[2,]          <- t(summary(gl_fit, times = seq(from=3, to=120, by=3), extend=TRUE)$std.err)
  gl_surv_diab0[2,]    <- t(summary(gl_fit_diab0, times = seq(from=3, to=120, by=3), extend=TRUE)$std.err)
  gl_surv_diab1[2,]    <- t(summary(gl_fit_diab1, times = seq(from=3, to=120, by=3), extend=TRUE)$std.err)
  gl_surv_white[2,]    <- t(summary(gl_fit_white, times = seq(from=3, to=120, by=3), extend=TRUE)$std.err)
  gl_surv_black[2,]    <- t(summary(gl_fit_black, times = seq(from=3, to=120, by=3), extend=TRUE)$std.err)
  gl_surv_hispanic[2,] <- t(summary(gl_fit_hispanic, times = seq(from=3, to=120, by=3), extend=TRUE)$std.err)
  gl_surv_other[2,]    <- t(summary(gl_fit_other, times = seq(from=3, to=120, by=3), extend=TRUE)$std.err)

  # Parametric Fits
  gl_fit_gompertz          <- flexsurv::flexsurvreg(Surv(time = months_to_censor, graft_loss==1) ~ 1, data=data_gs, dist="gompertz")
  gl_fit_gompertz_diab0    <- flexsurv::flexsurvreg(Surv(time = months_to_censor, graft_loss==1) ~ 1, data=data_tx_diab[[1]], dist="gompertz")
  gl_fit_gompertz_diab1    <- flexsurv::flexsurvreg(Surv(time = months_to_censor, graft_loss==1) ~ 1, data=data_tx_diab[[2]], dist="gompertz")
  gl_fit_gompertz_white    <- flexsurv::flexsurvreg(Surv(time = months_to_censor, graft_loss==1) ~ 1, data=data_tx_race[[1]], dist="gompertz")
  gl_fit_gompertz_black    <- flexsurv::flexsurvreg(Surv(time = months_to_censor, graft_loss==1) ~ 1, data=data_tx_race[[2]], dist="gompertz")
  gl_fit_gompertz_hispanic <- flexsurv::flexsurvreg(Surv(time = months_to_censor, graft_loss==1) ~ 1, data=data_tx_race[[3]], dist="gompertz")
  gl_fit_gompertz_other    <- flexsurv::flexsurvreg(Surv(time = months_to_censor, graft_loss==1) ~ 1, data=data_tx_race[[4]], dist="gompertz")

  gl_shape[1,] <- c(gl_fit_gompertz$coefficients,
                    gl_fit_gompertz_diab0$coefficients,
                    gl_fit_gompertz_diab1$coefficients,
                    gl_fit_gompertz_white$coefficients,
                    gl_fit_gompertz_black$coefficients,
                    gl_fit_gompertz_hispanic$coefficients,
                    gl_fit_gompertz_other$coefficients)
  gl_shape[2:3,] <- c(gl_fit_gompertz$cov,
                      gl_fit_gompertz_diab0$cov,
                      gl_fit_gompertz_diab1$cov,
                      gl_fit_gompertz_white$cov,
                      gl_fit_gompertz_black$cov,
                      gl_fit_gompertz_hispanic$cov,
                      gl_fit_gompertz_other$cov)

  # Death after graft loss
  gl_mort_fit          <- survival::survfit(Surv(time = months_to_censor, death==1) ~ 1, data = data_gl)
  gl_mort_fit_diab0    <- survival::survfit(Surv(time = months_to_censor, death==1) ~ 1, data = data_gl_mort_diab[[1]])
  gl_mort_fit_diab1    <- survival::survfit(Surv(time = months_to_censor, death==1) ~ 1, data = data_gl_mort_diab[[2]])
  gl_mort_fit_white    <- survival::survfit(Surv(time = months_to_censor, death==1) ~ 1, data = data_gl_mort_race[[1]])
  gl_mort_fit_black    <- survival::survfit(Surv(time = months_to_censor, death==1) ~ 1, data = data_gl_mort_race[[2]])
  gl_mort_fit_hispanic <- survival::survfit(Surv(time = months_to_censor, death==1) ~ 1, data = data_gl_mort_race[[3]])
  gl_mort_fit_other    <- survival::survfit(Surv(time = months_to_censor, death==1) ~ 1, data = data_gl_mort_race[[4]])

  gl_mort_surv[1,]          <- t(summary(gl_mort_fit, times = seq(from=3, to=60, by=3), extend=TRUE)$surv)
  gl_mort_surv_diab0[1,]    <- t(summary(gl_mort_fit_diab0, times = seq(from=3, to=60, by=3), extend=TRUE)$surv)
  gl_mort_surv_diab1[1,]    <- t(summary(gl_mort_fit_diab1, times = seq(from=3, to=60, by=3), extend=TRUE)$surv)
  gl_mort_surv_white[1,]    <- t(summary(gl_mort_fit_white, times = seq(from=3, to=60, by=3), extend=TRUE)$surv)
  gl_mort_surv_black[1,]    <- t(summary(gl_mort_fit_black, times = seq(from=3, to=60, by=3), extend=TRUE)$surv)
  gl_mort_surv_hispanic[1,] <- t(summary(gl_mort_fit_hispanic, times = seq(from=3, to=60, by=3), extend=TRUE)$surv)
  gl_mort_surv_other[1,]    <- t(summary(gl_mort_fit_other, times = seq(from=3, to=60, by=3), extend=TRUE)$surv)

  gl_mort_surv[2,]          <- t(summary(gl_mort_fit, times = seq(from=3, to=60, by=3), extend=TRUE)$std.err)
  gl_mort_surv_diab0[2,]    <- t(summary(gl_mort_fit_diab0, times = seq(from=3, to=60, by=3), extend=TRUE)$std.err)
  gl_mort_surv_diab1[2,]    <- t(summary(gl_mort_fit_diab1, times = seq(from=3, to=60, by=3), extend=TRUE)$std.err)
  gl_mort_surv_white[2,]    <- t(summary(gl_mort_fit_white, times = seq(from=3, to=60, by=3), extend=TRUE)$std.err)
  gl_mort_surv_black[2,]    <- t(summary(gl_mort_fit_black, times = seq(from=3, to=60, by=3), extend=TRUE)$std.err)
  gl_mort_surv_hispanic[2,] <- t(summary(gl_mort_fit_hispanic, times = seq(from=3, to=60, by=3), extend=TRUE)$std.err)
  gl_mort_surv_other[2,]    <- t(summary(gl_mort_fit_other, times = seq(from=3, to=60, by=3), extend=TRUE)$std.err)

  # Parametric Fits
  gl_mort_fit_gompertz          <- flexsurv::flexsurvreg(Surv(time = months_to_censor, death==1) ~ 1, data=data_gl_weibull, dist="weibull")
  gl_mort_fit_gompertz_diab0    <- flexsurv::flexsurvreg(Surv(time = months_to_censor, death==1) ~ 1, data=data_gl_weibull_diab[[1]], dist="weibull")
  gl_mort_fit_gompertz_diab1    <- flexsurv::flexsurvreg(Surv(time = months_to_censor, death==1) ~ 1, data=data_gl_weibull_diab[[2]], dist="weibull")
  gl_mort_fit_gompertz_white    <- flexsurv::flexsurvreg(Surv(time = months_to_censor, death==1) ~ 1, data=data_gl_weibull_race[[1]], dist="weibull")
  gl_mort_fit_gompertz_black    <- flexsurv::flexsurvreg(Surv(time = months_to_censor, death==1) ~ 1, data=data_gl_weibull_race[[2]], dist="weibull")
  gl_mort_fit_gompertz_hispanic <- flexsurv::flexsurvreg(Surv(time = months_to_censor, death==1) ~ 1, data=data_gl_weibull_race[[3]], dist="weibull")
  gl_mort_fit_gompertz_other    <- flexsurv::flexsurvreg(Surv(time = months_to_censor, death==1) ~ 1, data=data_gl_weibull_race[[4]], dist="weibull")

  gl_mort_shape[1,] <- c(gl_mort_fit_gompertz$coefficients,
                         gl_mort_fit_gompertz_diab0$coefficients,
                         gl_mort_fit_gompertz_diab1$coefficients,
                         gl_mort_fit_gompertz_white$coefficients,
                         gl_mort_fit_gompertz_black$coefficients,
                         gl_mort_fit_gompertz_hispanic$coefficients,
                         gl_mort_fit_gompertz_other$coefficients)
  gl_mort_shape[2:3,] <- c(gl_mort_fit_gompertz$cov,
                           gl_mort_fit_gompertz_diab0$cov,
                           gl_mort_fit_gompertz_diab1$cov,
                           gl_mort_fit_gompertz_white$cov,
                           gl_mort_fit_gompertz_black$cov,
                           gl_mort_fit_gompertz_hispanic$cov,
                           gl_mort_fit_gompertz_other$cov)

  ### Death same day as graft loss
  prop_gl_mort <- matrix(nrow = 7, ncol = 2,
                          dimnames = list(c("All", "Diab0", "Diab1", "White",
                                            "Black", "Hispanic", "Other"),
                                          c("Shape1", "Shape2")))
  prop_gl_mort[1,]   <- rev(as.vector(table(data_gl$death_at_gl)))
  prop_gl_mort[2:3,] <- table(data_gl$diab_stat, data_gl$death_at_gl)[,2:1]
  prop_gl_mort[4:7,] <- table(data_gl$race_ethnic, data_gl$death_at_gl)[,2:1]

  ### Return results
  return(list(alphas_30day, gs_mort_surv, gs_mort_surv_diab0,
              gs_mort_surv_diab1, gs_mort_surv_white, gs_mort_surv_black,
              gs_mort_surv_hispanic, gs_mort_surv_other, gs_mort_shape,
              gl_surv, gl_surv_diab0, gl_surv_diab1, gl_surv_white,
              gl_surv_black, gl_surv_hispanic, gl_surv_other,
              gl_shape, gl_mort_surv, gl_mort_surv_diab0,
              gl_mort_surv_diab1, gl_mort_surv_white, gl_mort_surv_black,
              gl_mort_surv_hispanic, gl_mort_surv_other, gl_mort_shape,
              prop_gl_mort))
}


tx_targets <- get_likelihood_targets_post()

for(i in 1:length(tx_targets)){
  write.csv(tx_targets[i], paste0(out_path, "Recalibration/Stage2/target_", i+33, ".csv"))
}

