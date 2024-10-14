## Load packages
library(doParallel)
library(doRNG)
library(foreach)
library(dplyr)
library(doSNOW)
library(data.table)


n <- 100000

ddtx_coef <- numeric()
ldtx_coef <- numeric()
mort_coef <- numeric()
remove_coef <- numeric()
ddtx_shape <- numeric()
ldtx_shape <- numeric()
mort_shape <- numeric()
remove_shape <- numeric()

for(i in 1:n) {
  param_set <- psa_list()

  ## Store coefficients for waitlist functions
  ddtx_coef   <- rbind(ddtx_coef, t(param_set$ddtx_coef_psa))
  ldtx_coef   <- rbind(ldtx_coef, t(param_set$ldtx_coef_psa))
  mort_coef   <- rbind(mort_coef, t(param_set$mort_coef_psa))
  remove_coef <- rbind(remove_coef, t(param_set$remove_coef_psa))

  ddtx_shape     <- rbind(ddtx_shape, param_set$ddtx_shape_psa)
  ldtx_shape     <- rbind(ldtx_shape, param_set$ldtx_shape_psa)
  mort_shape     <- rbind(mort_shape, param_set$mort_shape_psa)
  remove_shape   <- rbind(remove_shape, param_set$remove_shape_psa)
}

## Set up to run in parallel
comb <- function(...) {
  mapply('rbind', ..., SIMPLIFY=FALSE)
}

parallel::detectCores()
n.cores <- 50

my.cluster <- parallel::makeCluster(
  n.cores,
  type = "PSOCK"
  )


registerDoSNOW(my.cluster)

seed_num <- sample(1:999999, 1)
# doParallel::registerDoParallel(cl = my.cluster)
doRNG::registerDoRNG(seed = seed_num)

pb <- txtProgressBar(min=1, max=n, style=3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress=progress)

out <- foreach(i = 1:n, .combine='comb', .multicombine=TRUE, .options.snow=opts,
               .packages = c("MASS", "copula", "psych", "dplyr", "survival",
                             "flexsurv")) %dopar% {
                 ## Run Waitlist Simulation
                               list_sim_df <- list.simulation(n=100000,
                                                              seed = round(runif(1, min=1, max=99999)),
                                 data = sample_df,
                                 ddtx_coef = ddtx_coef[i,],
                                 ddtx_p = ddtx_shape[i,1],
                                 ldtx_coef = ldtx_coef[i,],
                                 ldtx_gamma = ldtx_shape[i,1],
                                 mort_coef = mort_coef[i,],
                                 mort_shape = mort_shape[i,1],
                                 remove_coef = remove_coef[i,],
                                 remove_p = remove_shape[i,1])

                 ## Clean and run survival models
                 list_sim_df <- list_sim_df %>%
                   mutate(rec_age_at_tx_c = can_age_at_listing_c + (months_to_event/12)) %>%
                   mutate(years_dial = baseline_yrs_dial + (months_to_event/12)) %>%
                   mutate(tx_year_c = list_year_c + (months_to_event/12)) %>%
                   mutate(tx_year_c = ifelse(tx_year_c > 10, 10, tx_year_c))

                 ### All Candidates
                 # Kaplan-Meier
                 ddtx_km_fit    <- survfit(Surv(time = months_to_event, ddtx==1) ~ 1, data = list_sim_df)
                 ldtx_km_fit    <- survfit(Surv(time = months_to_event, ldtx==1) ~ 1, data = list_sim_df)
                 mort_km_fit    <- survfit(Surv(time = months_to_event, mort==1) ~ 1, data = list_sim_df)
                 remove_km_fit  <- survfit(Surv(time = months_to_event, remove==1) ~ 1, data = list_sim_df)

                 ddtx_surv   <- summary(ddtx_km_fit, times = seq(from=3, to=180, by=3), extend=TRUE)$surv
                 ldtx_surv   <- summary(ldtx_km_fit, times = seq(from=3, to=180, by=3), extend=TRUE)$surv
                 mort_surv   <- summary(mort_km_fit, times = seq(from=3, to=180, by=3), extend=TRUE)$surv
                 remove_surv <- summary(remove_km_fit, times = seq(from=3, to=180, by=3), extend=TRUE)$surv

                 # Parametric Fits
                 ddtx_weibull_fit = flexsurvreg(Surv(time = months_to_event, ddtx==1)~1,
                                                data=list_sim_df, dist="weibull")
                 ldtx_gompertz_fit = flexsurvreg(Surv(time = months_to_event, ldtx==1)~1,
                                                 data=list_sim_df, dist="gompertz")
                 mort_weibull_fit = flexsurvreg(Surv(time = months_to_event, mort==1)~1,
                                                data=list_sim_df, dist="weibull")
                 remove_weibull_fit = flexsurvreg(Surv(time = months_to_event, remove==1)~1,
                                                  data=list_sim_df, dist="weibull")

                 ddtx_shape_all <- ddtx_weibull_fit$coefficients

                 ldtx_shape_all <- ldtx_gompertz_fit$coefficients

                 mort_shape_all <- mort_weibull_fit$coefficients

                 remove_shape_all <- remove_weibull_fit$coefficients

                 ### By Diabetes Status
                 df_diab <- split(list_sim_df, list_sim_df$diab_stat)

                 ## No Diabetes
                 # Kaplan-Meier
                 ddtx_km_fit_diab0    <- survfit(Surv(time = months_to_event, ddtx==1) ~ 1, data = df_diab[[1]])
                 ldtx_km_fit_diab0    <- survfit(Surv(time = months_to_event, ldtx==1) ~ 1, data = df_diab[[1]])
                 mort_km_fit_diab0    <- survfit(Surv(time = months_to_event, mort==1) ~ 1, data = df_diab[[1]])
                 remove_km_fit_diab0  <- survfit(Surv(time = months_to_event, remove==1) ~ 1, data = df_diab[[1]])

                 ddtx_surv_diab0   <- summary(ddtx_km_fit_diab0, times = seq(from=3, to=180, by=3), extend=TRUE)$surv
                 ldtx_surv_diab0   <- summary(ldtx_km_fit_diab0, times = seq(from=3, to=180, by=3), extend=TRUE)$surv
                 mort_surv_diab0   <- summary(mort_km_fit_diab0, times = seq(from=3, to=180, by=3), extend=TRUE)$surv
                 remove_surv_diab0 <- summary(remove_km_fit_diab0, times = seq(from=3, to=180, by=3), extend=TRUE)$surv

                 # Parametric Fits
                 ddtx_weibull_fit_diab0 = flexsurvreg(Surv(time = months_to_event, ddtx==1)~1,
                                                data=df_diab[[1]], dist="weibull")
                 ldtx_gompertz_fit_diab0 = flexsurvreg(Surv(time = months_to_event, ldtx==1)~1,
                                                 data=df_diab[[1]], dist="gompertz")
                 mort_weibull_fit_diab0 = flexsurvreg(Surv(time = months_to_event, mort==1)~1,
                                                data=df_diab[[1]], dist="weibull")
                 remove_weibull_fit_diab0 = flexsurvreg(Surv(time = months_to_event, remove==1)~1,
                                                  data=df_diab[[1]], dist="weibull")

                 ddtx_shape_diab0 <- ddtx_weibull_fit_diab0$coefficients

                 ldtx_shape_diab0 <- ldtx_gompertz_fit_diab0$coefficients

                 mort_shape_diab0 <- mort_weibull_fit_diab0$coefficients

                 remove_shape_diab0 <- remove_weibull_fit_diab0$coefficients

                 ## History of Diabetes
                 # Kaplan-Meier
                 ddtx_km_fit_diab1    <- survfit(Surv(time = months_to_event, ddtx==1) ~ 1, data = df_diab[[2]])
                 ldtx_km_fit_diab1    <- survfit(Surv(time = months_to_event, ldtx==1) ~ 1, data = df_diab[[2]])
                 mort_km_fit_diab1    <- survfit(Surv(time = months_to_event, mort==1) ~ 1, data = df_diab[[2]])
                 remove_km_fit_diab1  <- survfit(Surv(time = months_to_event, remove==1) ~ 1, data = df_diab[[2]])

                 ddtx_surv_diab1   <- summary(ddtx_km_fit_diab1, times = seq(from=3, to=180, by=3), extend=TRUE)$surv
                 ldtx_surv_diab1   <- summary(ldtx_km_fit_diab1, times = seq(from=3, to=180, by=3), extend=TRUE)$surv
                 mort_surv_diab1   <- summary(mort_km_fit_diab1, times = seq(from=3, to=180, by=3), extend=TRUE)$surv
                 remove_surv_diab1 <- summary(remove_km_fit_diab1, times = seq(from=3, to=180, by=3), extend=TRUE)$surv

                 # Parametric Fits
                 ddtx_weibull_fit_diab1 = flexsurvreg(Surv(time = months_to_event, ddtx==1)~1,
                                                      data=df_diab[[2]], dist="weibull")
                 ldtx_gompertz_fit_diab1 = flexsurvreg(Surv(time = months_to_event, ldtx==1)~1,
                                                       data=df_diab[[2]], dist="gompertz")
                 mort_weibull_fit_diab1 = flexsurvreg(Surv(time = months_to_event, mort==1)~1,
                                                      data=df_diab[[2]], dist="weibull")
                 remove_weibull_fit_diab1 = flexsurvreg(Surv(time = months_to_event, remove==1)~1,
                                                        data=df_diab[[2]], dist="weibull")

                 ddtx_shape_diab1 <- ddtx_weibull_fit_diab1$coefficients

                 ldtx_shape_diab1 <- ldtx_gompertz_fit_diab1$coefficients

                 mort_shape_diab1 <- mort_weibull_fit_diab1$coefficients

                 remove_shape_diab1 <- remove_weibull_fit_diab1$coefficients


                 ### By Race/Ethnicity
                 df_race <- split(list_sim_df, list_sim_df$race)

                 ## White Candidates
                 # Kaplan-Meier
                 ddtx_km_fit_white    <- survfit(Surv(time = months_to_event, ddtx==1) ~ 1, data = df_race[[1]])
                 ldtx_km_fit_white    <- survfit(Surv(time = months_to_event, ldtx==1) ~ 1, data = df_race[[1]])
                 mort_km_fit_white    <- survfit(Surv(time = months_to_event, mort==1) ~ 1, data = df_race[[1]])
                 remove_km_fit_white  <- survfit(Surv(time = months_to_event, remove==1) ~ 1, data = df_race[[1]])

                 ddtx_surv_white   <- summary(ddtx_km_fit_white, times = seq(from=3, to=180, by=3), extend=TRUE)$surv
                 ldtx_surv_white   <- summary(ldtx_km_fit_white, times = seq(from=3, to=180, by=3), extend=TRUE)$surv
                 mort_surv_white   <- summary(mort_km_fit_white, times = seq(from=3, to=180, by=3), extend=TRUE)$surv
                 remove_surv_white <- summary(remove_km_fit_white, times = seq(from=3, to=180, by=3), extend=TRUE)$surv

                 # Parametric Fits
                 ddtx_weibull_fit_white = flexsurvreg(Surv(time = months_to_event, ddtx==1)~1,
                                                      data=df_race[[1]], dist="weibull")
                 ldtx_gompertz_fit_white = flexsurvreg(Surv(time = months_to_event, ldtx==1)~1,
                                                       data=df_race[[1]], dist="gompertz")
                 mort_weibull_fit_white = flexsurvreg(Surv(time = months_to_event, mort==1)~1,
                                                      data=df_race[[1]], dist="weibull")
                 remove_weibull_fit_white = flexsurvreg(Surv(time = months_to_event, remove==1)~1,
                                                        data=df_race[[1]], dist="weibull")

                 ddtx_shape_white <- ddtx_weibull_fit_white$coefficients

                 ldtx_shape_white <- ldtx_gompertz_fit_white$coefficients

                 mort_shape_white <- mort_weibull_fit_white$coefficients

                 remove_shape_white <- remove_weibull_fit_white$coefficients

                 ## Black Candidates
                 # Kaplan-Meier
                 ddtx_km_fit_black    <- survfit(Surv(time = months_to_event, ddtx==1) ~ 1, data = df_race[[2]])
                 ldtx_km_fit_black    <- survfit(Surv(time = months_to_event, ldtx==1) ~ 1, data = df_race[[2]])
                 mort_km_fit_black    <- survfit(Surv(time = months_to_event, mort==1) ~ 1, data = df_race[[2]])
                 remove_km_fit_black  <- survfit(Surv(time = months_to_event, remove==1) ~ 1, data = df_race[[2]])

                 ddtx_surv_black   <- summary(ddtx_km_fit_black, times = seq(from=3, to=180, by=3), extend=TRUE)$surv
                 ldtx_surv_black   <- summary(ldtx_km_fit_black, times = seq(from=3, to=180, by=3), extend=TRUE)$surv
                 mort_surv_black   <- summary(mort_km_fit_black, times = seq(from=3, to=180, by=3), extend=TRUE)$surv
                 remove_surv_black <- summary(remove_km_fit_black, times = seq(from=3, to=180, by=3), extend=TRUE)$surv

                 # Parametric Fits
                 ddtx_weibull_fit_black = flexsurvreg(Surv(time = months_to_event, ddtx==1)~1,
                                                      data=df_race[[2]], dist="weibull")
                 ldtx_gompertz_fit_black = flexsurvreg(Surv(time = months_to_event, ldtx==1)~1,
                                                       data=df_race[[2]], dist="gompertz")
                 mort_weibull_fit_black = flexsurvreg(Surv(time = months_to_event, mort==1)~1,
                                                      data=df_race[[2]], dist="weibull")
                 remove_weibull_fit_black = flexsurvreg(Surv(time = months_to_event, remove==1)~1,
                                                        data=df_race[[2]], dist="weibull")

                 ddtx_shape_black <- ddtx_weibull_fit_black$coefficients

                 ldtx_shape_black <- ldtx_gompertz_fit_black$coefficients

                 mort_shape_black <- mort_weibull_fit_black$coefficients

                 remove_shape_black <- remove_weibull_fit_black$coefficients

                 ## Hispanic Candidates
                 # Kaplan-Meier
                 ddtx_km_fit_hispanic    <- survfit(Surv(time = months_to_event, ddtx==1) ~ 1, data = df_race[[3]])
                 ldtx_km_fit_hispanic    <- survfit(Surv(time = months_to_event, ldtx==1) ~ 1, data = df_race[[3]])
                 mort_km_fit_hispanic    <- survfit(Surv(time = months_to_event, mort==1) ~ 1, data = df_race[[3]])
                 remove_km_fit_hispanic  <- survfit(Surv(time = months_to_event, remove==1) ~ 1, data = df_race[[3]])

                 ddtx_surv_hispanic   <- summary(ddtx_km_fit_hispanic, times = seq(from=3, to=180, by=3), extend=TRUE)$surv
                 ldtx_surv_hispanic   <- summary(ldtx_km_fit_hispanic, times = seq(from=3, to=180, by=3), extend=TRUE)$surv
                 mort_surv_hispanic   <- summary(mort_km_fit_hispanic, times = seq(from=3, to=180, by=3), extend=TRUE)$surv
                 remove_surv_hispanic <- summary(remove_km_fit_hispanic, times = seq(from=3, to=180, by=3), extend=TRUE)$surv

                 # Parametric Fits
                 ddtx_weibull_fit_hispanic = flexsurvreg(Surv(time = months_to_event, ddtx==1)~1,
                                                      data=df_race[[3]], dist="weibull")
                 ldtx_gompertz_fit_hispanic = flexsurvreg(Surv(time = months_to_event, ldtx==1)~1,
                                                       data=df_race[[3]], dist="gompertz")
                 mort_weibull_fit_hispanic = flexsurvreg(Surv(time = months_to_event, mort==1)~1,
                                                      data=df_race[[3]], dist="weibull")
                 remove_weibull_fit_hispanic = flexsurvreg(Surv(time = months_to_event, remove==1)~1,
                                                        data=df_race[[3]], dist="weibull")

                 ddtx_shape_hispanic <- ddtx_weibull_fit_hispanic$coefficients

                 ldtx_shape_hispanic <- ldtx_gompertz_fit_hispanic$coefficients

                 mort_shape_hispanic <- mort_weibull_fit_hispanic$coefficients

                 remove_shape_hispanic <- remove_weibull_fit_hispanic$coefficients

                 ## Other Candidates
                 # Kaplan-Meier
                 ddtx_km_fit_other    <- survfit(Surv(time = months_to_event, ddtx==1) ~ 1, data = df_race[[4]])
                 ldtx_km_fit_other    <- survfit(Surv(time = months_to_event, ldtx==1) ~ 1, data = df_race[[4]])
                 mort_km_fit_other    <- survfit(Surv(time = months_to_event, mort==1) ~ 1, data = df_race[[4]])
                 remove_km_fit_other  <- survfit(Surv(time = months_to_event, remove==1) ~ 1, data = df_race[[4]])

                 ddtx_surv_other   <- summary(ddtx_km_fit_other, times = seq(from=3, to=180, by=3), extend=TRUE)$surv
                 ldtx_surv_other   <- summary(ldtx_km_fit_other, times = seq(from=3, to=180, by=3), extend=TRUE)$surv
                 mort_surv_other   <- summary(mort_km_fit_other, times = seq(from=3, to=180, by=3), extend=TRUE)$surv
                 remove_surv_other <- summary(remove_km_fit_other, times = seq(from=3, to=180, by=3), extend=TRUE)$surv

                 # Parametric Fits
                 ddtx_weibull_fit_other = flexsurvreg(Surv(time = months_to_event, ddtx==1)~1,
                                                      data=df_race[[4]], dist="weibull")
                 ldtx_gompertz_fit_other = flexsurvreg(Surv(time = months_to_event, ldtx==1)~1,
                                                       data=df_race[[4]], dist="gompertz")
                 mort_weibull_fit_other = flexsurvreg(Surv(time = months_to_event, mort==1)~1,
                                                      data=df_race[[4]], dist="weibull")
                 remove_weibull_fit_other = flexsurvreg(Surv(time = months_to_event, remove==1)~1,
                                                        data=df_race[[4]], dist="weibull")

                 ddtx_shape_other <- ddtx_weibull_fit_other$coefficients

                 ldtx_shape_other <- ldtx_gompertz_fit_other$coefficients

                 mort_shape_other <- mort_weibull_fit_other$coefficients

                 remove_shape_other <- remove_weibull_fit_other$coefficients

                 ### Table 1
                 df_tx <- list_sim_df %>%
                   mutate(can_blood_a = ifelse(can_blood_ab==0 & can_blood_b==0 & can_blood_o==0,
                                               1, 0)) %>%
                   mutate(rec_age_at_tx = rec_age_at_tx_c + 65) %>%
                   mutate(optn_reg1 = ifelse(optn_reg2==0 & optn_reg3==0 & optn_reg4==0 &
                                               optn_reg5==0 & optn_reg6==0 & optn_reg7==0 &
                                               optn_reg8==0 & optn_reg9==0 & optn_reg10==0 &
                                               optn_reg11==0, 1, 0)) %>%
                   mutate(white = ifelse(black==0 & hispanic==0 & other_race==0, 1, 0)) %>%
                   filter(ddtx==1)

                 sim_table1 <- matrix(data = NA, nrow = 1, ncol = 27)
                 sim_table1[,1] <- mean(df_tx$rec_age_at_tx)
                 sim_table1[,2] <- mean(df_tx$male)
                 sim_table1[,3] <- mean(df_tx$white)
                 sim_table1[,4] <- mean(df_tx$black)
                 sim_table1[,5] <- mean(df_tx$hispanic)
                 sim_table1[,6] <- mean(df_tx$other_race)
                 sim_table1[,7] <- mean(df_tx$can_blood_a)
                 sim_table1[,8] <- mean(df_tx$can_blood_ab)
                 sim_table1[,9] <- mean(df_tx$can_blood_b)
                 sim_table1[,10] <- mean(df_tx$can_blood_o)
                 sim_table1[,11] <- mean(df_tx$years_dial)
                 sim_table1[,12] <- mean(df_tx$diab_stat)
                 sim_table1[,13] <- mean(df_tx$copd)
                 sim_table1[,14] <- mean(df_tx$pvd)
                 sim_table1[,15] <- mean(df_tx$ang_cad)
                 sim_table1[,16] <- mean(df_tx$canhx_cpra)
                 sim_table1[,17] <- mean(df_tx$optn_reg1)
                 sim_table1[,18] <- mean(df_tx$optn_reg2)
                 sim_table1[,19] <- mean(df_tx$optn_reg3)
                 sim_table1[,20] <- mean(df_tx$optn_reg4)
                 sim_table1[,21] <- mean(df_tx$optn_reg5)
                 sim_table1[,22] <- mean(df_tx$optn_reg6)
                 sim_table1[,23] <- mean(df_tx$optn_reg7)
                 sim_table1[,24] <- mean(df_tx$optn_reg8)
                 sim_table1[,25] <- mean(df_tx$optn_reg9)
                 sim_table1[,26] <- mean(df_tx$optn_reg10)
                 sim_table1[,27] <- mean(df_tx$optn_reg11)

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

                 ### Export
                 list(t(ddtx_surv),
                      t(ldtx_surv),
                      t(mort_surv),
                      t(remove_surv),
                      t(ddtx_surv_diab0),
                      t(ldtx_surv_diab0),
                      t(mort_surv_diab0),
                      t(remove_surv_diab0),
                      t(ddtx_surv_diab1),
                      t(ldtx_surv_diab1),
                      t(mort_surv_diab1),
                      t(remove_surv_diab1),
                      t(ddtx_surv_white),
                      t(ldtx_surv_white),
                      t(mort_surv_white),
                      t(remove_surv_white),
                      t(ddtx_surv_black),
                      t(ldtx_surv_black),
                      t(mort_surv_black),
                      t(remove_surv_black),
                      t(ddtx_surv_hispanic),
                      t(ldtx_surv_hispanic),
                      t(mort_surv_hispanic),
                      t(remove_surv_hispanic),
                      t(ddtx_surv_other),
                      t(ldtx_surv_other),
                      t(mort_surv_other),
                      t(remove_surv_other),
                      m_ddtx_shape,
                      m_ldtx_shape,
                      m_mort_shape,
                      m_remove_shape,
                      sim_table1)
               }
parallel::stopCluster(cl = my.cluster)

for(i in 1:length(out)){
  write.csv(out[i], paste0(out_path, "Recalibration/waitlist/out_", i, ".csv"))
}

id <- seq(from = 1, to = n, by = 1)
m_ddtx_coef <- cbind(id, ddtx_coef, ddtx_shape)
m_ldtx_coef <- cbind(id, ldtx_coef, ldtx_shape)
m_mort_coef <- cbind(id, mort_coef, mort_shape)
m_remove_coef <- cbind(id, remove_coef, remove_shape)

write.csv(m_ddtx_coef, paste0(out_path, "Recalibration/waitlist/ddtx_coef.csv"))
write.csv(m_ldtx_coef, paste0(out_path, "Recalibration/waitlist/ldtx_coef.csv"))
write.csv(m_mort_coef, paste0(out_path, "Recalibration/waitlist/mort_coef.csv"))
write.csv(m_remove_coef, paste0(out_path, "Recalibration/waitlist/remove_coef.csv"))


l_target_names <- list("target_1.csv", "target_2.csv", "target_3.csv",
                       "target_4.csv", "target_5.csv", "target_6.csv",
                       "target_7.csv", "target_8.csv", "target_9.csv",
                       "target_10.csv", "target_11.csv", "target_12.csv",
                       "target_13.csv", "target_14.csv", "target_15.csv",
                       "target_16.csv", "target_17.csv", "target_18.csv",
                       "target_19.csv", "target_20.csv", "target_21.csv",
                       "target_22.csv", "target_23.csv", "target_24.csv",
                       "target_25.csv", "target_26.csv", "target_27.csv",
                       "target_28.csv", "target_29.csv", "target_30.csv",
                       "target_31.csv", "target_32.csv", "target_33.csv")

l_target <- list()
for (k in 1:length(l_target_names)){
  l_target[[k]] <- read.csv(paste0(out_path, "Recalibration/waitlist/", l_target_names[k]))[-1]
}

m_gof <- matrix(nrow = n, ncol = 4)
pb <- txtProgressBar(min = 0, max = n, style = 3)
for (i in 1:n) {
  m_gof[i,1] <- i
  l_gof <- get_gof_likelihood(targets = l_target, sim_results = out, iter = i)
  m_gof[i,2] <- l_gof[[1]]
  m_gof[i,3] <- l_gof[[2]]
  m_gof[i,4] <- l_gof[[3]]
  setTxtProgressBar(pb, i)
}

colnames(m_gof) <- c("id", "nll_surv", "nll_shape", "nll_table1")
write.csv(m_gof, paste0(out_path, "Recalibration/waitlist/m_gof.csv"))



df_gof <- as.data.frame(m_gof) %>%
  mutate(nll_surv_norm = nll_surv/1680) %>%
  mutate(nll_shape_norm = nll_shape/28) %>%
  mutate(gof_surv = nll_surv_norm + nll_table1) %>%
  mutate(gof_shape = nll_shape_norm + nll_table1) %>%
  mutate(rank_surv = rank(gof_surv)) %>%
  mutate(rank_shape = rank(gof_shape)) %>%
  mutate(w_surv = 1/exp(gof_surv)) %>%
  mutate(w_shape = 1/exp(gof_shape)) %>%
  mutate(w_surv = w_surv/sum(w_surv)) %>%
  mutate(w_shape = w_shape/sum(w_shape))
write.csv(df_gof, paste0(out_path, "Recalibration/waitlist/df_gof.csv"), row.names = FALSE)

### Load data
set.seed(462387)
n <- 1000
df_gof <- read.csv(paste0(out_path, "Recalibration/waitlist/df_gof.csv"))

m_ddtx_coef   <- as.matrix(read.csv(paste0(out_path, "Recalibration/waitlist/ddtx_coef.csv"))[,-1])
m_ldtx_coef   <- as.matrix(read.csv(paste0(out_path, "Recalibration/waitlist/ldtx_coef.csv"))[,-1])
m_mort_coef   <- as.matrix(read.csv(paste0(out_path, "Recalibration/waitlist/mort_coef.csv"))[,-1])
m_remove_coef <- as.matrix(read.csv(paste0(out_path, "Recalibration/waitlist/remove_coef.csv"))[,-1])

df_random_shape <- df_gof[sample(nrow(df_gof), size = n, replace = TRUE, prob = df_gof$w_shape),]
v_id_shape <- df_random_shape$id

df_random_surv <- df_gof[sample(nrow(df_gof), size = n, replace = TRUE, prob = df_gof$w_surv),]
v_id_surv <- df_random_surv$id

j <- 0
m_ddtx_coef_rand_shape <- matrix(nrow = length(v_id_shape), ncol = 26)
m_ldtx_coef_rand_shape <- matrix(nrow = length(v_id_shape), ncol = 26)
m_mort_coef_rand_shape <- matrix(nrow = length(v_id_shape), ncol = 26)
m_remove_coef_rand_shape <- matrix(nrow = length(v_id_shape), ncol = 26)

m_ddtx_shape_rand_shape <- matrix(nrow = length(v_id_shape), ncol = 1)
m_ldtx_shape_rand_shape <- matrix(nrow = length(v_id_shape), ncol = 1)
m_mort_shape_rand_shape <- matrix(nrow = length(v_id_shape), ncol = 1)
m_remove_shape_rand_shape <- matrix(nrow = length(v_id_shape), ncol = 1)

for(i in v_id_shape){
  j <- j + 1
  m_ddtx_coef_rand_shape[j,] <- m_ddtx_coef[i, -c(1, 28)]
  m_ldtx_coef_rand_shape[j,] <- m_ldtx_coef[i, -c(1, 28)]
  m_mort_coef_rand_shape[j,] <- m_mort_coef[i, -c(1, 28)]
  m_remove_coef_rand_shape[j,] <- m_remove_coef[i, -c(1, 28)]

  m_ddtx_shape_rand_shape[j] <- m_ddtx_coef[i, 28]
  m_ldtx_shape_rand_shape[j] <- m_ldtx_coef[i, 28]
  m_mort_shape_rand_shape[j] <- m_mort_coef[i, 28]
  m_remove_shape_rand_shape[j] <- m_remove_coef[i, 28]
}

j <- 0
m_ddtx_coef_rand_surv <- matrix(nrow = length(v_id_surv), ncol = 26)
m_ldtx_coef_rand_surv <- matrix(nrow = length(v_id_surv), ncol = 26)
m_mort_coef_rand_surv <- matrix(nrow = length(v_id_surv), ncol = 26)
m_remove_coef_rand_surv <- matrix(nrow = length(v_id_surv), ncol = 26)

m_ddtx_shape_rand_surv <- matrix(nrow = length(v_id_surv), ncol = 1)
m_ldtx_shape_rand_surv <- matrix(nrow = length(v_id_surv), ncol = 1)
m_mort_shape_rand_surv <- matrix(nrow = length(v_id_surv), ncol = 1)
m_remove_shape_rand_surv <- matrix(nrow = length(v_id_surv), ncol = 1)

for(i in v_id_surv){
  j <- j + 1
  m_ddtx_coef_rand_surv[j,] <- m_ddtx_coef[i, -c(1, 28)]
  m_ldtx_coef_rand_surv[j,] <- m_ldtx_coef[i, -c(1, 28)]
  m_mort_coef_rand_surv[j,] <- m_mort_coef[i, -c(1, 28)]
  m_remove_coef_rand_surv[j,] <- m_remove_coef[i, -c(1, 28)]

  m_ddtx_shape_rand_surv[j] <- m_ddtx_coef[i, 28]
  m_ldtx_shape_rand_surv[j] <- m_ldtx_coef[i, 28]
  m_mort_shape_rand_surv[j] <- m_mort_coef[i, 28]
  m_remove_shape_rand_surv[j] <- m_remove_coef[i, 28]
}

## Set up parallel sessions
comb <- function(...) {
  mapply('rbind', ..., SIMPLIFY=FALSE)
}
parallel::detectCores()
n.cores <- 50

my.cluster <- parallel::makeCluster(
  n.cores,
  type = "PSOCK"
)

registerDoSNOW(my.cluster)

seed_num <- sample(1:999999, 1)
# doParallel::registerDoParallel(cl = my.cluster)
doRNG::registerDoRNG(seed = seed_num)

pb <- txtProgressBar(min=1, max=n, style=3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress=progress)

out2 <- foreach(i = 1:n, .combine='comb', .multicombine=TRUE, .options.snow=opts,
               .packages = c("MASS", "copula", "psych", "dplyr", "survival",
                             "flexsurv")) %dopar% {
                               ## Run Waitlist Simulation (survival)
                               df_surv_rand <- list.simulation(n=100000,
                                                              seed = round(runif(1, min=1, max=99999)),
                                                              data = sample_df,
                                                              ddtx_coef = m_ddtx_coef_rand_surv[i,],
                                                              ddtx_p = m_ddtx_shape_rand_surv[i,1],
                                                              ldtx_coef = m_ldtx_coef_rand_surv[i,],
                                                              ldtx_gamma = m_ldtx_shape_rand_surv[i,1],
                                                              mort_coef = m_mort_coef_rand_surv[i,],
                                                              mort_shape = m_mort_shape_rand_surv[i,1],
                                                              remove_coef = m_remove_coef_rand_surv[i,],
                                                              remove_p = m_remove_shape_rand_surv[i,1])

                               ## Clean data
                               df_surv_rand <- df_surv_rand %>%
                                 mutate(rec_age_at_tx_c = can_age_at_listing_c + (months_to_event/12)) %>%
                                 mutate(years_dial = baseline_yrs_dial + (months_to_event/12)) %>%
                                 mutate(tx_year_c = list_year_c + (months_to_event/12)) %>%
                                 mutate(tx_year_c = ifelse(tx_year_c > 10, 10, tx_year_c))  %>%
                                 mutate(start_time = ifelse(can_age_at_listing_c >=65, 0,
                                                            abs(can_age_at_listing_c)-0.5))

                               ## Run Waitlist Simulation (survival)
                               df_shape_rand <- list.simulation(n=100000,
                                                                seed = round(runif(1, min=1, max=99999)),
                                                                data = sample_df,
                                                                ddtx_coef = m_ddtx_coef_rand_shape[i,],
                                                                ddtx_p = m_ddtx_shape_rand_shape[i,1],
                                                                ldtx_coef = m_ldtx_coef_rand_shape[i,],
                                                                ldtx_gamma = m_ldtx_shape_rand_shape[i,1],
                                                                mort_coef = m_mort_coef_rand_shape[i,],
                                                                mort_shape = m_mort_shape_rand_shape[i,1],
                                                                remove_coef = m_remove_coef_rand_shape[i,],
                                                                remove_p = m_remove_shape_rand_shape[i,1])

                               ## Clean data
                               df_shape_rand <- df_shape_rand %>%
                                 mutate(rec_age_at_tx_c = can_age_at_listing_c + (months_to_event/12)) %>%
                                 mutate(years_dial = baseline_yrs_dial + (months_to_event/12)) %>%
                                 mutate(tx_year_c = list_year_c + (months_to_event/12)) %>%
                                 mutate(tx_year_c = ifelse(tx_year_c > 10, 10, tx_year_c))  %>%
                                 mutate(start_time = ifelse(can_age_at_listing_c >=65, 0,
                                                            abs(can_age_at_listing_c)-0.5))

                               df_surv_rand$run <- i
                               df_shape_rand$run <- i

                               ## Combine
                               list(df_surv_rand, df_shape_rand)

                             }
parallel::stopCluster(cl = my.cluster)

df_surv_rand <- out2[[1]]
df_shape_rand <- out2[[2]]
remove(out2)
fwrite(df_surv_rand,
       file = paste0(out_path, "Recalibration/waitlist/df_surv_rand.csv"),
       row.names = FALSE)
fwrite(df_shape_rand,
       file = paste0(out_path, "Recalibration/waitlist/df_shape_rand.csv"),
       row.names = FALSE)

##### Start code here #####
df_surv_rand <- fread(file = paste0(out_path, "Recalibration/waitlist/df_surv_rand.csv"),
                      data.table = FALSE)
df_shape_rand <- fread(file = paste0(out_path, "Recalibration/waitlist/df_shape_rand.csv"),
                       data.table = FALSE)

obs_df  <- read.csv(file = paste0(in_path, "simulation_list3.csv"),
                    sep=",")
obs_df <- obs_df %>%
  mutate(ddtx = ifelse(event_cd==0, 1, 0)) %>%
  mutate(ldtx = ifelse(event_cd==1, 1, 0)) %>%
  mutate(mort = ifelse(event_cd==2, 1, 0)) %>%
  mutate(remove = ifelse(event_cd==3, 1, 0)) %>%
  dplyr::select(-c(event_cd, no_miss)) %>%
  mutate(rec_age_at_tx_c = can_age_at_listing_c + (months_to_event/12)) %>%
  mutate(years_dial = baseline_yrs_dial + (months_to_event/12)) %>%
  mutate(tx_year_c = list_year_c + (months_to_event/12)) %>%
  mutate(race = ifelse(black==0 & hispanic==0 & other_race==0, 0,
                       ifelse(black==1 & hispanic==0 & other_race==0, 1,
                              ifelse(black==0 & hispanic==1 & other_race==0, 2, 3)))) %>%
  mutate(race = factor(race, levels = c(0, 1, 2, 3),
                       labels = c("White", "Black",
                                  "Hispanic", "Other"))) %>%
  mutate(constant = rep(0, nrow(obs_df))) %>%
  dplyr::select("can_age_at_listing_c", "male", "black", "hispanic",
                "other_race", "baseline_yrs_dial", "list_year_c", "canhx_cpra",
                "can_blood_ab", "can_blood_b", "can_blood_o", "diab_stat",
                "copd", "pvd", "ang_cad", "optn_reg2", "optn_reg3", "optn_reg4",
                "optn_reg5", "optn_reg6", "optn_reg7", "optn_reg8", "optn_reg9",
                "optn_reg10", "optn_reg11", "constant", "event", "ddtx", "ldtx",
                "mort", "remove", "months_to_event", "race", "rec_age_at_tx_c",
                "years_dial", "tx_year_c", "start_time")

obs_df$ddtx[is.na(obs_df$ddtx)] <- 0
obs_df$ldtx[is.na(obs_df$ldtx)] <- 0
obs_df$mort[is.na(obs_df$mort)] <- 0
obs_df$remove[is.na(obs_df$remove)] <- 0


# all_df_surv <- rbind(obs_df, df_surv_rand)
# all_df_surv <- all_df_surv %>%
#   mutate(diab_stat = factor(diab_stat, levels = c(0, 1),
#                             labels = c("No", "Yes")))
#
# all_df_shape <- rbind(obs_df, df_shape_rand)
# all_df_shape <- all_df_shape %>%
#   mutate(diab_stat = factor(diab_stat, levels = c(0, 1),
#                             labels = c("No", "Yes")))

### Start KM curve code
library(survival)
library(survminer)
library(gridExtra)

# Sample data frames
v_keep <- sample(seq(from = 1, to = 1000, by = 1), size = 100, replace = FALSE)
v_keep <- sort(v_keep)
obs_df$run <- (max(v_keep)+1)

df_surv_rand_keep  <- df_surv_rand[df_surv_rand$run %in% v_keep, ]
df_shape_rand_keep <- df_shape_rand[df_shape_rand$run %in% v_keep,]

# Combine with observed data
all_df_surv <- rbind(obs_df, df_surv_rand_keep) %>%
  dplyr::select(c("start_time", "months_to_event", "ddtx", "ldtx", "mort",
                  "remove", "run", "diab_stat", "race")) %>%
  mutate(diab_stat = factor(diab_stat, levels = c(0, 1),
                            labels = c("No", "Yes")))
all_df_shape <- rbind(obs_df, df_shape_rand_keep) %>%
  dplyr::select(c("start_time", "months_to_event", "ddtx", "ldtx", "mort",
                  "remove", "run", "diab_stat", "race")) %>%
  mutate(diab_stat = factor(diab_stat, levels = c(0, 1),
                            labels = c("No", "Yes")))

top_group <- c(max(v_keep), max(v_keep)+1)
col_vector <- c(rep("gray", 100), "red")
names(col_vector) <- c(v_keep, max(v_keep)+1)

# score_surv    <- df_gof$w_surv
# score_shape   <- df_gof$w_shape

# score_surv  <- df_random_surv[v_keep,]$w_surv
# score_shape <- df_random_shape[v_keep,]$w_shape
#
# mean_score_surv  <- mean(score_surv)*2
# mean_score_shape <- mean(score_shape)*2
#
# v_weight_surv  <- c(score_surv, mean_score_surv)*40
# v_weight_shape <- c(score_shape, mean_score_shape)*40

g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)}

km_plot <- function(km_fit, outcome, v_weight = NULL, subgroup = "", n){
  require(survminer)

  if(subgroup=="") {
    tmp_plot <- ggsurvplot(km_fit,
                           censor = FALSE,
                           conf.int = FALSE,
                           title = outcome,
                           short.panel.labs = TRUE,
                           legend.title = "",
                           ylim = c(0, 1),
                           xlab = "Months") +
      theme(plot.title = element_text(size = 12)) +
      theme(legend.text = element_text(size = 6)) +
      theme(axis.text.x = element_text(size = 8)) +
      theme(axis.text.y = element_text(size = 8)) +
      theme(axis.title.x = element_text(size = 8)) +
      theme(axis.title.y = element_text(size = 8)) +
      guides(fill = guide_legend(nrow = 1)) +
      theme(legend.position='none')+
      theme(strip.text.x = element_text(size = 8)) +
      scale_color_manual(values = c(rep("gray", n), "red")) +
      guides(color = guide_none())
      # scale_size_manual(values = v_weight)
  } else {
    tmp_plot <- ggsurvplot(km_fit,
                           censor = FALSE,
                           conf.int = FALSE,
                           facet.by = subgroup,
                           title = outcome,
                           short.panel.labs = TRUE,
                           legend.title = "",
                           ylim = c(0, 1),
                           xlab = "Months",
                           palette = c(rep("gray", n), "red"),) +
      theme(plot.title = element_text(size = 12)) +
      theme(legend.text = element_text(size = 6)) +
      theme(axis.text.x = element_text(size = 8)) +
      theme(axis.text.y = element_text(size = 8)) +
      theme(axis.title.x = element_text(size = 8)) +
      theme(axis.title.y = element_text(size = 8)) +
      guides(fill = guide_legend(nrow = 1)) +
      theme(legend.position='none')+
      theme(strip.text.x = element_text(size = 8))
      # scale_size_manual(values = v_weight)
  }

  return(tmp_plot)
}

## KM fit, survival-based GOF
fit_ddtx_surv    <- survfit(Surv(time = start_time, time2 = months_to_event, ddtx) ~ run,
                       data = all_df_surv)  # Deceased Donor Tx

fit_ldtx_surv    <- survfit(Surv(time = start_time, time2 = months_to_event, ldtx) ~ run,
                       data = all_df_surv)  # Living Donor Tx

fit_mort_surv    <- survfit(Surv(time = start_time, time2 = months_to_event, mort) ~ run,
                       data = all_df_surv)  # Waitlist Mortality

fit_remove_surv  <- survfit(Surv(time = start_time, time2 = months_to_event, remove) ~ run,
                       data = all_df_surv)  # Other Removal

## KM fit, shape-based GOF
fit_ddtx_shape    <- survfit(Surv(time = start_time, time2 = months_to_event, ddtx) ~ run,
                            data = all_df_shape)  # Deceased Donor Tx

fit_ldtx_shape    <- survfit(Surv(time = start_time, time2 = months_to_event, ldtx) ~ run,
                            data = all_df_shape)  # Living Donor Tx

fit_mort_shape    <- survfit(Surv(time = start_time, time2 = months_to_event, mort) ~ run,
                            data = all_df_shape)  # Waitlist Mortality

fit_remove_shape  <- survfit(Surv(time = start_time, time2 = months_to_event, remove) ~ run,
                            data = all_df_shape)  # Other Removal


## KM curves, survival-based
ddtx_all_surv <- ggsurvplot(fit_ddtx_surv,
                       censor = FALSE,
                       conf.int = FALSE,
                       title = "Deceased Donor Tx",
                       short.panel.labs = TRUE,
                       legend.title = "",
                       ylim = c(0, 1),
                       xlab = "Months")
ddtx_all_surv$plot <- ddtx_all_surv$plot +
  theme(plot.title = element_text(size = 12)) +
  theme(legend.text = element_text(size = 6)) +
  theme(axis.text.x = element_text(size = 8)) +
  theme(axis.text.y = element_text(size = 8)) +
  theme(axis.title.x = element_text(size = 8)) +
  theme(axis.title.y = element_text(size = 8)) +
  guides(fill = guide_legend(nrow = 1)) +
  theme(legend.position='none')+
  theme(strip.text.x = element_text(size = 8)) +
  scale_color_manual(values = c(rep("gray", 100), "red")) +
  guides(color = guide_none()) +
  scale_size_manual(values = v_weight_surv)

ldtx_all_surv <- ggsurvplot(fit_ldtx_surv,
                       censor = FALSE,
                       conf.int = FALSE,
                       title = "Living Donor Tx",
                       short.panel.labs = TRUE,
                       legend.title = "",
                       ylim = c(0, 1),
                       xlab = "Months")
ldtx_all_surv$plot <- ldtx_all_surv$plot +
  theme(plot.title = element_text(size = 12)) +
  theme(legend.text = element_text(size = 6)) +
  theme(axis.text.x = element_text(size = 8)) +
  theme(axis.text.y = element_text(size = 8)) +
  theme(axis.title.x = element_text(size = 8)) +
  theme(axis.title.y = element_text(size = 8)) +
  guides(fill = guide_legend(nrow = 1)) +
  theme(legend.position='none')+
  theme(strip.text.x = element_text(size = 8)) +
  scale_color_manual(values = c(rep("gray", 100), "red")) +
  guides(color = guide_none())+
  scale_size_manual(values = v_weight_surv)

mort_all_surv <- ggsurvplot(fit_mort_surv,
                       censor = FALSE,
                       conf.int = FALSE,
                       title = "Waitlist Mortality",
                       short.panel.labs = TRUE,
                       legend.title = "",
                       ylim = c(0, 1),
                       xlab = "Months")
mort_all_surv$plot <- mort_all_surv$plot +
  theme(plot.title = element_text(size = 12)) +
  theme(legend.text = element_text(size = 6)) +
  theme(axis.text.x = element_text(size = 8)) +
  theme(axis.text.y = element_text(size = 8)) +
  theme(axis.title.x = element_text(size = 8)) +
  theme(axis.title.y = element_text(size = 8)) +
  guides(fill = guide_legend(nrow = 1)) +
  theme(legend.position='none')+
  theme(strip.text.x = element_text(size = 8)) +
  scale_color_manual(values = c(rep("gray", 100), "red")) +
  guides(color = guide_none())+
  scale_size_manual(values = v_weight_surv)

remove_all_surv <- ggsurvplot(fit_remove_surv,
                         censor = FALSE,
                         conf.int = FALSE,
                         title = "Other Waitlist Removal",
                         short.panel.labs = TRUE,
                         legend.title = "",
                         ylim = c(0, 1),
                         xlab = "Months")
remove_all_surv$plot <- remove_all_surv$plot +
  theme(plot.title = element_text(size = 12)) +
  theme(legend.text = element_text(size = 6)) +
  theme(axis.text.x = element_text(size = 8)) +
  theme(axis.text.y = element_text(size = 8)) +
  theme(axis.title.x = element_text(size = 8)) +
  theme(axis.title.y = element_text(size = 8)) +
  theme(legend.position='none') +
  theme(strip.text.x = element_text(size = 8)) +
  scale_color_manual(values = c(rep("gray", 100), "red")) +
  guides(color = guide_none())+
  scale_size_manual(values = v_weight_surv)


remove_legend <- ggsurvplot(fit_remove_surv,
                            censor = FALSE,
                            conf.int = FALSE,
                            title = "Other Waitlist Removal",
                            short.panel.labs = TRUE,
                            legend.title = "",
                            ylim = c(0, 1),
                            xlab = "Months")
remove_legend$plot <- remove_legend$plot +
  theme(plot.title = element_text(size = 12)) +
  theme(legend.text = element_text(size = 10)) +
  theme(axis.text.x = element_text(size = 8)) +
  theme(axis.text.y = element_text(size = 8)) +
  theme(axis.title.x = element_text(size = 8)) +
  theme(axis.title.y = element_text(size = 8)) +
  theme(legend.position='right') +
  theme(strip.text.x = element_text(size = 8)) +
  scale_color_manual(values = col_vector, breaks = top_group,
                     labels = c("Simulated", "Observed")) +
  theme(legend.position=c(.5,.5)) +
  guides(color = guide_legend("Type")) +
  scale_size_manual(values = v_weight_surv)


plots_all_surv <- grid.arrange(arrangeGrob(ddtx_all_surv$plot, ldtx_all_surv$plot,
                                      mort_all_surv$plot, remove_all_surv$plot +
                                        theme(legend.position="none"),
                                      nrow = 2),
                          nrow=2,heights=c(10, 1))
ggsave(filename = paste0(out_path, "waitlist_surv.png"),
       plot = plots_all_surv , width = 7, height = 4)

## KM curves, shape-based
ddtx_all_shape <- ggsurvplot(fit_ddtx_shape,
                            censor = FALSE,
                            conf.int = FALSE,
                            title = "Deceased Donor Tx",
                            short.panel.labs = TRUE,
                            legend.title = "",
                            ylim = c(0, 1),
                            xlab = "Months")
ddtx_all_shape$plot <- ddtx_all_shape$plot +
  theme(plot.title = element_text(size = 12)) +
  theme(legend.text = element_text(size = 6)) +
  theme(axis.text.x = element_text(size = 8)) +
  theme(axis.text.y = element_text(size = 8)) +
  theme(axis.title.x = element_text(size = 8)) +
  theme(axis.title.y = element_text(size = 8)) +
  guides(fill = guide_legend(nrow = 1)) +
  theme(legend.position='none')+
  theme(strip.text.x = element_text(size = 8)) +
  scale_color_manual(values = c(rep("gray", 100), "red")) +
  guides(color = guide_none()) +
  scale_size_manual(values = v_weight_shape)

ldtx_all_shape <- ggsurvplot(fit_ldtx_shape,
                            censor = FALSE,
                            conf.int = FALSE,
                            title = "Living Donor Tx",
                            short.panel.labs = TRUE,
                            legend.title = "",
                            ylim = c(0, 1),
                            xlab = "Months")
ldtx_all_shape$plot <- ldtx_all_shape$plot +
  theme(plot.title = element_text(size = 12)) +
  theme(legend.text = element_text(size = 6)) +
  theme(axis.text.x = element_text(size = 8)) +
  theme(axis.text.y = element_text(size = 8)) +
  theme(axis.title.x = element_text(size = 8)) +
  theme(axis.title.y = element_text(size = 8)) +
  guides(fill = guide_legend(nrow = 1)) +
  theme(legend.position='none')+
  theme(strip.text.x = element_text(size = 8)) +
  scale_color_manual(values = c(rep("gray", 100), "red")) +
  guides(color = guide_none())+
  scale_size_manual(values = v_weight_shape)

mort_all_shape <- ggsurvplot(fit_mort_shape,
                            censor = FALSE,
                            conf.int = FALSE,
                            title = "Waitlist Mortality",
                            short.panel.labs = TRUE,
                            legend.title = "",
                            ylim = c(0, 1),
                            xlab = "Months")
mort_all_shape$plot <- mort_all_shape$plot +
  theme(plot.title = element_text(size = 12)) +
  theme(legend.text = element_text(size = 6)) +
  theme(axis.text.x = element_text(size = 8)) +
  theme(axis.text.y = element_text(size = 8)) +
  theme(axis.title.x = element_text(size = 8)) +
  theme(axis.title.y = element_text(size = 8)) +
  guides(fill = guide_legend(nrow = 1)) +
  theme(legend.position='none')+
  theme(strip.text.x = element_text(size = 8)) +
  scale_color_manual(values = c(rep("gray", 100), "red")) +
  guides(color = guide_none())+
  scale_size_manual(values = v_weight_shape)

remove_all_shape <- ggsurvplot(fit_remove_shape,
                              censor = FALSE,
                              conf.int = FALSE,
                              title = "Other Waitlist Removal",
                              short.panel.labs = TRUE,
                              legend.title = "",
                              ylim = c(0, 1),
                              xlab = "Months")
remove_all_shape$plot <- remove_all_shape$plot +
  theme(plot.title = element_text(size = 12)) +
  theme(legend.text = element_text(size = 6)) +
  theme(axis.text.x = element_text(size = 8)) +
  theme(axis.text.y = element_text(size = 8)) +
  theme(axis.title.x = element_text(size = 8)) +
  theme(axis.title.y = element_text(size = 8)) +
  theme(legend.position='none') +
  theme(strip.text.x = element_text(size = 8)) +
  scale_color_manual(values = c(rep("gray", 100), "red")) +
  guides(color = guide_none())+
  scale_size_manual(values = v_weight_shape)


remove_legend <- ggsurvplot(fit_remove_shape,
                            censor = FALSE,
                            conf.int = FALSE,
                            title = "Other Waitlist Removal",
                            short.panel.labs = TRUE,
                            legend.title = "",
                            ylim = c(0, 1),
                            xlab = "Months")
remove_legend$plot <- remove_legend$plot +
  theme(plot.title = element_text(size = 12)) +
  theme(legend.text = element_text(size = 10)) +
  theme(axis.text.x = element_text(size = 8)) +
  theme(axis.text.y = element_text(size = 8)) +
  theme(axis.title.x = element_text(size = 8)) +
  theme(axis.title.y = element_text(size = 8)) +
  theme(legend.position='right') +
  theme(strip.text.x = element_text(size = 8)) +
  scale_color_manual(values = col_vector, breaks = top_group,
                     labels = c("Simulated", "Observed")) +
  theme(legend.position=c(.5,.5)) +
  guides(color = guide_legend("Type")) +
  scale_size_manual(values = v_weight_shape)


plots_all_shape <- grid.arrange(arrangeGrob(ddtx_all_shape$plot, ldtx_all_shape$plot,
                                           mort_all_shape$plot, remove_all_shape$plot +
                                             theme(legend.position="none"),
                                           nrow = 2),
                               nrow=2,heights=c(10, 1))
ggsave(filename = paste0(out_path, "waitlist_shape.png"),
       plot = plots_all_shape , width = 7, height = 4)

### By subgroup
## KM fits
# Race, survival-based GOF
fit_ddtx_surv_race    <- survfit(Surv(time = start_time, time2 = months_to_event, ddtx) ~ run + race,
                            data = all_df_surv)  # Deceased Donor Tx

fit_ldtx_surv_race    <- survfit(Surv(time = start_time, time2 = months_to_event, ldtx) ~ run + race,
                            data = all_df_surv)  # Living Donor Tx

fit_mort_surv_race    <- survfit(Surv(time = start_time, time2 = months_to_event, mort) ~ run + race,
                            data = all_df_surv)  # Waitlist Mortality

fit_remove_surv_race  <- survfit(Surv(time = start_time, time2 = months_to_event, remove) ~ run + race,
                            data = all_df_surv)  # Other Removal

# Race, shape-based GOF
fit_ddtx_shape_race    <- survfit(Surv(time = start_time, time2 = months_to_event, ddtx) ~ run + race,
                             data = all_df_shape)  # Deceased Donor Tx

fit_ldtx_shape_race    <- survfit(Surv(time = start_time, time2 = months_to_event, ldtx) ~ run + race,
                             data = all_df_shape)  # Living Donor Tx

fit_mort_shape_race    <- survfit(Surv(time = start_time, time2 = months_to_event, mort) ~ run + race,
                             data = all_df_shape)  # Waitlist Mortality

fit_remove_shape_race  <- survfit(Surv(time = start_time, time2 = months_to_event, remove) ~ run + race,
                             data = all_df_shape)  # Other Removal

# Diabetes, survival-based GOF
fit_ddtx_surv_diab    <- survfit(Surv(time = start_time, time2 = months_to_event, ddtx) ~ run + diab_stat,
                                 data = all_df_surv)  # Deceased Donor Tx

fit_ldtx_surv_diab    <- survfit(Surv(time = start_time, time2 = months_to_event, ldtx) ~ run + diab_stat,
                                 data = all_df_surv)  # Living Donor Tx

fit_mort_surv_diab    <- survfit(Surv(time = start_time, time2 = months_to_event, mort) ~ run + diab_stat,
                                 data = all_df_surv)  # Waitlist Mortality

fit_remove_surv_diab  <- survfit(Surv(time = start_time, time2 = months_to_event, remove) ~ run + diab_stat,
                                 data = all_df_surv)  # Other Removal

# Diabetes, shape-based GOF
fit_ddtx_shape_diab    <- survfit(Surv(time = start_time, time2 = months_to_event, ddtx) ~ run + diab_stat,
                                  data = all_df_shape)  # Deceased Donor Tx

fit_ldtx_shape_diab    <- survfit(Surv(time = start_time, time2 = months_to_event, ldtx) ~ run + diab_stat,
                                  data = all_df_shape)  # Living Donor Tx

fit_mort_shape_diab    <- survfit(Surv(time = start_time, time2 = months_to_event, mort) ~ run + diab_stat,
                                  data = all_df_shape)  # Waitlist Mortality

fit_remove_shape_diab  <- survfit(Surv(time = start_time, time2 = months_to_event, remove) ~ run + diab_stat,
                                  data = all_df_shape)  # Other Removal

# Make plots, survival-based GOF
ddtx_surv_race <- km_plot(km_fit = fit_ddtx_surv_race, outcome = "Deceased Donor Tx",
                     subgroup = "race", n = 100)
ldtx_surv_race <- km_plot(km_fit = fit_ldtx_surv_race, outcome = "Living Donor Tx",
                     subgroup = "race", n = 100)
mort_surv_race <- km_plot(km_fit = fit_mort_surv_race, outcome = "Waitlist Mortality",
                     subgroup = "race", n = 100)
remove_surv_race <- km_plot(km_fit = fit_remove_surv_race,
                       outcome = "Other Waitlist Removal",
                       subgroup = "race", n = 100)

plots_surv_race <- grid.arrange(arrangeGrob(ddtx_surv_race, ldtx_surv_race,
                                            mort_surv_race, remove_surv_race +
                                              theme(legend.position="none"),
                                            nrow = 2),
                                nrow=2,heights=c(10, 1))
ggsave(filename = paste0(out_path, "waitlist_surv_race.png"),
       plot = plots_surv_race , width = 7, height = 7)

ddtx_surv_diab <- km_plot(km_fit = fit_ddtx_surv_diab, outcome = "Deceased Donor Tx",
                          subgroup = "diab_stat", n = 100)
ldtx_surv_diab <- km_plot(km_fit = fit_ldtx_surv_diab, outcome = "Living Donor Tx",
                          subgroup = "diab_stat", n = 100)
mort_surv_diab <- km_plot(km_fit = fit_mort_surv_diab, outcome = "Waitlist Mortality",
                          subgroup = "diab_stat", n = 100)
remove_surv_diab <- km_plot(km_fit = fit_remove_surv_diab,
                            outcome = "Other Waitlist Removal",
                            subgroup = "diab_stat", n = 100)
plots_surv_diab <- grid.arrange(arrangeGrob(ddtx_surv_diab, ldtx_surv_diab,
                                            mort_surv_diab, remove_surv_diab +
                                              theme(legend.position="none"),
                                            nrow = 2),
                                nrow=2,heights=c(10, 1))
ggsave(filename = paste0(out_path, "waitlist_surv_diab.png"),
       plot = plots_surv_diab , width = 7, height = 7)
# Make plots, shape-based GOF
ddtx_shape_race <- km_plot(km_fit = fit_ddtx_shape_race, outcome = "Deceased Donor Tx",
                     subgroup = "race", n = 100)
ldtx_shape_race <- km_plot(km_fit = fit_ldtx_shape_race, outcome = "Living Donor Tx",
                     subgroup = "race", n = 100)
mort_shape_race <- km_plot(km_fit = fit_mort_shape_race, outcome = "Waitlist Mortality",
                     subgroup = "race", n = 100)
remove_shape_race <- km_plot(km_fit = fit_remove_shape_race,
                       outcome = "Other Waitlist Removal",
                       subgroup = "race", n = 100)
plots_shape_race <- grid.arrange(arrangeGrob(ddtx_shape_race, ldtx_shape_race,
                                            mort_shape_race, remove_shape_race +
                                              theme(legend.position="none"),
                                            nrow = 2),
                                nrow=2,heights=c(10, 1))
ggsave(filename = paste0(out_path, "waitlist_shape_race.png"),
       plot = plots_shape_race , width = 7, height = 7)

ddtx_shape_diab <- km_plot(km_fit = fit_ddtx_shape_diab, outcome = "Deceased Donor Tx",
                           subgroup = "diab_stat", n = 100)
ldtx_shape_diab <- km_plot(km_fit = fit_ldtx_shape_diab, outcome = "Living Donor Tx",
                           subgroup = "diab_stat", n = 100)
mort_shape_diab <- km_plot(km_fit = fit_mort_shape_diab, outcome = "Waitlist Mortality",
                           subgroup = "diab_stat", n = 100)
remove_shape_diab <- km_plot(km_fit = fit_remove_shape_diab,
                             outcome = "Other Waitlist Removal",
                             subgroup = "diab_stat", n = 100)
plots_shape_diab <- grid.arrange(arrangeGrob(ddtx_shape_diab, ldtx_shape_diab,
                                            mort_shape_diab, remove_shape_diab +
                                              theme(legend.position="none"),
                                            nrow = 2),
                                nrow=2,heights=c(10, 1))
ggsave(filename = paste0(out_path, "waitlist_shape_diab.png"),
       plot = plots_shape_diab , width = 7, height = 7)

### Table 1
library(modi)
library(tidyverse)
library(xtable)

table1_clean <- function(data){
  data <- data %>%
    mutate(can_blood_a = ifelse(can_blood_ab==0 & can_blood_b==0 & can_blood_o==0,
                                1, 0)) %>%
    mutate(rec_age_at_tx = rec_age_at_tx_c + 65) %>%
    mutate(optn_reg1 = ifelse(optn_reg2==0 & optn_reg3==0 & optn_reg4==0 &
                                optn_reg5==0 & optn_reg6==0 & optn_reg7==0 &
                                optn_reg8==0 & optn_reg9==0 & optn_reg10==0 &
                                optn_reg11==0, 1, 0)) %>%
    mutate(white = ifelse(black==0 & hispanic==0 & other_race==0, 1, 0)) %>%
    filter(ddtx==1)

  # label(data$rec_age_at_tx)   <- "Age at Transplant"
  # label(data$male)            <- "Sex"
  # label(data$race)            <- "Race/Ethnicity"
  # label(data$diab_stat)       <- "Diabetes Status"
  # label(data$copd)            <- "History of COPD"
  # label(data$pvd)             <- "History of PVD"
  # label(data$ang_cad)         <- "History of Angina/CAD"
  # label(data$dgf)             <- "History of DGF"
  # label(data$years_dial)      <- "Years on Dialysis Before Tx"
  # label(data$canhx_cpra)      <- "cPRA"
  # label(data$months_to_event) <- "Months To Tx"

  sim_table1 <- matrix(data = NA, nrow = 1000, ncol = 27)

  sim_table1[,1]  <- aggregate(data$rec_age_at_tx, list(data$run), FUN=mean)[,2]
  sim_table1[,2]  <- aggregate(data$male, list(data$run), FUN=mean)[,2]
  sim_table1[,3]  <- aggregate(data$white, list(data$run), FUN=mean)[,2]
  sim_table1[,4]  <- aggregate(data$black, list(data$run), FUN=mean)[,2]
  sim_table1[,5]  <- aggregate(data$hispanic, list(data$run), FUN=mean)[,2]
  sim_table1[,6]  <- aggregate(data$other_race, list(data$run), FUN=mean)[,2]
  sim_table1[,7]  <- aggregate(data$can_blood_a, list(data$run), FUN=mean)[,2]
  sim_table1[,8]  <- aggregate(data$can_blood_ab, list(data$run), FUN=mean)[,2]
  sim_table1[,9]  <- aggregate(data$can_blood_b, list(data$run), FUN=mean)[,2]
  sim_table1[,10] <- aggregate(data$can_blood_o, list(data$run), FUN=mean)[,2]
  sim_table1[,11] <- aggregate(data$years_dial, list(data$run), FUN=mean)[,2]
  sim_table1[,12] <- aggregate(data$diab_stat, list(data$run), FUN=mean)[,2]
  sim_table1[,13] <- aggregate(data$copd, list(data$run), FUN=mean)[,2]
  sim_table1[,14] <- aggregate(data$pvd, list(data$run), FUN=mean)[,2]
  sim_table1[,15] <- aggregate(data$ang_cad, list(data$run), FUN=mean)[,2]
  sim_table1[,16] <- aggregate(data$canhx_cpra, list(data$run), FUN=mean)[,2]
  sim_table1[,17] <- aggregate(data$optn_reg1, list(data$run), FUN=mean)[,2]
  sim_table1[,18] <- aggregate(data$optn_reg2, list(data$run), FUN=mean)[,2]
  sim_table1[,19] <- aggregate(data$optn_reg3, list(data$run), FUN=mean)[,2]
  sim_table1[,20] <- aggregate(data$optn_reg4, list(data$run), FUN=mean)[,2]
  sim_table1[,21] <- aggregate(data$optn_reg5, list(data$run), FUN=mean)[,2]
  sim_table1[,22] <- aggregate(data$optn_reg6, list(data$run), FUN=mean)[,2]
  sim_table1[,23] <- aggregate(data$optn_reg7, list(data$run), FUN=mean)[,2]
  sim_table1[,24] <- aggregate(data$optn_reg8, list(data$run), FUN=mean)[,2]
  sim_table1[,25] <- aggregate(data$optn_reg9, list(data$run), FUN=mean)[,2]
  sim_table1[,26] <- aggregate(data$optn_reg10, list(data$run), FUN=mean)[,2]
  sim_table1[,27] <- aggregate(data$optn_reg11, list(data$run), FUN=mean)[,2]

  sim_table1 <- as.data.frame(sim_table1)
  return(sim_table1)
}

weighted.table1 <- function(data=NULL, score = NULL) {

  means <- sapply(data, weighted.mean, w = score, na.rm=TRUE)
  means <- unlist(means)
  upper <- sapply(data, weighted.quantile, w = score, prob = 0.975, plot = FALSE)
  lower <- sapply(data, weighted.quantile, w = score, prob = 0.025, plot = FALSE)

  data_summary <- as.data.frame(rbind(means, upper, lower))

  return(data_summary)
}

nonweighted.table1 <- function(data=NULL) {
  means <- sapply(data, mean)
  means <- unlist(means)
  upper <- sapply(data, quantile, probs = c(0.975))
  lower <- sapply(data, quantile, probs = c(0.025))

  data_summary <- as.data.frame(rbind(means, upper, lower))

  return(data_summary)
}

# Table 1 - Survival-Based GOF
# df_surv_sim <- all_df_surv %>%
#   filter(run!=1001)

tx_df_surv <- table1_clean(data = df_surv_rand)
# table1_w_surv <- weighted.table1(tx_df_surv, score = df_random_surv$w_surv)
table1_surv <- nonweighted.table1(tx_df_surv)

# Table 1 - Shape-Based GOF
# df_shape_sim <- all_df_shape %>%
#   filter(run!=1001)

tx_df_shape <- table1_clean(data = df_shape_rand)
# table1_w_shape <- weighted.table1(tx_df_shape, score = df_random_shape$w_shape)
table1_shape <- nonweighted.table1(tx_df_shape)

# Observed Table 1
tx_df_obs <- table1_clean(data = obs_df)[1,]

# Difference Table 1
table1_all <- rbind(tx_df_obs, table1_surv, table1_shape)
row.names(table1_all) <- c("Observed",
                           "Survival-Based Mean", "Survival-Based 97.5", "Survival-Based 2.5",
                           "Shape-Based Mean", "Shape-Based 97.5", "Shape-Based 2.5")
colnames(table1_all) <- c("Age at Transplant", "Sex",
                          "White", "Black", "Hispanic", "Other Race/Ethnicity",
                          "Blood A", "Blood AB", "Blood B", "Blood O",
                          "Years on Dialysis Before Tx",
                          "History of Diabetes",
                          "History of COPD",
                          "History of PVD",
                          "History of Angina/CAD",
                          "cPRA",
                          "OPTN 1", "OPTN 2", "OPTN 3", "OPTN 4", "OPTN 5",
                          "OPTN 6", "OPTN 7", "OPTN 8", "OPTN 9", "OPTN 10",
                          "OPTN 11")
write.csv(table1_all, paste0(out_path, "Recalibration/waitlist/table1.csv"))


