## Load packages
library(doParallel)
library(doRNG)
library(foreach)
library(dtplyr)
library(dplyr)
library(doSNOW)
library(data.table)

## File paths
in_path  <- "//tsclient/Documents/Simulation Model/Model Inputs/"
out_path <- "//tsclient/Documents/Simulation Model/Model Outputs/Recalibration/"

## Load Waitlist Calibration Results
df_gof <- fread(paste0(out_path, "waitlist/df_gof.csv"))
m_ddtx_coef   <- as.matrix(fread(paste0(out_path, "waitlist/ddtx_coef.csv"))[,-1])
m_ldtx_coef   <- as.matrix(fread(paste0(out_path, "waitlist/ldtx_coef.csv"))[,-1])
m_mort_coef   <- as.matrix(fread(paste0(out_path, "waitlist/mort_coef.csv"))[,-1])
m_remove_coef <- as.matrix(fread(paste0(out_path, "waitlist/remove_coef.csv"))[,-1])


## Sample from waitlist calibration results
set.seed(126748)
n <- 50000

df_wait_sample <- df_gof[sample(nrow(df_gof), size = n, replace = TRUE, prob = df_gof$w_shape),]
v_wait_id <- df_wait_sample$id

j <- 0
m_ddtx_coef_samp   <- matrix(nrow = length(v_wait_id), ncol = 26)
m_ldtx_coef_samp   <- matrix(nrow = length(v_wait_id), ncol = 26)
m_mort_coef_samp   <- matrix(nrow = length(v_wait_id), ncol = 26)
m_remove_coef_samp <- matrix(nrow = length(v_wait_id), ncol = 26)

m_ddtx_shape_samp <- matrix(nrow = length(v_wait_id), ncol = 1)
m_ldtx_shape_samp <- matrix(nrow = length(v_wait_id), ncol = 1)
m_mort_shape_samp <- matrix(nrow = length(v_wait_id), ncol = 1)
m_remove_shape_samp <- matrix(nrow = length(v_wait_id), ncol = 1)

for(i in v_wait_id){
  j <- j + 1
  m_ddtx_coef_samp[j,] <- m_ddtx_coef[i, -c(1, 28)]
  m_ldtx_coef_samp[j,] <- m_ldtx_coef[i, -c(1, 28)]
  m_mort_coef_samp[j,] <- m_mort_coef[i, -c(1, 28)]
  m_remove_coef_samp[j,] <- m_remove_coef[i, -c(1, 28)]

  m_ddtx_shape_samp[j] <- m_ddtx_coef[i, 28]
  m_ldtx_shape_samp[j] <- m_ldtx_coef[i, 28]
  m_mort_shape_samp[j] <- m_mort_coef[i, 28]
  m_remove_shape_samp[j] <- m_remove_coef[i, 28]
}

id <- seq(from = 1, to = n, by = 1)
m_ddtx_coef <- cbind(id, m_ddtx_coef_samp, m_ddtx_shape_samp)
m_ldtx_coef <- cbind(id, m_ldtx_coef_samp, m_ldtx_shape_samp)
m_mort_coef <- cbind(id, m_mort_coef_samp, m_mort_shape_samp)
m_remove_coef <- cbind(id, m_remove_coef_samp, m_remove_shape_samp)

fwrite(m_ddtx_coef, paste0(out_path, "Stage2/part2/ddtx_coef.csv"))
fwrite(m_ldtx_coef, paste0(out_path, "Stage2/part2/ldtx_coef.csv"))
fwrite(m_mort_coef, paste0(out_path, "Stage2/part2/mort_coef.csv"))
fwrite(m_remove_coef, paste0(out_path, "Stage2/part2/remove_coef.csv"))

mlogit_coef <- numeric()
gs_mort_coef <- numeric()
gl_coef <- numeric()
dial_mort_coef <- numeric()
gl_mort_coef <- numeric()
gs_mort_shape <- numeric()
gl_shape <- numeric()
dial_mort_shape <- numeric()

pb1 <- txtProgressBar(min = 0, max = n, style = 3)

for(i in 1:n) {
  param_set <- psa_posttx()

  ## Store coefficients for post-tx functions
  mlogit_coef    <- rbind(mlogit_coef, t(param_set$mlogit_coef_psa))
  gs_mort_coef   <- rbind(gs_mort_coef, t(param_set$gs_mort_coef_psa))
  gl_coef        <- rbind(gl_coef, t(param_set$gl_coef_psa))
  dial_mort_coef <- rbind(dial_mort_coef, t(param_set$dial_mort_coef_psa))
  gl_mort_coef   <- rbind(gl_mort_coef, t(param_set$gl_mort_coef_psa))

  gs_mort_shape   <- rbind(gs_mort_shape, param_set$gs_mort_shape_psa)
  gl_shape        <- rbind(gl_shape, param_set$gl_shape_psa)
  dial_mort_shape <- rbind(dial_mort_shape, param_set$dial_mort_shape_psa)

  setTxtProgressBar(pb1, i)
}

gs_mort_coef_rand   <- cbind(gs_mort_coef, gs_mort_shape)
gl_coef_rand        <- cbind(gl_coef, gl_shape)
dial_mort_coef_rand <- cbind(dial_mort_coef, dial_mort_shape)

gs_mort_shape <- gs_mort_coef_rand[,ncol(gs_mort_coef_rand)]
gl_shape <- gl_coef_rand[,ncol(gl_coef_rand)]
dial_mort_shape <- dial_mort_coef_rand[,ncol(dial_mort_coef_rand)]

fwrite(mlogit_coef, paste0(out_path, "Stage2/part2/mlogit_coef.csv"))
fwrite(gs_mort_coef_rand, paste0(out_path, "Stage2/part2/gs_mort_coef.csv"))
fwrite(gl_coef_rand, paste0(out_path, "Stage2/part2/gl_coef.csv"))
fwrite(dial_mort_coef_rand, paste0(out_path, "Stage2/part2/dial_mort_coef.csv"))
fwrite(gl_mort_coef, paste0(out_path, "Stage2/part2/gl_mort_coef.csv"))

## Quick function for parametric fit using try catch
parametric.fit <- function(data, distribution, t, outcome) {
  tryCatch(
    expr = {
      fit <- flexsurv::flexsurvreg(Surv(time = t, outcome==1) ~ 1, data=data, dist=distribution)$coefficients
      return(fit)
    },
    error = function(e) {
      return(c(NA, NA))
      print(e)
    }
  )
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

doRNG::registerDoRNG(seed = seed_num)

pb <- txtProgressBar(min=1, max=n, style=3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress=progress)

out <- foreach(i = 1:n, .combine='comb', .multicombine=TRUE, .options.snow=opts,
               .packages = c("MASS", "copula", "psych", "dtplyr",
                             "survival", "flexsurv", "dplyr")) %dopar% {

                               ## Run Waitlist Simulation
                               list_sim_df <- list.simulation(
                                 n = 100000,
                                 seed = round(runif(1, min = 1, max = 99999)),
                                 data = sample_df,
                                 ddtx_coef = m_ddtx_coef_samp[i, ],
                                 ddtx_p = m_ddtx_shape_samp[i, 1],
                                 ldtx_coef = m_ldtx_coef_samp[i, ],
                                 ldtx_gamma = m_ldtx_shape_samp[i, 1],
                                 mort_coef = m_mort_coef_samp[i, ],
                                 mort_shape = m_mort_shape_samp[i, 1],
                                 remove_coef = m_remove_coef_samp[i, ],
                                 remove_p = m_remove_shape_samp[i, 1]
                               )

                               sim_df <- run.simulation(
                                 data = list_sim_df,
                                 seed = round(runif(i, min = 1, max = 99999)),
                                 gl30_coef = mlogit_coef[i, 1:18],
                                 dgf_coef = mlogit_coef[i, 19:36],
                                 mort30_coef = mlogit_coef[i, 37:54],
                                 gl_coef = gl_coef[i, ],
                                 gl_shape = gl_shape[i],
                                 gs_mort_coef = gs_mort_coef[i, ],
                                 gs_mort_shape = gs_mort_shape[i],
                                 dial_mort_coef = dial_mort_coef[i, ],
                                 dial_mort_shape = dial_mort_shape[i],
                                 gl_mort_coef = gl_mort_coef[i, ]
                               )


                               ### Assess waitlist fit
                               ## Clean and run survival models
                               list_sim_df <- list_sim_df %>%
                                 mutate(rec_age_at_tx_c = can_age_at_listing_c + (months_to_event/12)) %>%
                                 mutate(years_dial = baseline_yrs_dial + (months_to_event/12)) %>%
                                 mutate(tx_year_c = list_year_c + (months_to_event/12)) %>%
                                 mutate(tx_year_c = ifelse(tx_year_c > 10, 10, tx_year_c))

                               ### Split data by diabetes status and race
                               df_diab <- split(list_sim_df, list_sim_df$diab_stat)
                               df_race <- split(list_sim_df, list_sim_df$race)

                               ### Deceased donor transplant
                               ddtx_km_fit          <- survfit(Surv(time = months_to_event, ddtx==1) ~ 1, data = list_sim_df)
                               ddtx_km_fit_diab0    <- survfit(Surv(time = months_to_event, ddtx==1) ~ 1, data = df_diab[[1]])
                               ddtx_km_fit_diab1    <- survfit(Surv(time = months_to_event, ddtx==1) ~ 1, data = df_diab[[2]])
                               ddtx_km_fit_white    <- survfit(Surv(time = months_to_event, ddtx==1) ~ 1, data = df_race[[1]])
                               ddtx_km_fit_black    <- survfit(Surv(time = months_to_event, ddtx==1) ~ 1, data = df_race[[2]])
                               ddtx_km_fit_hispanic <- survfit(Surv(time = months_to_event, ddtx==1) ~ 1, data = df_race[[3]])
                               ddtx_km_fit_other    <- survfit(Surv(time = months_to_event, ddtx==1) ~ 1, data = df_race[[4]])

                               ddtx_surv          <- c(summary(ddtx_km_fit, times = seq(from=3, to=180, by=3), extend=TRUE)$surv)
                               ddtx_surv_diab0    <- c(summary(ddtx_km_fit_diab0, times = seq(from=3, to=180, by=3), extend=TRUE)$surv)
                               ddtx_surv_diab1    <- c(summary(ddtx_km_fit_diab1, times = seq(from=3, to=180, by=3), extend=TRUE)$surv)
                               ddtx_surv_white    <- c(summary(ddtx_km_fit_white, times = seq(from=3, to=180, by=3), extend=TRUE)$surv)
                               ddtx_surv_black    <- c(summary(ddtx_km_fit_black, times = seq(from=3, to=180, by=3), extend=TRUE)$surv)
                               ddtx_surv_hispanic <- c(summary(ddtx_km_fit_hispanic, times = seq(from=3, to=180, by=3), extend=TRUE)$surv)
                               ddtx_surv_other    <- c(summary(ddtx_km_fit_other, times = seq(from=3, to=180, by=3), extend=TRUE)$surv)

                               # Parametric Fits
                               ddtx_shape_all      <- flexsurvreg(Surv(time = months_to_event, ddtx==1)~1, data=list_sim_df, dist="weibull")$coefficients
                               ddtx_shape_diab0    <- flexsurvreg(Surv(time = months_to_event, ddtx==1)~1, data=df_diab[[1]], dist="weibull")$coefficients
                               ddtx_shape_diab1    <- flexsurvreg(Surv(time = months_to_event, ddtx==1)~1, data=df_diab[[2]], dist="weibull")$coefficients
                               ddtx_shape_white    <- flexsurvreg(Surv(time = months_to_event, ddtx==1)~1, data=df_race[[1]], dist="weibull")$coefficients
                               ddtx_shape_black    <- flexsurvreg(Surv(time = months_to_event, ddtx==1)~1, data=df_race[[2]], dist="weibull")$coefficients
                               ddtx_shape_hispanic <- flexsurvreg(Surv(time = months_to_event, ddtx==1)~1, data=df_race[[3]], dist="weibull")$coefficients
                               ddtx_shape_other    <- flexsurvreg(Surv(time = months_to_event, ddtx==1)~1, data=df_race[[4]], dist="weibull")$coefficients

                               ### Living donor transplant
                               ldtx_km_fit          <- survfit(Surv(time = months_to_event, ldtx==1) ~ 1, data = list_sim_df)
                               ldtx_km_fit_diab0    <- survfit(Surv(time = months_to_event, ldtx==1) ~ 1, data = df_diab[[1]])
                               ldtx_km_fit_diab1    <- survfit(Surv(time = months_to_event, ldtx==1) ~ 1, data = df_diab[[2]])
                               ldtx_km_fit_white    <- survfit(Surv(time = months_to_event, ldtx==1) ~ 1, data = df_race[[1]])
                               ldtx_km_fit_black    <- survfit(Surv(time = months_to_event, ldtx==1) ~ 1, data = df_race[[2]])
                               ldtx_km_fit_hispanic <- survfit(Surv(time = months_to_event, ldtx==1) ~ 1, data = df_race[[3]])
                               ldtx_km_fit_other    <- survfit(Surv(time = months_to_event, ldtx==1) ~ 1, data = df_race[[4]])

                               ldtx_surv          <- c(summary(ldtx_km_fit, times = seq(from=3, to=180, by=3), extend=TRUE)$surv)
                               ldtx_surv_diab0    <- c(summary(ldtx_km_fit_diab0, times = seq(from=3, to=180, by=3), extend=TRUE)$surv)
                               ldtx_surv_diab1    <- c(summary(ldtx_km_fit_diab1, times = seq(from=3, to=180, by=3), extend=TRUE)$surv)
                               ldtx_surv_white    <- c(summary(ldtx_km_fit_white, times = seq(from=3, to=180, by=3), extend=TRUE)$surv)
                               ldtx_surv_black    <- c(summary(ldtx_km_fit_black, times = seq(from=3, to=180, by=3), extend=TRUE)$surv)
                               ldtx_surv_hispanic <- c(summary(ldtx_km_fit_hispanic, times = seq(from=3, to=180, by=3), extend=TRUE)$surv)
                               ldtx_surv_other    <- c(summary(ldtx_km_fit_other, times = seq(from=3, to=180, by=3), extend=TRUE)$surv)

                               # Parametric Fits
                               ldtx_shape_all      <- flexsurvreg(Surv(time = months_to_event, ldtx==1)~1, data=list_sim_df, dist="gompertz")$coefficients
                               ldtx_shape_diab0    <- flexsurvreg(Surv(time = months_to_event, ldtx==1)~1, data=df_diab[[1]], dist="gompertz")$coefficients
                               ldtx_shape_diab1    <- flexsurvreg(Surv(time = months_to_event, ldtx==1)~1, data=df_diab[[2]], dist="gompertz")$coefficients
                               ldtx_shape_white    <- flexsurvreg(Surv(time = months_to_event, ldtx==1)~1, data=df_race[[1]], dist="gompertz")$coefficients
                               ldtx_shape_black    <- flexsurvreg(Surv(time = months_to_event, ldtx==1)~1, data=df_race[[2]], dist="gompertz")$coefficients
                               ldtx_shape_hispanic <- flexsurvreg(Surv(time = months_to_event, ldtx==1)~1, data=df_race[[3]], dist="gompertz")$coefficients
                               ldtx_shape_other    <- flexsurvreg(Surv(time = months_to_event, ldtx==1)~1, data=df_race[[4]], dist="gompertz")$coefficients

                               ### Waitlist mortality
                               mort_km_fit          <- survfit(Surv(time = months_to_event, mort==1) ~ 1, data = list_sim_df)
                               mort_km_fit_diab0    <- survfit(Surv(time = months_to_event, mort==1) ~ 1, data = df_diab[[1]])
                               mort_km_fit_diab1    <- survfit(Surv(time = months_to_event, mort==1) ~ 1, data = df_diab[[2]])
                               mort_km_fit_white    <- survfit(Surv(time = months_to_event, mort==1) ~ 1, data = df_race[[1]])
                               mort_km_fit_black    <- survfit(Surv(time = months_to_event, mort==1) ~ 1, data = df_race[[2]])
                               mort_km_fit_hispanic <- survfit(Surv(time = months_to_event, mort==1) ~ 1, data = df_race[[3]])
                               mort_km_fit_other    <- survfit(Surv(time = months_to_event, mort==1) ~ 1, data = df_race[[4]])

                               mort_surv          <- c(summary(mort_km_fit, times = seq(from=3, to=180, by=3), extend=TRUE)$surv)
                               mort_surv_diab0    <- c(summary(mort_km_fit_diab0, times = seq(from=3, to=180, by=3), extend=TRUE)$surv)
                               mort_surv_diab1    <- c(summary(mort_km_fit_diab1, times = seq(from=3, to=180, by=3), extend=TRUE)$surv)
                               mort_surv_white    <- c(summary(mort_km_fit_white, times = seq(from=3, to=180, by=3), extend=TRUE)$surv)
                               mort_surv_black    <- c(summary(mort_km_fit_black, times = seq(from=3, to=180, by=3), extend=TRUE)$surv)
                               mort_surv_hispanic <- c(summary(mort_km_fit_hispanic, times = seq(from=3, to=180, by=3), extend=TRUE)$surv)
                               mort_surv_other    <- c(summary(mort_km_fit_other, times = seq(from=3, to=180, by=3), extend=TRUE)$surv)

                               # Parametric Fits
                               mort_shape_all      <- flexsurvreg(Surv(time = months_to_event, mort==1)~1, data=list_sim_df, dist="weibull")$coefficients
                               mort_shape_diab0    <- flexsurvreg(Surv(time = months_to_event, mort==1)~1, data=df_diab[[1]], dist="weibull")$coefficients
                               mort_shape_diab1    <- flexsurvreg(Surv(time = months_to_event, mort==1)~1, data=df_diab[[2]], dist="weibull")$coefficients
                               mort_shape_white    <- flexsurvreg(Surv(time = months_to_event, mort==1)~1, data=df_race[[1]], dist="weibull")$coefficients
                               mort_shape_black    <- flexsurvreg(Surv(time = months_to_event, mort==1)~1, data=df_race[[2]], dist="weibull")$coefficients
                               mort_shape_hispanic <- flexsurvreg(Surv(time = months_to_event, mort==1)~1, data=df_race[[3]], dist="weibull")$coefficients
                               mort_shape_other    <- flexsurvreg(Surv(time = months_to_event, mort==1)~1, data=df_race[[4]], dist="weibull")$coefficients

                               ### Waitlist removal
                               remove_km_fit          <- survfit(Surv(time = months_to_event, remove==1) ~ 1, data = list_sim_df)
                               remove_km_fit_diab0    <- survfit(Surv(time = months_to_event, remove==1) ~ 1, data = df_diab[[1]])
                               remove_km_fit_diab1    <- survfit(Surv(time = months_to_event, remove==1) ~ 1, data = df_diab[[2]])
                               remove_km_fit_white    <- survfit(Surv(time = months_to_event, remove==1) ~ 1, data = df_race[[1]])
                               remove_km_fit_black    <- survfit(Surv(time = months_to_event, remove==1) ~ 1, data = df_race[[2]])
                               remove_km_fit_hispanic <- survfit(Surv(time = months_to_event, remove==1) ~ 1, data = df_race[[3]])
                               remove_km_fit_other    <- survfit(Surv(time = months_to_event, remove==1) ~ 1, data = df_race[[4]])

                               remove_surv          <- c(summary(remove_km_fit, times = seq(from=3, to=180, by=3), extend=TRUE)$surv)
                               remove_surv_diab0    <- c(summary(remove_km_fit_diab0, times = seq(from=3, to=180, by=3), extend=TRUE)$surv)
                               remove_surv_diab1    <- c(summary(remove_km_fit_diab1, times = seq(from=3, to=180, by=3), extend=TRUE)$surv)
                               remove_surv_white    <- c(summary(remove_km_fit_white, times = seq(from=3, to=180, by=3), extend=TRUE)$surv)
                               remove_surv_black    <- c(summary(remove_km_fit_black, times = seq(from=3, to=180, by=3), extend=TRUE)$surv)
                               remove_surv_hispanic <- c(summary(remove_km_fit_hispanic, times = seq(from=3, to=180, by=3), extend=TRUE)$surv)
                               remove_surv_other    <- c(summary(remove_km_fit_other, times = seq(from=3, to=180, by=3), extend=TRUE)$surv)

                               # Parametric Fits
                               remove_shape_all      <- flexsurvreg(Surv(time = months_to_event, remove==1)~1, data=list_sim_df, dist="weibull")$coefficients
                               remove_shape_diab0    <- flexsurvreg(Surv(time = months_to_event, remove==1)~1, data=df_diab[[1]], dist="weibull")$coefficients
                               remove_shape_diab1    <- flexsurvreg(Surv(time = months_to_event, remove==1)~1, data=df_diab[[2]], dist="weibull")$coefficients
                               remove_shape_white    <- flexsurvreg(Surv(time = months_to_event, remove==1)~1, data=df_race[[1]], dist="weibull")$coefficients
                               remove_shape_black    <- flexsurvreg(Surv(time = months_to_event, remove==1)~1, data=df_race[[2]], dist="weibull")$coefficients
                               remove_shape_hispanic <- flexsurvreg(Surv(time = months_to_event, remove==1)~1, data=df_race[[3]], dist="weibull")$coefficients
                               remove_shape_other    <- flexsurvreg(Surv(time = months_to_event, remove==1)~1, data=df_race[[4]], dist="weibull")$coefficients

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
                               v_ddtx_shape <- c(ddtx_shape_all, ddtx_shape_diab0,
                                                 ddtx_shape_diab1, ddtx_shape_white,
                                                 ddtx_shape_black, ddtx_shape_hispanic,
                                                 ddtx_shape_other)
                               v_ldtx_shape <- c(ldtx_shape_all, ldtx_shape_diab0, ldtx_shape_diab1,
                                                 ldtx_shape_white, ldtx_shape_black,
                                                 ldtx_shape_hispanic, ldtx_shape_other)
                               v_mort_shape <- c(mort_shape_all, mort_shape_diab0, mort_shape_diab1,
                                                 mort_shape_white, mort_shape_black,
                                                 mort_shape_hispanic, mort_shape_other)
                               v_remove_shape <- c(remove_shape_all, remove_shape_diab0, remove_shape_diab1,
                                                   remove_shape_white, remove_shape_black,
                                                   remove_shape_hispanic, remove_shape_other)


                               ### Assess post-tx fit
                               sim_df <- sim_df %>% filter(ddtx==1)

                               ### 30-day Outcomes
                               prop_30day_all <- rep(0, 4)
                               prop_30day_diab0 <- rep(0,4)
                               prop_30day_diab1 <- rep(0,4)
                               prop_30day_white <- rep(0,4)
                               prop_30day_black <- rep(0,4)
                               prop_30day_hispanic <- rep(0,4)
                               prop_30day_other <- rep(0,4)
                               for(i in 0:3) {
                                 prop_30day_all[i+1] <- nrow(sim_df[sim_df$tx_outcome==i,])
                                 prop_30day_diab0[i+1] <- nrow(sim_df[sim_df$tx_outcome==i & sim_df$diab_stat==0,])
                                 prop_30day_diab1[i+1] <- nrow(sim_df[sim_df$tx_outcome==i & sim_df$diab_stat==1,])
                                 prop_30day_white[i+1] <- nrow(sim_df[sim_df$tx_outcome==i & sim_df$race==0,])
                                 prop_30day_black[i+1] <- nrow(sim_df[sim_df$tx_outcome==i & sim_df$race==1,])
                                 prop_30day_hispanic[i+1] <- nrow(sim_df[sim_df$tx_outcome==i & sim_df$race==2,])
                                 prop_30day_other[i+1] <- nrow(sim_df[sim_df$tx_outcome==i & sim_df$race==3,])
                               }
                               prop_30day_all <- prop_30day_all/sum(prop_30day_all)
                               prop_30day_diab0 <- prop_30day_diab0/sum(prop_30day_diab0)
                               prop_30day_diab1 <- prop_30day_diab1/sum(prop_30day_diab1)
                               prop_30day_white <- prop_30day_white/sum(prop_30day_white)
                               prop_30day_black <- prop_30day_black/sum(prop_30day_black)
                               prop_30day_hispanic <- prop_30day_hispanic/sum(prop_30day_hispanic)
                               prop_30day_other <- prop_30day_other/sum(prop_30day_other)

                               # cols30 <- c("0", "1", "2", "3")
                               # prop_30day       <- matrix(nrow = 7, ncol = 4)
                               # if(length(table(sim_df$tx_outcome)) == 3){
                               #   names_prop30 <- names(table(sim_df$tx_outcome))
                               #   names_diff   <- setdiff(cols30, names_prop30)
                               #   all_prop30 <- c(table(sim_df$tx_outcome),
                               #                   names_diff = 0)
                               #   names(all_prop30) <- c(names_prop30, names_diff)
                               #   all_prop30 <- all_prop30[c("0", "1", "2", "3")]
                               #   prop_30day[1,]   <- as.vector(all_prop30)
                               #
                               #   diab_prop30 <- cbind(table(sim_df$diab_stat, sim_df$tx_outcome),
                               #                        names_diff = c(0,0))
                               #   colnames(diab_prop30) <- c(names_prop30, names_diff)
                               #   diab_prop30 <- diab_prop30[, c("0", "1", "2", "3")]
                               #   prop_30day[2:3,] <- diab_prop30
                               #
                               #   race_prop30 <- cbind(table(sim_df$race, sim_df$tx_outcome),
                               #                        names_diff = c(0,0,0,0))
                               #   colnames(race_prop30) <- c(names_prop30, names_diff)
                               #   race_prop30 <- race_prop30[, c("0", "1", "2", "3")]
                               #   prop_30day[4:7,] <- race_prop30
                               # } else {
                               #   prop_30day[1,]   <- as.vector(table(sim_df$tx_outcome))
                               #   prop_30day[2:3,] <- table(sim_df$diab_stat, sim_df$tx_outcome)
                               #   prop_30day[4:7,] <- table(sim_df$race, sim_df$tx_outcome)
                               # }

                               # make rows proportions
                               # prop_30day <- t(apply(prop_30day,1, function(x) x/sum(x)))
                               # prop_30day <- c(prop_30day)
                               prop_30day <- c(prop_30day_all,
                                               prop_30day_diab0,
                                               prop_30day_diab1,
                                               prop_30day_white,
                                               prop_30day_black,
                                               prop_30day_hispanic,
                                               prop_30day_other)

                               ## Split data
                               df_diab <- split(sim_df, sim_df$diab_stat)
                               df_race <- split(sim_df, sim_df$race)

                               ### Death with function
                               # KM Curves
                               gs_mort_fit          <- survival::survfit(Surv(time = gs_death_time, gs_death==1) ~ 1, data = sim_df)
                               gs_mort_fit_diab0    <- survival::survfit(Surv(time = gs_death_time, gs_death==1) ~ 1, data = df_diab[[1]])
                               gs_mort_fit_diab1    <- survival::survfit(Surv(time = gs_death_time, gs_death==1) ~ 1, data = df_diab[[2]])
                               gs_mort_fit_white    <- survival::survfit(Surv(time = gs_death_time, gs_death==1) ~ 1, data = df_race[[1]])
                               gs_mort_fit_black    <- survival::survfit(Surv(time = gs_death_time, gs_death==1) ~ 1, data = df_race[[2]])
                               gs_mort_fit_hispanic <- survival::survfit(Surv(time = gs_death_time, gs_death==1) ~ 1, data = df_race[[3]])
                               gs_mort_fit_other    <- survival::survfit(Surv(time = gs_death_time, gs_death==1) ~ 1, data = df_race[[4]])

                               gs_mort_surv          <- c(summary(gs_mort_fit, times = seq(from=3, to=120, by=3), extend=TRUE)$surv)
                               gs_mort_surv_diab0    <- c(summary(gs_mort_fit_diab0, times = seq(from=3, to=120, by=3), extend=TRUE)$surv)
                               gs_mort_surv_diab1    <- c(summary(gs_mort_fit_diab1, times = seq(from=3, to=120, by=3), extend=TRUE)$surv)
                               gs_mort_surv_white    <- c(summary(gs_mort_fit_white, times = seq(from=3, to=120, by=3), extend=TRUE)$surv)
                               gs_mort_surv_black    <- c(summary(gs_mort_fit_black, times = seq(from=3, to=120, by=3), extend=TRUE)$surv)
                               gs_mort_surv_hispanic <- c(summary(gs_mort_fit_hispanic, times = seq(from=3, to=120, by=3), extend=TRUE)$surv)
                               gs_mort_surv_other    <- c(summary(gs_mort_fit_other, times = seq(from=3, to=120, by=3), extend=TRUE)$surv)

                               # Parametric Fits
                               gs_mort_shape_all      <- flexsurv::flexsurvreg(Surv(time = gs_death_time, gs_death==1) ~ 1,
                                                                               inits = c(shape = 0.001, rate = 0.001), data=sim_df, dist="gompertz")$coefficients
                               gs_mort_shape_diab0    <- flexsurv::flexsurvreg(Surv(time = gs_death_time, gs_death==1) ~ 1,
                                                                               inits = c(shape = 0.001, rate = 0.001), data=df_diab[[1]], dist="gompertz")$coefficients
                               gs_mort_shape_diab1    <- flexsurv::flexsurvreg(Surv(time = gs_death_time, gs_death==1) ~ 1,
                                                                               inits = c(shape = 0.001, rate = 0.001), data=df_diab[[2]], dist="gompertz")$coefficients
                               gs_mort_shape_white    <- flexsurv::flexsurvreg(Surv(time = gs_death_time, gs_death==1) ~ 1,
                                                                               inits = c(shape = 0.001, rate = 0.001), data=df_race[[1]], dist="gompertz")$coefficients
                               gs_mort_shape_black    <- flexsurv::flexsurvreg(Surv(time = gs_death_time, gs_death==1) ~ 1,
                                                                               inits = c(shape = 0.001, rate = 0.001), data=df_race[[2]], dist="gompertz")$coefficients
                               gs_mort_shape_hispanic <- flexsurv::flexsurvreg(Surv(time = gs_death_time, gs_death==1) ~ 1,
                                                                               inits = c(shape = 0.001, rate = 0.001), data=df_race[[3]], dist="gompertz")$coefficients
                               gs_mort_shape_other    <- flexsurv::flexsurvreg(Surv(time = gs_death_time, gs_death==1) ~ 1,
                                                                               inits = c(shape = 0.001, rate = 0.001), data=df_race[[4]], dist="gompertz")$coefficients
                               v_gs_mort_shape <- c(gs_mort_shape_all,
                                                    gs_mort_shape_diab0,
                                                    gs_mort_shape_diab1,
                                                    gs_mort_shape_white,
                                                    gs_mort_shape_black,
                                                    gs_mort_shape_hispanic,
                                                    gs_mort_shape_other)
                               ### Graft loss
                               gl_fit          <- survival::survfit(Surv(time = gl_time, graft_loss==1) ~ 1, data = sim_df)
                               gl_fit_diab0    <- survival::survfit(Surv(time = gl_time, graft_loss==1) ~ 1, data = df_diab[[1]])
                               gl_fit_diab1    <- survival::survfit(Surv(time = gl_time, graft_loss==1) ~ 1, data = df_diab[[2]])
                               gl_fit_white    <- survival::survfit(Surv(time = gl_time, graft_loss==1) ~ 1, data = df_race[[1]])
                               gl_fit_black    <- survival::survfit(Surv(time = gl_time, graft_loss==1) ~ 1, data = df_race[[2]])
                               gl_fit_hispanic <- survival::survfit(Surv(time = gl_time, graft_loss==1) ~ 1, data = df_race[[3]])
                               gl_fit_other    <- survival::survfit(Surv(time = gl_time, graft_loss==1) ~ 1, data = df_race[[4]])

                               gl_surv          <- c(summary(gl_fit, times = seq(from=3, to=120, by=3), extend=TRUE)$surv)
                               gl_surv_diab0    <- c(summary(gl_fit_diab0, times = seq(from=3, to=120, by=3), extend=TRUE)$surv)
                               gl_surv_diab1    <- c(summary(gl_fit_diab1, times = seq(from=3, to=120, by=3), extend=TRUE)$surv)
                               gl_surv_white    <- c(summary(gl_fit_white, times = seq(from=3, to=120, by=3), extend=TRUE)$surv)
                               gl_surv_black    <- c(summary(gl_fit_black, times = seq(from=3, to=120, by=3), extend=TRUE)$surv)
                               gl_surv_hispanic <- c(summary(gl_fit_hispanic, times = seq(from=3, to=120, by=3), extend=TRUE)$surv)
                               gl_surv_other    <- c(summary(gl_fit_other, times = seq(from=3, to=120, by=3), extend=TRUE)$surv)

                               # Parametric Fits
                               gl_shape_all      <- flexsurv::flexsurvreg(Surv(time = gl_time, graft_loss==1) ~ 1,
                                                                          inits = c(shape = 0.001, rate = 0.001), data=sim_df, dist="gompertz")$coefficients
                               gl_shape_diab0    <- flexsurv::flexsurvreg(Surv(time = gl_time, graft_loss==1) ~ 1,
                                                                          inits = c(shape = 0.001, rate = 0.001), data=df_diab[[1]], dist="gompertz")$coefficients
                               gl_shape_diab1    <- flexsurv::flexsurvreg(Surv(time = gl_time, graft_loss==1) ~ 1,
                                                                          inits = c(shape = 0.001, rate = 0.001), data=df_diab[[2]], dist="gompertz")$coefficients
                               gl_shape_white    <- flexsurv::flexsurvreg(Surv(time = gl_time, graft_loss==1) ~ 1,
                                                                          inits = c(shape = 0.001, rate = 0.001), data=df_race[[1]], dist="gompertz")$coefficients
                               gl_shape_black    <- flexsurv::flexsurvreg(Surv(time = gl_time, graft_loss==1) ~ 1,
                                                                          inits = c(shape = 0.001, rate = 0.001), data=df_race[[2]], dist="gompertz")$coefficients
                               gl_shape_hispanic <- flexsurv::flexsurvreg(Surv(time = gl_time, graft_loss==1) ~ 1,
                                                                          inits = c(shape = 0.001, rate = 0.001), data=df_race[[3]], dist="gompertz")$coefficients
                               gl_shape_other    <- flexsurv::flexsurvreg(Surv(time = gl_time, graft_loss==1) ~ 1,
                                                                          inits = c(shape = 0.001, rate = 0.001), data=df_race[[4]], dist="gompertz")$coefficients
                               v_gl_shape <- c(gl_shape_all,
                                               gl_shape_diab0,
                                               gl_shape_diab1,
                                               gl_shape_white,
                                               gl_shape_black,
                                               gl_shape_hispanic,
                                               gl_shape_other)

                               # Death after graft loss
                               tryCatch(expr = {
                                 gl_mort_fit  <- survival::survfit(Surv(time = gl_death_time, gl_death==1) ~ 1, data = sim_df)
                                 gl_mort_surv <- c(summary(gl_mort_fit, times = seq(from=3, to=60, by=3), extend=TRUE)$surv)
                               }, error = function(e) {gl_mort_surv <- rep(NA, 20)}
                               )

                               tryCatch(expr = {
                                 gl_mort_fit_diab0  <- survival::survfit(Surv(time = gl_death_time, gl_death==1) ~ 1, data = df_diab[[1]])
                                 gl_mort_surv_diab0 <- c(summary(gl_mort_fit_diab0, times = seq(from=3, to=60, by=3), extend=TRUE)$surv)
                               }, error = function(e) {gl_mort_surv_diab0 <- rep(NA, 20)}
                               )

                               tryCatch(expr = {
                                 gl_mort_fit_diab1  <- survival::survfit(Surv(time = gl_death_time, gl_death==1) ~ 1, data = df_diab[[2]])
                                 gl_mort_surv_diab1 <- c(summary(gl_mort_fit_diab0, times = seq(from=3, to=60, by=3), extend=TRUE)$surv)
                               }, error = function(e) {gl_mort_surv_diab1 <- rep(NA, 20)}
                               )

                               tryCatch(expr = {
                                 gl_mort_fit_white  <- survival::survfit(Surv(time = gl_death_time, gl_death==1) ~ 1, data = df_race[[1]])
                                 gl_mort_surv_white <- c(summary(gl_mort_fit_white, times = seq(from=3, to=60, by=3), extend=TRUE)$surv)
                               }, error = function(e) {gl_mort_surv_white <- rep(NA, 20)}
                               )

                               tryCatch(expr = {
                                 gl_mort_fit_black  <- survival::survfit(Surv(time = gl_death_time, gl_death==1) ~ 1, data = df_race[[2]])
                                 gl_mort_surv_black <- c(summary(gl_mort_fit_white, times = seq(from=3, to=60, by=3), extend=TRUE)$surv)
                               }, error = function(e) {gl_mort_surv_black <- rep(NA, 20)}
                               )

                               tryCatch(expr = {
                                 gl_mort_fit_hispanic  <- survival::survfit(Surv(time = gl_death_time, gl_death==1) ~ 1, data = df_race[[3]])
                                 gl_mort_surv_hispanic <- c(summary(gl_mort_fit_white, times = seq(from=3, to=60, by=3), extend=TRUE)$surv)
                               }, error = function(e) {gl_mort_surv_hispanic <- rep(NA, 20)}
                               )

                               tryCatch(expr = {
                                 gl_mort_fit_other  <- survival::survfit(Surv(time = gl_death_time, gl_death==1) ~ 1, data = df_race[[4]])
                                 gl_mort_surv_other <- c(summary(gl_mort_fit_white, times = seq(from=3, to=60, by=3), extend=TRUE)$surv)
                               }, error = function(e) {gl_mort_surv_other <- rep(NA, 20)}
                               )

                               # Parametric Fits
                               tryCatch(expr = {
                                 sim_df_weibull <- sim_df %>%
                                   filter(graft_loss==1 & gl_death_time > 0)
                                 df_diab_weibull <- split(sim_df_weibull, sim_df_weibull$diab_stat)
                                 df_race_weibull <- split(sim_df_weibull, sim_df_weibull$race)
                               }, error = function(e) {v_gl_mort_shape <- rep(NA, 14)}
                               )


                               gl_mort_shape_all <- parametric.fit(data = sim_df_weibull,
                                                                   distribution = "weibull",
                                                                   t = sim_df_weibull$gl_death_time,
                                                                   outcome = sim_df_weibull$gl_death)
                               gl_mort_shape_diab0 <- parametric.fit(data = df_diab_weibull[[1]],
                                                                     distribution = "weibull",
                                                                     t = df_diab_weibull[[1]]$gl_death_time,
                                                                     outcome = df_diab_weibull[[1]]$gl_death)
                               gl_mort_shape_diab1 <- parametric.fit(data = df_diab_weibull[[2]],
                                                                     distribution = "weibull",
                                                                     t = df_diab_weibull[[2]]$gl_death_time,
                                                                     outcome = df_diab_weibull[[2]]$gl_death)
                               gl_mort_shape_white <- parametric.fit(data = df_race_weibull[[1]],
                                                                     distribution = "weibull",
                                                                     t = df_race_weibull[[1]]$gl_death_time,
                                                                     outcome = df_race_weibull[[1]]$gl_death)
                               gl_mort_shape_black <- parametric.fit(data = df_race_weibull[[2]],
                                                                     distribution = "weibull",
                                                                     t = df_race_weibull[[2]]$gl_death_time,
                                                                     outcome = df_race_weibull[[2]]$gl_death)
                               gl_mort_shape_hispanic <- parametric.fit(data = df_race_weibull[[3]],
                                                                        distribution = "weibull",
                                                                        t = df_race_weibull[[3]]$gl_death_time,
                                                                        outcome = df_race_weibull[[3]]$gl_death)
                               gl_mort_shape_other <- parametric.fit(data = df_race_weibull[[4]],
                                                                     distribution = "weibull",
                                                                     t = df_race_weibull[[4]]$gl_death_time,
                                                                     outcome = df_race_weibull[[4]]$gl_death)

                               v_gl_mort_shape <- c(gl_mort_shape_all,
                                                    gl_mort_shape_diab0,
                                                    gl_mort_shape_diab1,
                                                    gl_mort_shape_white,
                                                    gl_mort_shape_black,
                                                    gl_mort_shape_hispanic,
                                                    gl_mort_shape_other)

                               ### Death same day as graft loss
                               tryCatch(expr = {
                                 df_gl <- sim_df %>%
                                   filter(graft_loss==1) %>%
                                   mutate(death_at_gl = ifelse(gl_death==1 & gl_death_time==0, 1, 0))

                                 # Create empty vector (with names)
                                 prop_gl_mort <- rep(0, 7)
                                 # names(prop_gl_mort) <- c("All", "No Diabetes", "Yes Diabetes",
                                 #                          "White", "Black", "Hispanic", 'Other')

                                 # All graft loss patients
                                 prop_gl_mort[1] <- nrow(df_gl[df_gl$death_at_gl==1,])/nrow(df_gl)
                                 # By diabetes status
                                 prop_gl_mort[2] <- nrow(df_gl[df_gl$death_at_gl==1 &
                                                                 df_gl$diab_stat==0,])/nrow(df_gl)
                                 prop_gl_mort[3] <- nrow(df_gl[df_gl$death_at_gl==1 &
                                                                 df_gl$diab_stat==1,])/nrow(df_gl)
                                 # By race/ethnicity
                                 prop_gl_mort[4] <- nrow(df_gl[df_gl$death_at_gl==1 &
                                                                 df_gl$race==0,])/nrow(df_gl)
                                 prop_gl_mort[5] <- nrow(df_gl[df_gl$death_at_gl==1 &
                                                                 df_gl$race==1,])/nrow(df_gl)
                                 prop_gl_mort[6] <- nrow(df_gl[df_gl$death_at_gl==1 &
                                                                 df_gl$race==2,])/nrow(df_gl)
                                 prop_gl_mort[7] <- nrow(df_gl[df_gl$death_at_gl==1 &
                                                                 df_gl$race==3,])/nrow(df_gl)
                               },
                               error = function(e) {prop_gl_mort <- rep(NA, times = 7)}
                               )

                               # Combine to export
                               list(ddtx_surv, ldtx_surv,
                                    mort_surv, remove_surv,
                                    ddtx_surv_diab0, ldtx_surv_diab0,
                                    mort_surv_diab0, remove_surv_diab0,
                                    ddtx_surv_diab1, ldtx_surv_diab1,
                                    mort_surv_diab1, remove_surv_diab1,
                                    ddtx_surv_white, ldtx_surv_white,
                                    mort_surv_white, remove_surv_white,
                                    ddtx_surv_black, ldtx_surv_black,
                                    mort_surv_black, remove_surv_black,
                                    ddtx_surv_hispanic, ldtx_surv_hispanic,
                                    mort_surv_hispanic, remove_surv_hispanic,
                                    ddtx_surv_other, ldtx_surv_other,
                                    mort_surv_other, remove_surv_other,
                                    v_ddtx_shape, v_ldtx_shape,
                                    v_mort_shape, v_remove_shape,
                                    sim_table1, prop_30day,
                                    gs_mort_surv, gs_mort_surv_diab0,
                                    gs_mort_surv_diab1, gs_mort_surv_white,
                                    gs_mort_surv_black, gs_mort_surv_hispanic,
                                    gs_mort_surv_other, gl_surv,
                                    gl_surv_diab0, gl_surv_diab1,
                                    gl_surv_white, gl_surv_black,
                                    gl_surv_hispanic, gl_surv_other,
                                    gl_mort_surv, gl_mort_surv_diab0,
                                    gl_mort_surv_diab1, gl_mort_surv_white,
                                    gl_mort_surv_black, gl_mort_surv_hispanic,
                                    gl_mort_surv_other, v_gs_mort_shape,
                                    v_gl_shape, v_gl_mort_shape,
                                    prop_gl_mort)

                             }
parallel::stopCluster(cl = my.cluster)

for(i in 1:length(out)){
  fwrite(out[[i]], paste0(out_path, "Stage2/part2/out_", i, ".csv"))
}

l_target_names_wait <- list("target_1.csv", "target_2.csv", "target_3.csv",
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
l_target_names_post <- list("target_34.csv", "target_35.csv", "target_36.csv",
                            "target_37.csv", "target_38.csv", "target_39.csv",
                            "target_40.csv", "target_41.csv", "target_42.csv",
                            "target_43.csv", "target_44.csv", "target_45.csv",
                            "target_46.csv", "target_47.csv", "target_48.csv",
                            "target_49.csv", "target_50.csv", "target_51.csv",
                            "target_52.csv", "target_53.csv", "target_54.csv",
                            "target_55.csv", "target_56.csv", "target_57.csv",
                            "target_58.csv", "target_59.csv")


l_target <- list()
for (k in 1:length(l_target_names_wait)){
  l_target[[k]] <- read.csv(paste0(out_path, "waitlist/", l_target_names_wait[[k]]))[-1]
}
j <- 0
for (k in 34:59){
  j <- j+1
  l_target[[k]] <- read.csv(paste0(out_path, "Stage2/", l_target_names_post[[j]]))[-1]
}


m_gof <- matrix(nrow = n, ncol = 14)
colnames(m_gof) <- c("id", "gof_wait_surv", "gof_wait_shape", "gof_table1",
                     "gof_gs_mort_surv", "gof_gl_surv", "gof_gl_mort_surv",
                     "gof_gs_mort_shape", "gof_gl_shape", "gof_gl_mort_shape",
                     "gof_post_surv", "gof_post_shape", "gof_30", "gof_gl_mort_prop")
pb <- txtProgressBar(min = 0, max = n, style = 3)
for (i in 1:n) {
  m_gof[i,1] <- i
  l_gof <- get_gof_all_likelihood(targets = l_target, sim_results = out, iter = i)
  m_gof[i,2] <- l_gof[[1]]
  m_gof[i,3] <- l_gof[[2]]
  m_gof[i,4] <- l_gof[[3]]
  m_gof[i,5] <- l_gof[[4]]
  m_gof[i,6] <- l_gof[[5]]
  m_gof[i,7] <- l_gof[[6]]
  m_gof[i,8] <- l_gof[[7]]
  m_gof[i,9] <- l_gof[[8]]
  m_gof[i,10] <- l_gof[[9]]
  m_gof[i,11] <- l_gof[[10]]
  m_gof[i,12] <- l_gof[[11]]
  m_gof[i,13] <- l_gof[[12]]
  m_gof[i,14] <- l_gof[[13]]

  setTxtProgressBar(pb, i)
}

fwrite(m_gof, paste0(out_path, "Stage2/part2/m_gof.csv"))






", "gof_wait_surv", "gof_wait_shape", "gof_table1",
                     "gof_gs_mort_surv", "gof_gl_surv", "gof_gl_mort_surv",
                     "gof_gs_mort_shape", "gof_gl_shape", "gof_gl_mort_shape",
                     "gof_post_surv", "gof_post_shape", "gof_30", "gof_gl_mort_prop")
pb <- txtProgressBar(min = 0, max = 1000, style = 3)
for (i in 1:1000) {
  m_gof[i,1] <- i
  l_gof <- get_gof_all_likelihood(targets = l_target, sim_results = out, iter = i)
  m_gof[i,2] <- l_gof[[1]]
  m_gof[i,3] <- l_gof[[2]]
  m_gof[i,4] <- l_gof[[3]]
  m_gof[i,5] <- l_gof[[4]]
  m_gof[i,6] <- l_gof[[5]]
  m_gof[i,7] <- l_gof[[6]]
  m_gof[i,8] <- l_gof[[7]]
  m_gof[i,9] <- l_gof[[8]]
  m_gof[i,10] <- l_gof[[9]]
  m_gof[i,11] <- l_gof[[10]]
  m_gof[i,12] <- l_gof[[11]]
  m_gof[i,13] <- l_gof[[12]]
  m_gof[i,14] <- l_gof[[13]]

  setTxtProgressBar(pb, i)
}






