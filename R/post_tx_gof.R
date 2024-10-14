get_gof_all_likelihood <- function(targets, sim_results, iter) {

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

  ### 30-Day Outcomes
  gof_30_all      <- LaplacesDemon::ddirichlet(sim_results[[34]][iter,1:4], as.numeric(targets[[34]][1,]), log = TRUE)
  gof_30_diab0    <- LaplacesDemon::ddirichlet(sim_results[[34]][iter,5:8], as.numeric(targets[[34]][2,]), log = TRUE)
  gof_30_diab1    <- LaplacesDemon::ddirichlet(sim_results[[34]][iter,9:12], as.numeric(targets[[34]][3,]), log = TRUE)
  gof_30_white    <- LaplacesDemon::ddirichlet(sim_results[[34]][iter,13:16], as.numeric(targets[[34]][4,]), log = TRUE)
  gof_30_black    <- LaplacesDemon::ddirichlet(sim_results[[34]][iter,17:20], as.numeric(targets[[34]][5,]), log = TRUE)
  gof_30_hispanic <- LaplacesDemon::ddirichlet(sim_results[[34]][iter,21:24], as.numeric(targets[[34]][6,]), log = TRUE)
  gof_30_other    <- LaplacesDemon::ddirichlet(sim_results[[34]][iter,25:28], as.numeric(targets[[34]][7,]), log = TRUE)
  gof_30 <- -1 * sum(gof_30_all, gof_30_diab0, gof_30_diab1, gof_30_white,
                     gof_30_black, gof_30_hispanic, gof_30_other)
  ### Death with Function
  ## Survival-based GOF
  gof_gs_mort_all <- dnorm(sim_results[[35]][iter,],
                           mean = as.numeric(targets[[35]][1,]),
                           sd = as.numeric(targets[[35]][2,]), log = T)
  gof_gs_mort_diab0 <- dnorm(sim_results[[36]][iter,],
                             mean = as.numeric(targets[[36]][1,]),
                             sd = as.numeric(targets[[36]][2,]), log = T)
  gof_gs_mort_diab1 <- dnorm(sim_results[[37]][iter,],
                             mean = as.numeric(targets[[37]][1,]),
                             sd = as.numeric(targets[[37]][2,]), log = T)
  gof_gs_mort_white <- dnorm(sim_results[[38]][iter,],
                             mean = as.numeric(targets[[38]][1,]),
                             sd = as.numeric(targets[[38]][2,]), log = T)
  gof_gs_mort_black <- dnorm(sim_results[[39]][iter,],
                             mean = as.numeric(targets[[39]][1,]),
                             sd = as.numeric(targets[[39]][2,]), log = T)
  gof_gs_mort_hispanic <- dnorm(sim_results[[40]][iter,],
                                mean = as.numeric(targets[[40]][1,]),
                                sd = as.numeric(targets[[40]][2,]), log = T)
  gof_gs_mort_other <- dnorm(sim_results[[41]][iter,],
                             mean = as.numeric(targets[[41]][1,]),
                             sd = as.numeric(targets[[41]][2,]), log = T)

  ## Shape-based GOF
  gof_gs_mort_shape_all <- dmvnorm(sim_results[[56]][iter,1:2],
                                   mean = as.numeric(targets[[42]][1,1:2]),
                                   sigma = as.matrix(targets[[42]][2:3,1:2]), log = T)
  gof_gs_mort_shape_diab0 <- dmvnorm(sim_results[[56]][iter,3:4],
                                     mean = as.numeric(targets[[42]][1,3:4]),
                                     sigma = as.matrix(targets[[42]][2:3,3:4]), log = T)
  gof_gs_mort_shape_diab1 <- dmvnorm(sim_results[[56]][iter,5:6],
                                     mean = as.numeric(targets[[42]][1,5:6]),
                                     sigma = as.matrix(targets[[42]][2:3,5:6]), log = T)
  gof_gs_mort_shape_white <- dmvnorm(sim_results[[56]][iter,7:8],
                                     mean = as.numeric(targets[[42]][1,7:8]),
                                     sigma = as.matrix(targets[[42]][2:3,7:8]), log = T)
  gof_gs_mort_shape_black <- dmvnorm(sim_results[[56]][iter,9:10],
                                     mean = as.numeric(targets[[42]][1,9:10]),
                                     sigma = as.matrix(targets[[42]][2:3,9:10]), log = T)
  gof_gs_mort_shape_hispanic <- dmvnorm(sim_results[[56]][iter,11:12],
                                        mean = as.numeric(targets[[42]][1,11:12]),
                                        sigma = as.matrix(targets[[42]][2:3,11:12]), log = T)
  gof_gs_mort_shape_other <- dmvnorm(sim_results[[56]][iter,13:14],
                                     mean = as.numeric(targets[[42]][1,13:14]),
                                     sigma = as.matrix(targets[[42]][2:3,13:14]), log = T)

  ### Graft Loss
  ## Survival-based GOF
  gof_gl_all <- dnorm(sim_results[[42]][iter,],
                      mean = as.numeric(targets[[43]][1,]),
                      sd = as.numeric(targets[[43]][2,]), log = T)
  gof_gl_diab0 <- dnorm(sim_results[[43]][iter,],
                        mean = as.numeric(targets[[44]][1,]),
                        sd = as.numeric(targets[[44]][2,]), log = T)
  gof_gl_diab1 <- dnorm(sim_results[[44]][iter,],
                        mean = as.numeric(targets[[45]][1,]),
                        sd = as.numeric(targets[[45]][2,]), log = T)
  gof_gl_white <- dnorm(sim_results[[45]][iter,],
                        mean = as.numeric(targets[[46]][1,]),
                        sd = as.numeric(targets[[46]][2,]), log = T)
  gof_gl_black <- dnorm(sim_results[[46]][iter,],
                        mean = as.numeric(targets[[47]][1,]),
                        sd = as.numeric(targets[[47]][2,]), log = T)
  gof_gl_hispanic <- dnorm(sim_results[[47]][iter,],
                           mean = as.numeric(targets[[48]][1,]),
                           sd = as.numeric(targets[[48]][2,]), log = T)
  gof_gl_other <- dnorm(sim_results[[48]][iter,],
                        mean = as.numeric(targets[[49]][1,]),
                        sd = as.numeric(targets[[49]][2,]), log = T)

  ## Shape-based GOF
  gof_gl_shape_all <- dmvnorm(sim_results[[57]][iter,1:2],
                              mean = as.numeric(targets[[50]][1,1:2]),
                              sigma = as.matrix(targets[[50]][2:3,1:2]), log = T)
  gof_gl_shape_diab0 <- dmvnorm(sim_results[[57]][iter,3:4],
                                mean = as.numeric(targets[[50]][1,3:4]),
                                sigma = as.matrix(targets[[50]][2:3,3:4]), log = T)
  gof_gl_shape_diab1 <- dmvnorm(sim_results[[57]][iter,5:6],
                                mean = as.numeric(targets[[50]][1,5:6]),
                                sigma = as.matrix(targets[[50]][2:3,5:6]), log = T)
  gof_gl_shape_white <- dmvnorm(sim_results[[57]][iter,7:8],
                                mean = as.numeric(targets[[50]][1,7:8]),
                                sigma = as.matrix(targets[[50]][2:3,7:8]), log = T)
  gof_gl_shape_black <- dmvnorm(sim_results[[57]][iter,9:10],
                                mean = as.numeric(targets[[50]][1,9:10]),
                                sigma = as.matrix(targets[[50]][2:3,9:10]), log = T)
  gof_gl_shape_hispanic <- dmvnorm(sim_results[[57]][iter,11:12],
                                   mean = as.numeric(targets[[50]][1,11:12]),
                                   sigma = as.matrix(targets[[50]][2:3,11:12]), log = T)
  gof_gl_shape_other <- dmvnorm(sim_results[[57]][iter,13:14],
                                mean = as.numeric(targets[[50]][1,13:14]),
                                sigma = as.matrix(targets[[50]][2:3,13:14]), log = T)

  ### Death after Graft Loss
  ## Survival-based GOF
  gof_gl_mort_all <- dnorm(sim_results[[49]][iter,],
                           mean = as.numeric(targets[[51]][1,]),
                           sd = as.numeric(targets[[51]][2,]), log = T)
  gof_gl_mort_diab0 <- dnorm(sim_results[[50]][iter,],
                             mean = as.numeric(targets[[52]][1,]),
                             sd = as.numeric(targets[[52]][2,]), log = T)
  gof_gl_mort_diab1 <- dnorm(sim_results[[51]][iter,],
                             mean = as.numeric(targets[[53]][1,]),
                             sd = as.numeric(targets[[53]][2,]), log = T)
  gof_gl_mort_white <- dnorm(sim_results[[52]][iter,],
                             mean = as.numeric(targets[[54]][1,]),
                             sd = as.numeric(targets[[54]][2,]), log = T)
  gof_gl_mort_black <- dnorm(sim_results[[53]][iter,],
                             mean = as.numeric(targets[[55]][1,]),
                             sd = as.numeric(targets[[55]][2,]), log = T)
  gof_gl_mort_hispanic <- dnorm(sim_results[[54]][iter,],
                                mean = as.numeric(targets[[56]][1,]),
                                sd = as.numeric(targets[[56]][2,]), log = T)
  gof_gl_mort_other <- dnorm(sim_results[[55]][iter,],
                             mean = as.numeric(targets[[57]][1,]),
                             sd = as.numeric(targets[[57]][2,]), log = T)

  ## Shape-based GOF
  gof_gl_mort_shape_all <- dmvnorm(sim_results[[58]][iter,1:2],
                                   mean = as.numeric(targets[[58]][1,1:2]),
                                   sigma = as.matrix(targets[[58]][2:3,1:2]), log = T)
  gof_gl_mort_shape_diab0 <- dmvnorm(sim_results[[58]][iter,3:4],
                                     mean = as.numeric(targets[[58]][1,3:4]),
                                     sigma = as.matrix(targets[[58]][2:3,3:4]), log = T)
  gof_gl_mort_shape_diab1 <- dmvnorm(sim_results[[58]][iter,5:6],
                                     mean = as.numeric(targets[[58]][1,5:6]),
                                     sigma = as.matrix(targets[[58]][2:3,5:6]), log = T)
  gof_gl_mort_shape_white <- dmvnorm(sim_results[[58]][iter,7:8],
                                     mean = as.numeric(targets[[58]][1,7:8]),
                                     sigma = as.matrix(targets[[58]][2:3,7:8]), log = T)
  gof_gl_mort_shape_black <- dmvnorm(sim_results[[58]][iter,9:10],
                                     mean = as.numeric(targets[[58]][1,9:10]),
                                     sigma = as.matrix(targets[[58]][2:3,9:10]), log = T)
  gof_gl_mort_shape_hispanic <- dmvnorm(sim_results[[58]][iter,11:12],
                                        mean = as.numeric(targets[[58]][1,11:12]),
                                        sigma = as.matrix(targets[[58]][2:3,11:12]), log = T)
  gof_gl_mort_shape_other <- dmvnorm(sim_results[[58]][iter,13:14],
                                     mean = as.numeric(targets[[58]][1,13:14]),
                                     sigma = as.matrix(targets[[58]][2:3,13:14]), log = T)

  ### Death same day as graft loss
  gof_prop_all      <- dbeta(sim_results[[59]][iter,1], shape1=targets[[59]][1,1], shape2=targets[[59]][1,2], log = T)
  gof_prop_diab0    <- dbeta(sim_results[[59]][iter,2], shape1=targets[[59]][2,1], shape2=targets[[59]][2,2], log = T)
  gof_prop_diab1    <- dbeta(sim_results[[59]][iter,3], shape1=targets[[59]][3,1], shape2=targets[[59]][3,2], log = T)
  gof_prop_white    <- dbeta(sim_results[[59]][iter,4], shape1=targets[[59]][4,1], shape2=targets[[59]][4,2], log = T)
  gof_prop_black    <- dbeta(sim_results[[59]][iter,5], shape1=targets[[59]][5,1], shape2=targets[[59]][5,2], log = T)
  gof_prop_hispanic <- dbeta(sim_results[[59]][iter,6], shape1=targets[[59]][6,1], shape2=targets[[59]][6,2], log = T)
  gof_prop_other    <- dbeta(sim_results[[59]][iter,7], shape1=targets[[59]][7,1], shape2=targets[[59]][7,2], log = T)
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
