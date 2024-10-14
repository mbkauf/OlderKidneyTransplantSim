## Load packages
library(doParallel)
library(doRNG)
library(foreach)
library(dtplyr)
library(dplyr)
library(doSNOW)
library(data.table)
library(survminer)
library(ggpubr)
library(gridExtra)

## File paths
in_path  <- "//tsclient/Documents/Simulation Model/Model Outputs/Recalibration/"
samp_path  <- "//tsclient/Documents/Simulation Model/Model Inputs/"
out_path <- "//tsclient/Documents/Simulation Model/Results/"
serv_path <- "C:/Users/mbkauf/Documents/Calibration/Recalibration/"

## Load model coefficients to test
m_ddtx_coef      <- fread(paste0(serv_path, "ddtx_coef.csv"))
m_ldtx_coef      <- fread(paste0(serv_path, "ldtx_coef.csv"))
m_mort_coef      <- fread(paste0(serv_path, "mort_coef.csv"))
m_remove_coef    <- fread(paste0(serv_path, "remove_coef.csv"))
m_mlogit_coef    <- fread(paste0(serv_path, "mlogit_coef.csv"))
m_gs_mort_coef   <- fread(paste0(serv_path, "gs_mort_coef.csv"))
m_gl_coef        <- fread(paste0(serv_path, "gl_coef.csv"))
m_dial_mort_coef <- fread(paste0(serv_path, "dial_mort_coef.csv"))
m_gl_mort_coef   <- fread(paste0(serv_path, "gl_mort_coef.csv"))

## Load GOF matrix
df_gof <- fread(paste0(in_path, "Stage2/part3/df_gof.csv"))

## Load data set
##### Data set #####
sample_df <- as.matrix(read.csv(file = paste0(samp_path, "sample_df_2.csv"),
                                sep=","))[,-1]

## Load model functions
source("//tsclient/Documents/R/OlderKidneyTransplantSim/R/02_model_functions.R", echo = TRUE)

## Get top 250 parameter sets
set.seed(126748)

m_gof_test <- df_gof %>%
  filter((rank_surv_post <= 1000 & rank_shape_post <= 1000) |
           (rank_surv_post <= 1000 & is.na(gof_shape_post))) %>%
  dplyr::select(id, rank_surv_post, rank_shape_post) %>%
  arrange(rank_surv_post)
v_id <- as.numeric(unlist(m_gof_test[,1]))
n <- length(v_id)

# m_gof_250 <- as.matrix(subset(df_gof, rank_surv_post <= n, select = c(id, rank_shape_post)))
# m_gof_250 <- m_gof_250[order(m_gof_250[,2], decreasing = FALSE),]
# v_id <- m_gof_250[,1]

# Keep best fitting parameter sets
m_ddtx_samp      <- m_ddtx_coef %>% filter(id %in% v_id) %>% arrange(factor(id, levels = v_id))
m_ldtx_samp      <- m_ldtx_coef %>% filter(id %in% v_id) %>% arrange(factor(id, levels = v_id))
m_mort_samp      <- m_mort_coef %>% filter(id %in% v_id) %>% arrange(factor(id, levels = v_id))
m_remove_samp    <- m_remove_coef %>% filter(id %in% v_id) %>% arrange(factor(id, levels = v_id))
m_mlogit_samp    <- m_mlogit_coef %>% filter(id %in% v_id) %>% arrange(factor(id, levels = v_id))
m_gs_mort_samp   <- m_gs_mort_coef %>% filter(id %in% v_id) %>% arrange(factor(id, levels = v_id))
m_gl_samp        <- m_gl_coef %>% filter(id %in% v_id) %>% arrange(factor(id, levels = v_id))
m_dial_mort_samp <- m_dial_mort_coef %>% filter(id %in% v_id) %>% arrange(factor(id, levels = v_id))
m_gl_mort_samp   <- m_gl_mort_coef %>% filter(id %in% v_id) %>% arrange(factor(id, levels = v_id))

m_ddtx_coef_samp      <- m_ddtx_samp   %>% dplyr::select(V2:V27) %>% as.matrix()
m_ldtx_coef_samp      <- m_ldtx_samp   %>% dplyr::select(V2:V27) %>% as.matrix()
m_mort_coef_samp      <- m_mort_samp   %>% dplyr::select(V2:V27) %>% as.matrix()
m_remove_coef_samp    <- m_remove_samp %>% dplyr::select(V2:V27) %>% as.matrix()
m_mlogit_coef_samp    <- m_mlogit_samp %>% dplyr::select(V1:V54) %>% as.matrix()
m_gs_mort_coef_samp   <- m_gs_mort_samp %>% dplyr::select(V1:V18) %>% as.matrix()
m_gl_coef_samp        <- m_gl_samp %>% dplyr::select(V1:V19) %>% as.matrix()
m_dial_mort_coef_samp <- m_dial_mort_samp %>% dplyr::select(V1:V20) %>% as.matrix()
m_gl_mort_coef_samp   <- m_gl_mort_samp %>% dplyr::select(V1:V20) %>% as.matrix()

m_ddtx_shape_samp      <- m_ddtx_samp      %>% dplyr::select(V28) %>% as.matrix()
m_ldtx_shape_samp      <- m_ldtx_samp      %>% dplyr::select(V28) %>% as.matrix()
m_mort_shape_samp      <- m_mort_samp      %>% dplyr::select(V28) %>% as.matrix()
m_remove_shape_samp    <- m_remove_samp    %>% dplyr::select(V28) %>% as.matrix()
m_gs_mort_shape_samp    <- m_gs_mort_samp   %>% dplyr::select(V19) %>% as.matrix()
m_gl_shape_samp         <- m_gl_samp        %>% dplyr::select(V20) %>% as.matrix()
m_dial_mort_shape_samp  <- m_dial_mort_samp %>% dplyr::select(V21) %>% as.matrix()

## Set up to run in parallel
comb <- function(...) {
  mapply('rbind', ..., SIMPLIFY=FALSE)
}

parallel::detectCores()
n.cores <- 75

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

out <- foreach(i = 1:n, .combine='rbind', .multicombine=TRUE, .options.snow=opts,
               .packages = c("MASS", "copula", "psych", "dtplyr",
                             "survival", "flexsurv", "dplyr")) %dopar% {

                               ## Run Waitlist Simulation
                               list_sim_df <- list.simulation(
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
                                 gl30_coef = m_mlogit_coef_samp[i, 1:18],
                                 dgf_coef = m_mlogit_coef_samp[i, 19:36],
                                 mort30_coef = m_mlogit_coef_samp[i, 37:54],
                                 gl_coef = m_gl_coef_samp[i, ],
                                 gl_shape = m_gl_shape_samp[i,1],
                                 gs_mort_coef = m_gs_mort_coef_samp[i, ],
                                 gs_mort_shape = m_gs_mort_shape_samp[i,1],
                                 dial_mort_coef = m_dial_mort_coef_samp[i, ],
                                 dial_mort_shape = m_dial_mort_shape_samp[i,1],
                                 gl_mort_coef = m_gl_mort_coef_samp[i, ]
                               )

                               ## Add id and rank
                               sim_df$id <- i
                               sim_df

                             }
parallel::stopCluster(cl = my.cluster)

## Load observed data
obs_df  <- fread(file = paste0(samp_path, "simulation_list3.csv"))

obs_df <- lazy_dt(obs_df) %>%
  mutate(ddtx = ifelse(event_cd==0, 1, 0)) %>%
  mutate(ldtx = ifelse(event_cd==1, 1, 0)) %>%
  mutate(mort = ifelse(event_cd==2, 1, 0)) %>%
  mutate(remove = ifelse(event_cd==3, 1, 0)) %>%
  mutate(sim = rep(0, nrow(obs_df))) %>%
  mutate(sim = factor(sim, levels = c(0, 1),
                      labels = c("Observed", "Simulated"))) %>%
  mutate(rec_age_at_tx_c = can_age_at_listing_c + (months_to_event/12)) %>%
  mutate(years_dial = baseline_yrs_dial + (months_to_event/12)) %>%
  mutate(tx_year_c = list_year_c + (months_to_event/12)) %>%
  mutate(race = ifelse(black==0 & hispanic==0 & other_race==0, 0,
                       ifelse(black==1 & hispanic==0 & other_race==0, 1,
                              ifelse(black==0 & hispanic==1 & other_race==0, 2, 3)))) %>%
  mutate(sim = rep(0, nrow(obs_df))) %>%
  mutate(race = factor(race, levels = c(0, 1, 2, 3),
                       labels = c("White", "Black",
                                  "Hispanic", "Other"))) %>%
  mutate(id = rep(n + 1, nrow(obs_df))) %>%
  mutate(ddtx = ifelse(is.na(ddtx), 0, ddtx),
         ldtx = ifelse(is.na(ldtx), 0, ldtx),
         mort = ifelse(is.na(mort), 0, mort),
         remove = ifelse(is.na(remove), 0, remove)) %>%
  dplyr::select(c(ddtx, ldtx, mort, remove, months_to_event, race, id, diab_stat,
                  graft_loss, months_to_gl, gs_death, months_to_gs_mort, gl_death,
                  months_to_gl_mort, tx_outcome, start_time)) %>%
  as_tibble()

sim_df <- lazy_dt(out) %>%
  rename(months_to_gs_mort = gs_death_time,
         months_to_gl = gl_time,
         months_to_gl_mort = gl_death_time
  ) %>%
  mutate(start_time = ifelse(can_age_at_listing_c >=65, 0,
                             abs(can_age_at_listing_c)-0.5)) %>%
  mutate(race = factor(race, levels = c(0, 1, 2, 3),
                       labels = c("White", "Black",
                                  "Hispanic", "Other"))) %>%
  dplyr::select(c(ddtx, ldtx, mort, remove, months_to_event, race, id, diab_stat,
                  graft_loss, months_to_gl, gs_death, months_to_gs_mort, gl_death,
                  months_to_gl_mort, tx_outcome, start_time)) %>%
  as_tibble()

all_df <- rbind(obs_df, sim_df)
all_df <- lazy_dt(all_df) %>%
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
                                    120, months_to_gl_mort)) %>%
  as_tibble()

## Graph functions
g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)}

surv_plot <- function(fit, title_txt, n_top) {
  tmp_plot <- ggsurvplot(fit,
                         censor = FALSE,
                         conf.int = FALSE,
                         title = title_txt,
                         short.panel.labs = TRUE,
                         legend.title = "",
                         ylim = c(0, 1),
                         xlab = "Months")
  tmp_plot$plot <- tmp_plot$plot +
    theme(plot.title = element_text(size = 12)) +
    theme(legend.text = element_text(size = 6)) +
    theme(axis.text.x = element_text(size = 8)) +
    theme(axis.text.y = element_text(size = 8)) +
    theme(axis.title.x = element_text(size = 8)) +
    theme(axis.title.y = element_text(size = 8)) +
    guides(fill = guide_legend(nrow = 1)) +
    theme(legend.position='none')+
    theme(strip.text.x = element_text(size = 8)) +
    scale_color_manual(values = c(rep("gray", n_top), "red")) +
    guides(color = guide_none())

  return(tmp_plot)
}

## Top 250
n.top <- n

all_df_list <- lazy_dt(all_df) %>%
  filter(id<=n.top | id==n + 1) %>%
  as_tibble()

all_df_post <- lazy_dt(all_df) %>%
  filter(ddtx==1) %>%
  filter(id<=n.top | id==n + 1) %>%
  as_tibble()

top_group <- c(n.top, n + 1)
col_vector <- c(rep("gray", n.top), "red")
names(col_vector) <- seq(1, n + 1)

# fit_ddtx    <- survfit(Surv(time = start_time, time2 = months_to_event, ddtx) ~ id,
#                        data = all_df_list)  # Deceased Donor Tx
#
# fit_ldtx    <- survfit(Surv(time = start_time, time2 = months_to_event, ldtx) ~ id,
#                        data = all_df_list)  # Living Donor Tx
#
# fit_mort    <- survfit(Surv(time = start_time, time2 = months_to_event, mort) ~ id,
#                        data = all_df_list)  # Waitlist Mortality
#
# fit_remove  <- survfit(Surv(time = start_time, time2 = months_to_event, remove) ~ id,
#                        data = all_df_list)  # Other Removal
#
# ddtx_all <- surv_plot(fit = fit_ddtx, title_txt = "Deceased Donor Tx", n_top = n.top)
# ldtx_all <- surv_plot(fit = fit_ldtx, title_txt = "Living Donor Tx", n_top = n.top)
# mort_all <- surv_plot(fit = fit_mort, title_txt = "Waitlist Mortality", n_top = n.top)
# remove_all <- surv_plot(fit = fit_remove, title_txt = "Other Waitlist Removal", n_top = n.top)
#
# remove_legend <- ggsurvplot(fit_remove,
#                             censor = FALSE,
#                             conf.int = FALSE,
#                             title = "Other Waitlist Removal",
#                             short.panel.labs = TRUE,
#                             legend.title = "",
#                             ylim = c(0, 1),
#                             xlab = "Months")
# remove_legend$plot <- remove_legend$plot +
#   theme(plot.title = element_text(size = 12)) +
#   theme(legend.text = element_text(size = 10)) +
#   theme(axis.text.x = element_text(size = 8)) +
#   theme(axis.text.y = element_text(size = 8)) +
#   theme(axis.title.x = element_text(size = 8)) +
#   theme(axis.title.y = element_text(size = 8)) +
#   theme(legend.position='right') +
#   theme(strip.text.x = element_text(size = 8)) +
#   scale_color_manual(values = col_vector, breaks = top_group,
#                      labels = c("Simulated", "Observed")) +
#   theme(legend.position=c(.5,.5)) +
#   guides(color = guide_legend("Type"))
#
#
# plots_list_250 <- grid.arrange(arrangeGrob(ddtx_all$plot, ldtx_all$plot,
#                                       mort_all$plot, remove_all$plot +
#                                         theme(legend.position="none"),
#                                       nrow = 2),
#                           nrow=2,heights=c(10, 1))
# ggsave(filename = paste0(out_path, "recalibration_list_", n.top, ".pdf"),
#        plot = plots_all , width = 7, height = 4)


fit_gl <- survfit(Surv(time = months_to_gl, graft_loss) ~ id,
                  data = all_df_post)  # Graft loss

fit_gs_death <- survfit(Surv(time = months_to_gs_mort, gs_death) ~ id,
                        data = all_df_post)  # Death w/ function

fit_gl_death <- survfit(Surv(time = months_to_gl_mort, gl_death) ~ id,
                        data = all_df_post)  # Death after graft loss

gl_all <- surv_plot(fit = fit_gl, title_txt = "Graft Loss", n_top = n.top)
gs_death_all <- surv_plot(fit = fit_gs_death, title_txt = "Death w/ Function", n_top = n.top)
gl_death_all <- surv_plot(fit = fit_gl_death, title_txt = "Death after Graft Loss", n_top = n.top)



plots_post_250 <- grid.arrange(arrangeGrob(gl_all$plot, gs_death_all$plot,
                                       gl_death_all$plot +
                                         theme(legend.position="none"),
                                       nrow = 2),
                           nrow=2,heights=c(10, 1))

# ggsave(filename = paste0(out_path, "recalibration_post_", n.top, ".pdf"),
#        plot = plots_post, width = 7, height = 4)


## Top 100
n.top <- 100

all_df_list <- lazy_dt(all_df) %>%
  filter(id<=n.top | id==n + 1) %>%
  as_tibble()

all_df_post <- lazy_dt(all_df) %>%
  filter(ddtx==1) %>%
  filter(id<=n.top | id==n + 1) %>%
  as_tibble()

top_group <- c(n.top, n + 1)
col_vector <- c(rep("gray", n.top), "red")
# names(col_vector) <- seq(1, n + 1)

# fit_ddtx    <- survfit(Surv(time = start_time, time2 = months_to_event, ddtx) ~ id,
#                        data = all_df_list)  # Deceased Donor Tx
#
# fit_ldtx    <- survfit(Surv(time = start_time, time2 = months_to_event, ldtx) ~ id,
#                        data = all_df_list)  # Living Donor Tx
#
# fit_mort    <- survfit(Surv(time = start_time, time2 = months_to_event, mort) ~ id,
#                        data = all_df_list)  # Waitlist Mortality
#
# fit_remove  <- survfit(Surv(time = start_time, time2 = months_to_event, remove) ~ id,
#                        data = all_df_list)  # Other Removal
#
# ddtx_all <- surv_plot(fit = fit_ddtx, title_txt = "Deceased Donor Tx", n_top = n.top)
# ldtx_all <- surv_plot(fit = fit_ldtx, title_txt = "Living Donor Tx", n_top = n.top)
# mort_all <- surv_plot(fit = fit_mort, title_txt = "Waitlist Mortality", n_top = n.top)
# remove_all <- surv_plot(fit = fit_remove, title_txt = "Other Waitlist Removal", n_top = n.top)
#
# remove_legend <- ggsurvplot(fit_remove,
#                             censor = FALSE,
#                             conf.int = FALSE,
#                             title = "Other Waitlist Removal",
#                             short.panel.labs = TRUE,
#                             legend.title = "",
#                             ylim = c(0, 1),
#                             xlab = "Months")
# remove_legend$plot <- remove_legend$plot +
#   theme(plot.title = element_text(size = 12)) +
#   theme(legend.text = element_text(size = 10)) +
#   theme(axis.text.x = element_text(size = 8)) +
#   theme(axis.text.y = element_text(size = 8)) +
#   theme(axis.title.x = element_text(size = 8)) +
#   theme(axis.title.y = element_text(size = 8)) +
#   theme(legend.position='right') +
#   theme(strip.text.x = element_text(size = 8)) +
#   scale_color_manual(values = col_vector, breaks = top_group,
#                      labels = c("Simulated", "Observed")) +
#   theme(legend.position=c(.5,.5)) +
#   guides(color = guide_legend("Type"))
#
#
# plots_list_100 <- grid.arrange(arrangeGrob(ddtx_all$plot, ldtx_all$plot,
#                                       mort_all$plot, remove_all$plot +
#                                         theme(legend.position="none"),
#                                       nrow = 2),
#                           nrow=2,heights=c(10, 1))
# ggsave(filename = paste0(out_path, "recalibration_list_", n.top, ".pdf"),
#        plot = plots_all , width = 7, height = 4)


fit_gl <- survfit(Surv(time = months_to_gl, graft_loss) ~ id,
                  data = all_df_post)  # Graft loss

fit_gs_death <- survfit(Surv(time = months_to_gs_mort, gs_death) ~ id,
                        data = all_df_post)  # Death w/ function

fit_gl_death <- survfit(Surv(time = months_to_gl_mort, gl_death) ~ id,
                        data = all_df_post)  # Death after graft loss

gl_all <- surv_plot(fit = fit_gl, title_txt = "Graft Loss", n_top = n.top)
gs_death_all <- surv_plot(fit = fit_gs_death, title_txt = "Death w/ Function", n_top = n.top)
gl_death_all <- surv_plot(fit = fit_gl_death, title_txt = "Death after Graft Loss", n_top = n.top)



plots_post_100 <- grid.arrange(arrangeGrob(gl_all$plot, gs_death_all$plot,
                                       gl_death_all$plot +
                                         theme(legend.position="none"),
                                       nrow = 2),
                           nrow=2,heights=c(10, 1))

# ggsave(filename = paste0(out_path, "recalibration_post_", n.top, ".pdf"),
#        plot = plots_post , width = 7, height = 4)

## Top 50
n.top <- 50

all_df_list <- lazy_dt(all_df) %>%
  filter(id<=n.top | id==n + 1) %>%
  as_tibble()

all_df_post <- lazy_dt(all_df) %>%
  filter(ddtx==1) %>%
  filter(id<=n.top | id==n + 1) %>%
  as_tibble()

top_group <- c(n.top, n + 1)
col_vector <- c(rep("gray", n.top), "red")
# names(col_vector) <- seq(1, n + 1)

# fit_ddtx    <- survfit(Surv(time = start_time, time2 = months_to_event, ddtx) ~ id,
#                        data = all_df_list)  # Deceased Donor Tx
#
# fit_ldtx    <- survfit(Surv(time = start_time, time2 = months_to_event, ldtx) ~ id,
#                        data = all_df_list)  # Living Donor Tx
#
# fit_mort    <- survfit(Surv(time = start_time, time2 = months_to_event, mort) ~ id,
#                        data = all_df_list)  # Waitlist Mortality
#
# fit_remove  <- survfit(Surv(time = start_time, time2 = months_to_event, remove) ~ id,
#                        data = all_df_list)  # Other Removal
#
# ddtx_all <- surv_plot(fit = fit_ddtx, title_txt = "Deceased Donor Tx", n_top = n.top)
# ldtx_all <- surv_plot(fit = fit_ldtx, title_txt = "Living Donor Tx", n_top = n.top)
# mort_all <- surv_plot(fit = fit_mort, title_txt = "Waitlist Mortality", n_top = n.top)
# remove_all <- surv_plot(fit = fit_remove, title_txt = "Other Waitlist Removal", n_top = n.top)
#
# remove_legend <- ggsurvplot(fit_remove,
#                             censor = FALSE,
#                             conf.int = FALSE,
#                             title = "Other Waitlist Removal",
#                             short.panel.labs = TRUE,
#                             legend.title = "",
#                             ylim = c(0, 1),
#                             xlab = "Months")
# remove_legend$plot <- remove_legend$plot +
#   theme(plot.title = element_text(size = 12)) +
#   theme(legend.text = element_text(size = 10)) +
#   theme(axis.text.x = element_text(size = 8)) +
#   theme(axis.text.y = element_text(size = 8)) +
#   theme(axis.title.x = element_text(size = 8)) +
#   theme(axis.title.y = element_text(size = 8)) +
#   theme(legend.position='right') +
#   theme(strip.text.x = element_text(size = 8)) +
#   scale_color_manual(values = col_vector, breaks = top_group,
#                      labels = c("Simulated", "Observed")) +
#   theme(legend.position=c(.5,.5)) +
#   guides(color = guide_legend("Type"))
#
#
# plots_list_50 <- grid.arrange(arrangeGrob(ddtx_all$plot, ldtx_all$plot,
#                                       mort_all$plot, remove_all$plot +
#                                         theme(legend.position="none"),
#                                       nrow = 2),
#                           nrow=2,heights=c(10, 1))
# ggsave(filename = paste0(out_path, "recalibration_list_", n.top, ".pdf"),
#        plot = plots_all , width = 7, height = 4)


fit_gl <- survfit(Surv(time = months_to_gl, graft_loss) ~ id,
                  data = all_df_post)  # Graft loss

fit_gs_death <- survfit(Surv(time = months_to_gs_mort, gs_death) ~ id,
                        data = all_df_post)  # Death w/ function

fit_gl_death <- survfit(Surv(time = months_to_gl_mort, gl_death) ~ id,
                        data = all_df_post)  # Death after graft loss

gl_all <- surv_plot(fit = fit_gl, title_txt = "Graft Loss", n_top = n.top)
gs_death_all <- surv_plot(fit = fit_gs_death, title_txt = "Death w/ Function", n_top = n.top)
gl_death_all <- surv_plot(fit = fit_gl_death, title_txt = "Death after Graft Loss", n_top = n.top)



plots_post_50 <- grid.arrange(arrangeGrob(gl_all$plot, gs_death_all$plot,
                                       gl_death_all$plot +
                                         theme(legend.position="none"),
                                       nrow = 2),
                           nrow=2,heights=c(10, 1))

# ggsave(filename = paste0(out_path, "recalibration_post_", n.top, ".pdf"),
#        plot = plots_post , width = 7, height = 4)

## Top 10
n.top <- 10

all_df_list <- lazy_dt(all_df) %>%
  filter(id<=n.top | id==n + 1) %>%
  as_tibble()

all_df_post <- lazy_dt(all_df) %>%
  filter(ddtx==1) %>%
  filter(id<=n.top | id==n + 1) %>%
  as_tibble()

top_group <- c(n.top, n + 1)
col_vector <- c(rep("gray", n.top), "red")
# names(col_vector) <- seq(1, n + 1)

# fit_ddtx    <- survfit(Surv(time = start_time, time2 = months_to_event, ddtx) ~ id,
#                        data = all_df_list)  # Deceased Donor Tx
#
# fit_ldtx    <- survfit(Surv(time = start_time, time2 = months_to_event, ldtx) ~ id,
#                        data = all_df_list)  # Living Donor Tx
#
# fit_mort    <- survfit(Surv(time = start_time, time2 = months_to_event, mort) ~ id,
#                        data = all_df_list)  # Waitlist Mortality
#
# fit_remove  <- survfit(Surv(time = start_time, time2 = months_to_event, remove) ~ id,
#                        data = all_df_list)  # Other Removal
#
# ddtx_all <- surv_plot(fit = fit_ddtx, title_txt = "Deceased Donor Tx", n_top = n.top)
# ldtx_all <- surv_plot(fit = fit_ldtx, title_txt = "Living Donor Tx", n_top = n.top)
# mort_all <- surv_plot(fit = fit_mort, title_txt = "Waitlist Mortality", n_top = n.top)
# remove_all <- surv_plot(fit = fit_remove, title_txt = "Other Waitlist Removal", n_top = n.top)
#
# remove_legend <- ggsurvplot(fit_remove,
#                             censor = FALSE,
#                             conf.int = FALSE,
#                             title = "Other Waitlist Removal",
#                             short.panel.labs = TRUE,
#                             legend.title = "",
#                             ylim = c(0, 1),
#                             xlab = "Months")
# remove_legend$plot <- remove_legend$plot +
#   theme(plot.title = element_text(size = 12)) +
#   theme(legend.text = element_text(size = 10)) +
#   theme(axis.text.x = element_text(size = 8)) +
#   theme(axis.text.y = element_text(size = 8)) +
#   theme(axis.title.x = element_text(size = 8)) +
#   theme(axis.title.y = element_text(size = 8)) +
#   theme(legend.position='right') +
#   theme(strip.text.x = element_text(size = 8)) +
#   scale_color_manual(values = col_vector, breaks = top_group,
#                      labels = c("Simulated", "Observed")) +
#   theme(legend.position=c(.5,.5)) +
#   guides(color = guide_legend("Type"))
#
#
# plots_list_10 <- grid.arrange(arrangeGrob(ddtx_all$plot, ldtx_all$plot,
#                                       mort_all$plot, remove_all$plot +
#                                         theme(legend.position="none"),
#                                       nrow = 2),
#                           nrow=2,heights=c(10, 1))
# ggsave(filename = paste0(out_path, "recalibration_list_", n.top, ".pdf"),
#        plot = plots_all , width = 7, height = 4)


fit_gl <- survfit(Surv(time = months_to_gl, graft_loss) ~ id,
                  data = all_df_post)  # Graft loss

fit_gs_death <- survfit(Surv(time = months_to_gs_mort, gs_death) ~ id,
                        data = all_df_post)  # Death w/ function

fit_gl_death <- survfit(Surv(time = months_to_gl_mort, gl_death) ~ id,
                        data = all_df_post)  # Death after graft loss

gl_all <- surv_plot(fit = fit_gl, title_txt = "Graft Loss", n_top = n.top)
gs_death_all <- surv_plot(fit = fit_gs_death, title_txt = "Death w/ Function", n_top = n.top)
gl_death_all <- surv_plot(fit = fit_gl_death, title_txt = "Death after Graft Loss", n_top = n.top)



plots_post_10 <- grid.arrange(arrangeGrob(gl_all$plot, gs_death_all$plot,
                                       gl_death_all$plot +
                                         theme(legend.position="none"),
                                       nrow = 2),
                           nrow=2,heights=c(10, 1))

# ggsave(filename = paste0(out_path, "recalibration_post_", n.top, ".pdf"),
#        plot = plots_post , width = 7, height = 4)

##
library(stringr)

df_mlogit_250 <- data.frame(unclass(table(all_df_250$id, all_df_250$tx_outcome)))
df_mlogit_100 <- data.frame(unclass(table(all_df_100$id, all_df_100$tx_outcome)))
df_mlogit_50 <- data.frame(unclass(table(all_df_50$id, all_df_50$tx_outcome)))

sample_summary <- function(data=NULL) {
  means <- sapply(data, mean)
  means <- unlist(means)
  upper <- sapply(data, quantile, prob = 0.975)
  lower <- sapply(data, quantile, prob = 0.025)
  outcome <- c("GS", "GL", "DGF", "Mort")

  data_summary <- as.data.frame(rbind(means, upper, lower))
  colnames(data_summary) <- outcome

  data_summary$type <- c("target", "upper", "lower")

  return(data_summary)
}

sample_summary_glmort <- function(data=NULL) {
  means <- sapply(data, mean)
  means <- unlist(means)
  upper <- sapply(data, quantile, prob = 0.975)
  lower <- sapply(data, quantile, prob = 0.025)
  outcome <- c("Death_at_GL")

  data_summary <- as.data.frame(rbind(means, upper, lower))
  colnames(data_summary) <- outcome

  data_summary$type <- c("target", "upper", "lower")

  return(data_summary)
}

summary_df_mlogit_long <- function(data=NULL, n) {
  t_mlogit <- table(data$id, data$tx_outcome)
  df_mlogit <- data.frame(unclass(prop.table(t_mlogit, margin = 1)))
  colnames(df_mlogit) <- c("GS", "GL", "DGF", "Mort")
  df_obs <- df_mlogit[(n+1),]
  df_sim <- df_mlogit[1:n,]
  df_sim_summary <- sample_summary(data = df_sim)
  df_sim_summary$sim <- "Simulated"
  df_obs$type <- "target"
  df_obs$sim <- "Observed"
  na_obs_upper <- c(NA, NA, NA, NA, "upper", "Observed")
  na_obs_lower <- c(NA, NA, NA, NA, "lower", "Observed")
  df_obs <- rbind(df_obs, na_obs_upper, na_obs_lower)
  df_all <- rbind(df_sim_summary, df_obs)

  df_means <- df_all %>% dplyr::filter(type=="target")
  df_upper <- df_all %>% dplyr::filter(type=="upper")
  df_lower <- df_all %>% dplyr::filter(type=="lower")

  df_means_long <- gather(df_means, "outcome", "mean", -c(sim, type))
  df_upper_long <- gather(df_upper, "outcome", "upper", -c(sim, type))
  df_lower_long <- gather(df_lower, "outcome", "lower", -c(sim, type))
  upper <- as.numeric(df_upper_long$upper)
  lower <- as.numeric(df_lower_long$lower)
  df_all_long <- cbind(df_means_long, upper, lower)
  df_all_long$mean <- as.numeric(df_all_long$mean)
  df_all_long$upper <- as.numeric(df_all_long$upper)
  df_all_long$lower <- as.numeric(df_all_long$lower)

  return(df_all_long)
}

summary_df_glmort_long <- function(data=NULL, n) {
  data <- lazy_dt(data) %>%
    mutate(death_at_gl = ifelse(gl_death==1 & months_to_gl_mort==0, 1,
                                ifelse(is.na(months_to_gl_mort), NaN, 0))) %>%
    as_tibble()

  t_glmort  <- table(data$id, data$death_at_gl)
  df_glmort <- data.frame(unclass(prop.table(t_glmort, margin = 1)))
  colnames(df_glmort) <- c("No", "Death_at_GL")
  df_obs <- df_glmort[(n+1),2]
  df_sim <- as.data.frame(matrix(df_glmort[1:n,2], ncol = 1))
  df_sim_summary <- sample_summary_glmort(data = df_sim)
  df_sim_summary$sim <- "Simulated"
  df_obs <- matrix(c(df_obs, "target"), ncol = 2)
  df_obs <- as.data.frame(df_obs)
  df_obs$sim <- "Observed"
  na_obs_upper <- c(NA, "upper", "Observed")
  na_obs_lower <- c(NA, "lower", "Observed")
  df_obs <- rbind(df_obs, na_obs_upper, na_obs_lower)
  colnames(df_obs) <- c("Death_at_GL", "type", "sim")
  df_all <- rbind(df_sim_summary, df_obs)

  df_means <- df_all %>% dplyr::filter(type=="target")
  df_upper <- df_all %>% dplyr::filter(type=="upper")
  df_lower <- df_all %>% dplyr::filter(type=="lower")

  df_means_long <- gather(df_means, "outcome", "mean", -c(sim, type))
  df_upper_long <- gather(df_upper, "outcome", "upper", -c(sim, type))
  df_lower_long <- gather(df_lower, "outcome", "lower", -c(sim, type))
  upper <- as.numeric(df_upper_long$upper)
  lower <- as.numeric(df_lower_long$lower)
  df_all_long <- cbind(df_means_long, upper, lower)
  df_all_long$mean <- as.numeric(df_all_long$mean)
  df_all_long$upper <- as.numeric(df_all_long$upper)
  df_all_long$lower <- as.numeric(df_all_long$lower)
  df_all_long[2,5:6] <- NA

  return(df_all_long)

}

all_df_ddtx_250 <-lazy_dt(all_df) %>%
  filter(ddtx==1) %>%
  filter(id<=n | id==n + 1) %>%
  as_tibble()
all_df_ddtx_100 <-lazy_dt(all_df) %>%
  filter(ddtx==1) %>%
  filter(id<=100 | id==n + 1) %>%
  as_tibble()
all_df_ddtx_50 <-lazy_dt(all_df) %>%
  filter(ddtx==1) %>%
  filter(id<=50 | id==n + 1) %>%
  as_tibble()
all_df_ddtx_10 <-lazy_dt(all_df) %>%
  filter(ddtx==1) %>%
  filter(id<=10 | id==n + 1) %>%
  as_tibble()

df_mlogit_250_all_long <- summary_df_mlogit_long(data = all_df_ddtx_250, n = n)
df_mlogit_100_all_long <- summary_df_mlogit_long(data = all_df_ddtx_100, n = 100)
df_mlogit_50_all_long  <- summary_df_mlogit_long(data = all_df_ddtx_50, n = 50)
df_mlogit_10_all_long  <- summary_df_mlogit_long(data = all_df_ddtx_10, n = 10)

df_glmort_250_all_long <- summary_df_glmort_long(data = all_df_ddtx_250, n = n)
df_glmort_100_all_long <- summary_df_glmort_long(data = all_df_ddtx_100, n = 100)
df_glmort_50_all_long  <- summary_df_glmort_long(data = all_df_ddtx_50, n = 50)
df_glmort_10_all_long  <- summary_df_glmort_long(data = all_df_ddtx_10, n = 10)

## Create plots
mlogit_plot <- function(data=NULL, n) {
  pd <- position_dodge(0.8)
  tmp_plot <- ggplot(data = data, mapping = aes(x=outcome, colour=sim, y=mean)) +
    geom_point(position = pd, size=2) +
    geom_errorbar(data = data, aes(x=outcome, ymin=lower, ymax=upper, colour=sim), linewidth=1,
                  size = 0.8, position = pd) +
    theme_bw() +
    scale_color_manual(values = c("#440154FF", "#E69F00")) +
    theme(legend.position="right",
          strip.text.x = element_text(size = 10),
          axis.text.x = element_text(size = 10),
          axis.text.y = element_text(size = 10),
          axis.title.x = element_blank()) +
    # ylim(0, 1) +
    ggtitle(paste("30-Day Outcomes: Top", n))


  return(tmp_plot)
}

glmort_plot <- function(data=NULL, n) {
  pd <- position_dodge(0.8)
  tmp_plot <- ggplot(data = data, mapping = aes(x=outcome, colour=sim, y=mean)) +
    geom_point(position = pd, size=2) +
    geom_errorbar(data = data, aes(x=outcome, ymin=lower, ymax=upper, colour=sim), linewidth=1,
                  size = 0.8, position = pd) +
    theme_bw() +
    scale_color_manual(values = c("#440154FF", "#E69F00")) +
    theme(legend.position="right",
          strip.text.x = element_text(size = 10),
          axis.text.x = element_text(size = 10),
          axis.text.y = element_text(size = 10),
          axis.title.x = element_blank()) +
    # ylim(0, 1) +
    ggtitle(paste("Death at Graft Loss: Top", n))


  return(tmp_plot)
}

mlogit_250_plot <- mlogit_plot(data = df_mlogit_250_all_long, n = n)
mlogit_100_plot <- mlogit_plot(data = df_mlogit_100_all_long, n = 100)
mlogit_50_plot  <- mlogit_plot(data = df_mlogit_50_all_long, n = 50)
mlogit_10_plot  <- mlogit_plot(data = df_mlogit_10_all_long, n = 10)

glmort_250_plot <- glmort_plot(data = df_glmort_250_all_long, n = n)
glmort_100_plot <- glmort_plot(data = df_glmort_100_all_long, n = 100)
glmort_50_plot  <- glmort_plot(data = df_glmort_50_all_long, n = 50)
glmort_10_plot  <- glmort_plot(data = df_glmort_10_all_long, n = 10)

gl_250 <- list(mlogit_250_plot, glmort_250_plot)
gl_100 <- list(mlogit_100_plot, glmort_100_plot)
gl_50  <- list(mlogit_50_plot, glmort_50_plot)
gl_10  <- list(mlogit_10_plot, glmort_10_plot)

lay <- rbind(c(1, 1, 1, 2, 2),
             c(1, 1, 1, 2, 2))

plot_250 <- gridExtra::grid.arrange(grobs=gl_250,
                                    layout_matrix = lay)
plot_100 <- gridExtra::grid.arrange(grobs=gl_100,
                                    layout_matrix = lay)
plot_50  <- gridExtra::grid.arrange(grobs=gl_50,
                                    layout_matrix = lay)
plot_10  <- gridExtra::grid.arrange(grobs=gl_10,
                                    layout_matrix = lay)

ggsave(filename = paste0(serv_path, "recalibration_prop_694.pdf"),
       plot = plot_250 , width = 5, height = 4)
ggsave(filename = paste0(serv_path, "recalibration_prop_100.pdf"),
       plot = plot_100 , width = 5, height = 4)
ggsave(filename = paste0(serv_path, "recalibration_prop_50.pdf"),
       plot = plot_50 , width = 5, height = 4)
ggsave(filename = paste0(serv_path, "recalibration_prop_10.pdf"),
       plot = plot_10 , width = 5, height = 4)

### Compare table 1
library(modi)
library(tidyverse)
library(xtable)

# Reload and clean observed data
## Load observed data
obs_df  <- fread(file = paste0(samp_path, "simulation_list3.csv"))

obs_df <- lazy_dt(obs_df) %>%
  mutate(ddtx = ifelse(event_cd==0, 1, 0)) %>%
  mutate(ldtx = ifelse(event_cd==1, 1, 0)) %>%
  mutate(mort = ifelse(event_cd==2, 1, 0)) %>%
  mutate(remove = ifelse(event_cd==3, 1, 0)) %>%
  mutate(sim = rep(0, nrow(obs_df))) %>%
  mutate(sim = factor(sim, levels = c(0, 1),
                      labels = c("Observed", "Simulated"))) %>%
  mutate(rec_age_at_tx_c = can_age_at_listing_c + (months_to_event/12)) %>%
  mutate(years_dial = baseline_yrs_dial + (months_to_event/12)) %>%
  mutate(tx_year_c = list_year_c + (months_to_event/12)) %>%
  mutate(race = ifelse(black==0 & hispanic==0 & other_race==0, 0,
                       ifelse(black==1 & hispanic==0 & other_race==0, 1,
                              ifelse(black==0 & hispanic==1 & other_race==0, 2, 3)))) %>%
  mutate(sim = rep(0, nrow(obs_df))) %>%
  mutate(race = factor(race, levels = c(0, 1, 2, 3),
                       labels = c("White", "Black",
                                  "Hispanic", "Other"))) %>%
  mutate(id = rep(n + 1, nrow(obs_df))) %>%
  mutate(ddtx = ifelse(is.na(ddtx), 0, ddtx),
         ldtx = ifelse(is.na(ldtx), 0, ldtx),
         mort = ifelse(is.na(mort), 0, mort),
         remove = ifelse(is.na(remove), 0, remove)) %>%
  as_tibble()

# Functions
table1_clean <- function(data){
  data <- lazy_dt(data) %>%
    mutate(can_blood_a = ifelse(can_blood_ab==0 & can_blood_b==0 & can_blood_o==0,
                                1, 0)) %>%
    mutate(rec_age_at_tx_c = can_age_at_listing_c + (months_to_event/12)) %>%
    mutate(rec_age_at_tx = rec_age_at_tx_c + 65) %>%
    mutate(years_dial = baseline_yrs_dial + (months_to_event/12)) %>%
    mutate(tx_year_c = list_year_c + (months_to_event/12)) %>%
    mutate(optn_reg1 = ifelse(optn_reg2==0 & optn_reg3==0 & optn_reg4==0 &
                                optn_reg5==0 & optn_reg6==0 & optn_reg7==0 &
                                optn_reg8==0 & optn_reg9==0 & optn_reg10==0 &
                                optn_reg11==0, 1, 0)) %>%
    mutate(white = ifelse(black==0 & hispanic==0 & other_race==0, 1, 0)) %>%
    filter(ddtx==1) %>%
    as_tibble()

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

  sim_table1 <- matrix(data = NA, nrow = n, ncol = 27)

  sim_table1[,1]  <- aggregate(data$rec_age_at_tx, list(data$id), FUN=mean)[,2]
  sim_table1[,2]  <- aggregate(data$male, list(data$id), FUN=mean)[,2]
  sim_table1[,3]  <- aggregate(data$white, list(data$id), FUN=mean)[,2]
  sim_table1[,4]  <- aggregate(data$black, list(data$id), FUN=mean)[,2]
  sim_table1[,5]  <- aggregate(data$hispanic, list(data$id), FUN=mean)[,2]
  sim_table1[,6]  <- aggregate(data$other_race, list(data$id), FUN=mean)[,2]
  sim_table1[,7]  <- aggregate(data$can_blood_a, list(data$id), FUN=mean)[,2]
  sim_table1[,8]  <- aggregate(data$can_blood_ab, list(data$id), FUN=mean)[,2]
  sim_table1[,9]  <- aggregate(data$can_blood_b, list(data$id), FUN=mean)[,2]
  sim_table1[,10] <- aggregate(data$can_blood_o, list(data$id), FUN=mean)[,2]
  sim_table1[,11] <- aggregate(data$years_dial, list(data$id), FUN=mean)[,2]
  sim_table1[,12] <- aggregate(data$diab_stat, list(data$id), FUN=mean)[,2]
  sim_table1[,13] <- aggregate(data$copd, list(data$id), FUN=mean)[,2]
  sim_table1[,14] <- aggregate(data$pvd, list(data$id), FUN=mean)[,2]
  sim_table1[,15] <- aggregate(data$ang_cad, list(data$id), FUN=mean)[,2]
  sim_table1[,16] <- aggregate(data$canhx_cpra, list(data$id), FUN=mean)[,2]
  sim_table1[,17] <- aggregate(data$optn_reg1, list(data$id), FUN=mean)[,2]
  sim_table1[,18] <- aggregate(data$optn_reg2, list(data$id), FUN=mean)[,2]
  sim_table1[,19] <- aggregate(data$optn_reg3, list(data$id), FUN=mean)[,2]
  sim_table1[,20] <- aggregate(data$optn_reg4, list(data$id), FUN=mean)[,2]
  sim_table1[,21] <- aggregate(data$optn_reg5, list(data$id), FUN=mean)[,2]
  sim_table1[,22] <- aggregate(data$optn_reg6, list(data$id), FUN=mean)[,2]
  sim_table1[,23] <- aggregate(data$optn_reg7, list(data$id), FUN=mean)[,2]
  sim_table1[,24] <- aggregate(data$optn_reg8, list(data$id), FUN=mean)[,2]
  sim_table1[,25] <- aggregate(data$optn_reg9, list(data$id), FUN=mean)[,2]
  sim_table1[,26] <- aggregate(data$optn_reg10, list(data$id), FUN=mean)[,2]
  sim_table1[,27] <- aggregate(data$optn_reg11, list(data$id), FUN=mean)[,2]

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


tx_df <- table1_clean(data = out)
# table1_w_surv <- weighted.table1(tx_df_surv, score = df_random_surv$w_surv)
table1 <- nonweighted.table1(tx_df)

# Observed Table 1
tx_df_obs <- table1_clean(data = obs_df)[1,]

# Difference Table 1
table1_all <- rbind(tx_df_obs, table1)
row.names(table1_all) <- c("Observed",
                           "Simulated Mean", "Simulated 97.5",
                           "Simulated 2.5")
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
write.csv(table1_all, paste0(serv_path, "table1.csv"))

