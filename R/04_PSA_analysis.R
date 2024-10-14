library(data.table)
library(dplyr)
library(dtplyr)
library(stringr)
library(ggplot2)
library(dampack)
library(ggpubr)
library(ggsci)

df_base      <- matrix(nrow = (64*100), ncol = 92)
df_policy_1  <- matrix(nrow = (64*100), ncol = 92)
df_policy_2  <- matrix(nrow = (64*100), ncol = 92)
df_policy_3  <- matrix(nrow = (64*100), ncol = 92)
df_policy_4  <- matrix(nrow = (64*100), ncol = 92)
df_policy_5  <- matrix(nrow = (64*100), ncol = 92)
df_policy_6  <- matrix(nrow = (64*100), ncol = 92)
df_policy_7  <- matrix(nrow = (64*100), ncol = 92)
df_policy_8  <- matrix(nrow = (64*100), ncol = 92)
df_policy_9  <- matrix(nrow = (64*100), ncol = 92)
df_policy_10 <- matrix(nrow = (64*100), ncol = 92)
df_policy_11 <- matrix(nrow = (64*100), ncol = 92)
df_policy_12 <- matrix(nrow = (64*100), ncol = 92)
df_policy_13 <- matrix(nrow = (64*100), ncol = 92)
df_policy_14 <- matrix(nrow = (64*100), ncol = 92)
df_policy_15 <- matrix(nrow = (64*100), ncol = 92)
df_policy_16 <- matrix(nrow = (64*100), ncol = 92)
df_policy_17 <- matrix(nrow = (64*100), ncol = 92)
df_policy_18 <- matrix(nrow = (64*100), ncol = 92)
df_policy_19 <- matrix(nrow = (64*100), ncol = 92)
df_policy_20 <- matrix(nrow = (64*100), ncol = 92)

for (i in 1:100) {
  tmp_base_df            <- fread(paste0("results/PSA/", i,"/base_df.csv"))
  tmp_policy_1_cq_counts <- fread(paste0("results/PSA/", i,"/policy_1_cq_counts.csv"))
  tmp_policy_2_cq_counts <- fread(paste0("results/PSA/", i,"/policy_2_cq_counts.csv"))
  tmp_policy_3_cq_counts <- fread(paste0("results/PSA/", i,"/policy_3_cq_counts.csv"))
  tmp_policy_4_cq_counts <- fread(paste0("results/PSA/", i,"/policy_4_cq_counts.csv"))
  tmp_policy_5_cq_counts <- fread(paste0("results/PSA/", i,"/policy_5_cq_counts.csv"))
  tmp_policy_6_cq_counts <- fread(paste0("results/PSA/", i,"/policy_6_cq_counts.csv"))
  tmp_policy_7_cq_counts <- fread(paste0("results/PSA/", i,"/policy_7_cq_counts.csv"))
  tmp_policy_8_cq_counts <- fread(paste0("results/PSA/", i,"/policy_8_cq_counts.csv"))
  tmp_policy_9_cq_counts <- fread(paste0("results/PSA/", i,"/policy_9_cq_counts.csv"))
  tmp_policy_10_cq_counts <- fread(paste0("results/PSA/", i,"/policy_10_cq_counts.csv"))
  tmp_policy_11_cq_counts <- fread(paste0("results/PSA/", i,"/policy_11_cq_counts.csv"))
  tmp_policy_12_cq_counts <- fread(paste0("results/PSA/", i,"/policy_12_cq_counts.csv"))
  tmp_policy_13_cq_counts <- fread(paste0("results/PSA/", i,"/policy_13_cq_counts.csv"))
  tmp_policy_14_cq_counts <- fread(paste0("results/PSA/", i,"/policy_14_cq_counts.csv"))
  tmp_policy_15_cq_counts <- fread(paste0("results/PSA/", i,"/policy_15_cq_counts.csv"))
  tmp_policy_16_cq_counts <- fread(paste0("results/PSA/", i,"/policy_16_cq_counts.csv"))
  tmp_policy_17_cq_counts <- fread(paste0("results/PSA/", i,"/policy_17_cq_counts.csv"))
  tmp_policy_18_cq_counts <- fread(paste0("results/PSA/", i,"/policy_18_cq_counts.csv"))
  tmp_policy_19_cq_counts <- fread(paste0("results/PSA/", i,"/policy_19_cq_counts.csv"))
  tmp_policy_20_cq_counts <- fread(paste0("results/PSA/", i,"/policy_20_cq_counts.csv"))

  df_base[(1 + (i - 1)*64):(i*64), 1:91]      <- as.matrix(tmp_base_df)
  df_policy_1[(1 + (i - 1)*64):(i*64), 1:91]  <- as.matrix(tmp_policy_1_cq_counts)
  df_policy_2[(1 + (i - 1)*64):(i*64), 1:91]  <- as.matrix(tmp_policy_2_cq_counts)
  df_policy_3[(1 + (i - 1)*64):(i*64), 1:91]  <- as.matrix(tmp_policy_3_cq_counts)
  df_policy_4[(1 + (i - 1)*64):(i*64), 1:91]  <- as.matrix(tmp_policy_4_cq_counts)
  df_policy_5[(1 + (i - 1)*64):(i*64), 1:91]  <- as.matrix(tmp_policy_5_cq_counts)
  df_policy_6[(1 + (i - 1)*64):(i*64), 1:91]  <- as.matrix(tmp_policy_6_cq_counts)
  df_policy_7[(1 + (i - 1)*64):(i*64), 1:91]  <- as.matrix(tmp_policy_7_cq_counts)
  df_policy_8[(1 + (i - 1)*64):(i*64), 1:91]  <- as.matrix(tmp_policy_8_cq_counts)
  df_policy_9[(1 + (i - 1)*64):(i*64), 1:91]  <- as.matrix(tmp_policy_9_cq_counts)
  df_policy_10[(1 + (i - 1)*64):(i*64), 1:91] <- as.matrix(tmp_policy_10_cq_counts)
  df_policy_11[(1 + (i - 1)*64):(i*64), 1:91] <- as.matrix(tmp_policy_11_cq_counts)
  df_policy_12[(1 + (i - 1)*64):(i*64), 1:91] <- as.matrix(tmp_policy_12_cq_counts)
  df_policy_13[(1 + (i - 1)*64):(i*64), 1:91] <- as.matrix(tmp_policy_13_cq_counts)
  df_policy_14[(1 + (i - 1)*64):(i*64), 1:91] <- as.matrix(tmp_policy_14_cq_counts)
  df_policy_15[(1 + (i - 1)*64):(i*64), 1:91] <- as.matrix(tmp_policy_15_cq_counts)
  df_policy_16[(1 + (i - 1)*64):(i*64), 1:91] <- as.matrix(tmp_policy_16_cq_counts)
  df_policy_17[(1 + (i - 1)*64):(i*64), 1:91] <- as.matrix(tmp_policy_17_cq_counts)
  df_policy_18[(1 + (i - 1)*64):(i*64), 1:91] <- as.matrix(tmp_policy_18_cq_counts)
  df_policy_19[(1 + (i - 1)*64):(i*64), 1:91] <- as.matrix(tmp_policy_19_cq_counts)
  df_policy_20[(1 + (i - 1)*64):(i*64), 1:91] <- as.matrix(tmp_policy_20_cq_counts)

  df_base[(1 + (i - 1)*64):(i*64), 92] <- rep(i, 64)
  df_policy_1[(1 + (i - 1)*64):(i*64), 92] <- rep(i, 64)
  df_policy_2[(1 + (i - 1)*64):(i*64), 92] <- rep(i, 64)
  df_policy_3[(1 + (i - 1)*64):(i*64), 92] <- rep(i, 64)
  df_policy_4[(1 + (i - 1)*64):(i*64), 92] <- rep(i, 64)
  df_policy_5[(1 + (i - 1)*64):(i*64), 92] <- rep(i, 64)
  df_policy_6[(1 + (i - 1)*64):(i*64), 92] <- rep(i, 64)
  df_policy_7[(1 + (i - 1)*64):(i*64), 92] <- rep(i, 64)
  df_policy_8[(1 + (i - 1)*64):(i*64), 92] <- rep(i, 64)
  df_policy_9[(1 + (i - 1)*64):(i*64), 92] <- rep(i, 64)
  df_policy_10[(1 + (i - 1)*64):(i*64), 92] <- rep(i, 64)
  df_policy_11[(1 + (i - 1)*64):(i*64), 92] <- rep(i, 64)
  df_policy_12[(1 + (i - 1)*64):(i*64), 92] <- rep(i, 64)
  df_policy_13[(1 + (i - 1)*64):(i*64), 92] <- rep(i, 64)
  df_policy_14[(1 + (i - 1)*64):(i*64), 92] <- rep(i, 64)
  df_policy_15[(1 + (i - 1)*64):(i*64), 92] <- rep(i, 64)
  df_policy_16[(1 + (i - 1)*64):(i*64), 92] <- rep(i, 64)
  df_policy_17[(1 + (i - 1)*64):(i*64), 92] <- rep(i, 64)
  df_policy_18[(1 + (i - 1)*64):(i*64), 92] <- rep(i, 64)
  df_policy_19[(1 + (i - 1)*64):(i*64), 92] <- rep(i, 64)
  df_policy_20[(1 + (i - 1)*64):(i*64), 92] <- rep(i, 64)

  print(i)
}

## Save results
fwrite(df_base, file = paste0("results/PSA/df_base.csv"))
fwrite(df_policy_1, file = paste0("results/PSA/df_policy_1.csv"))
fwrite(df_policy_2, file = paste0("results/PSA/df_policy_2.csv"))
fwrite(df_policy_3, file = paste0("results/PSA/df_policy_3.csv"))
fwrite(df_policy_4, file = paste0("results/PSA/df_policy_4.csv"))
fwrite(df_policy_5, file = paste0("results/PSA/df_policy_5.csv"))
fwrite(df_policy_6, file = paste0("results/PSA/df_policy_6.csv"))
fwrite(df_policy_7, file = paste0("results/PSA/df_policy_7.csv"))
fwrite(df_policy_8, file = paste0("results/PSA/df_policy_8.csv"))
fwrite(df_policy_9, file = paste0("results/PSA/df_policy_9.csv"))
fwrite(df_policy_10, file = paste0("results/PSA/df_policy_10.csv"))
fwrite(df_policy_11, file = paste0("results/PSA/df_policy_11.csv"))
fwrite(df_policy_12, file = paste0("results/PSA/df_policy_12.csv"))
fwrite(df_policy_13, file = paste0("results/PSA/df_policy_13.csv"))
fwrite(df_policy_14, file = paste0("results/PSA/df_policy_14.csv"))
fwrite(df_policy_15, file = paste0("results/PSA/df_policy_15.csv"))
fwrite(df_policy_16, file = paste0("results/PSA/df_policy_16.csv"))
fwrite(df_policy_17, file = paste0("results/PSA/df_policy_17.csv"))
fwrite(df_policy_18, file = paste0("results/PSA/df_policy_18.csv"))
fwrite(df_policy_19, file = paste0("results/PSA/df_policy_19.csv"))
fwrite(df_policy_20, file = paste0("results/PSA/df_policy_20.csv"))

### Start Here ###
df_base      <- fread(file = paste0("results/PSA/df_base.csv"))
df_policy_1  <- fread(file = paste0("results/PSA/df_policy_1.csv"))
df_policy_2  <- fread(file = paste0("results/PSA/df_policy_2.csv"))
df_policy_3  <- fread(file = paste0("results/PSA/df_policy_3.csv"))
df_policy_4  <- fread(file = paste0("results/PSA/df_policy_4.csv"))
df_policy_5  <- fread(file = paste0("results/PSA/df_policy_5.csv"))
df_policy_6  <- fread(file = paste0("results/PSA/df_policy_6.csv"))
df_policy_7  <- fread(file = paste0("results/PSA/df_policy_7.csv"))
df_policy_8  <- fread(file = paste0("results/PSA/df_policy_8.csv"))
df_policy_9  <- fread(file = paste0("results/PSA/df_policy_9.csv"))
df_policy_10 <- fread(file = paste0("results/PSA/df_policy_10.csv"))
df_policy_11 <- fread(file = paste0("results/PSA/df_policy_11.csv"))
df_policy_12 <- fread(file = paste0("results/PSA/df_policy_12.csv"))
df_policy_13 <- fread(file = paste0("results/PSA/df_policy_13.csv"))
df_policy_14 <- fread(file = paste0("results/PSA/df_policy_14.csv"))
df_policy_15 <- fread(file = paste0("results/PSA/df_policy_15.csv"))
df_policy_16 <- fread(file = paste0("results/PSA/df_policy_16.csv"))
df_policy_17 <- fread(file = paste0("results/PSA/df_policy_17.csv"))
df_policy_18 <- fread(file = paste0("results/PSA/df_policy_18.csv"))
df_policy_19 <- fread(file = paste0("results/PSA/df_policy_19.csv"))
df_policy_20 <- fread(file = paste0("results/PSA/df_policy_20.csv"))

l_policy_df <- list(df_policy_1, df_policy_2, df_policy_3,
                    df_policy_4, df_policy_5, df_policy_6,
                    df_policy_7, df_policy_8, df_policy_9,
                    df_policy_10, df_policy_11,
                    df_policy_12, df_policy_13,
                    df_policy_14, df_policy_15,
                    df_policy_16, df_policy_17,
                    df_policy_18, df_policy_19,
                    df_policy_20)
names(l_policy_df) <-  c("df_policy_1",  "df_policy_2",
                         "df_policy_3",  "df_policy_4",
                         "df_policy_5",  "df_policy_6",
                         "df_policy_7",  "df_policy_8",
                         "df_policy_9",  "df_policy_10",
                         "df_policy_11", "df_policy_12",
                         "df_policy_13", "df_policy_14",
                         "df_policy_15", "df_policy_16",
                         "df_policy_17", "df_policy_18",
                         "df_policy_19", "df_policy_20")

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

scale_all   <- function(x) (x / (n_sample * 10))
scale_diab0 <- function(x) (x / (n_diab_0 * 10))
scale_diab1 <- function(x) (x / (sum(df_sample$diab_stat) * 10))
scale_white <- function(x) (x / (n_white * 10))
scale_black <- function(x) (x / (sum(df_sample$black) * 10))
scale_hisp  <- function(x) (x / (sum(df_sample$hispanic) * 10))
scale_other <- function(x) (x / (sum(df_sample$other) * 10))
scale_65    <- function(x) (x / (sum(df_sample$can_age_at_listing_c < 5) * 10))
scale_70    <- function(x) {
  x / (sum(df_sample$can_age_at_listing_c >= 5 &
             df_sample$can_age_at_listing_c < 10) * 10)
}
scale_75    <- function(x) {
  x / (sum(df_sample$can_age_at_listing_c >= 10) * 10)
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
    as.matrix()

  return(df_pp)
}

v_colnames <- c("id", v_colnames, "run")
base_df_pp <- get_cea_results(df = df_base, v_names = v_colnames,
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
v_scenarios <- c("Base Case", "5% Increase", "10% Increase", "15% Increase",
                 "20% Increase", "25% Increase")

# df_inc_shift_05 <- data.frame(rbind(rep(0, 50), l_inc_means[[1]],
#                                     l_inc_means[[5]], l_inc_means[[9]],
#                                     l_inc_means[[13]], l_inc_means[[17]]))
# df_inc_shift_10 <- data.frame(rbind(rep(0, 50), l_inc_means[[2]],
#                                     l_inc_means[[6]], l_inc_means[[10]],
#                                     l_inc_means[[14]], l_inc_means[[18]]))
# df_inc_shift_15 <- data.frame(rbind(rep(0, 50), l_inc_means[[3]],
#                                     l_inc_means[[7]], l_inc_means[[11]],
#                                     l_inc_means[[15]], l_inc_means[[19]]))
# df_inc_shift_20 <- data.frame(rbind(rep(0, 50), l_inc_means[[4]],
#                                     l_inc_means[[8]], l_inc_means[[12]],
#                                     l_inc_means[[16]], l_inc_means[[20]]))
#
# df_inc_shift_05$v_scenarios <- v_scenarios
# df_inc_shift_10$v_scenarios <- v_scenarios
# df_inc_shift_15$v_scenarios <- v_scenarios
# df_inc_shift_20$v_scenarios <- v_scenarios

## Create dataframe for PSA analysis
base_df_pp <- as.data.frame(base_df_pp)
l_policy <- lapply(l_policy_df_clean, function(x) as.data.frame(x))

df_psa_cost_health <- data.frame(cbind(base_df_pp$tot_cost_health,
                                l_policy[[1]]$tot_cost_health,
                                l_policy[[5]]$tot_cost_health,
                                l_policy[[9]]$tot_cost_health,
                                l_policy[[13]]$tot_cost_health,
                                l_policy[[17]]$tot_cost_health))
colnames(df_psa_cost_health) <- c("cost_health_base", "cost_health_5",
                                  "cost_health_10", "cost_health_15",
                                  "cost_health_20", "cost_health_25")

df_psa_cost_soc <- data.frame(cbind(base_df_pp$tot_cost_societal,
                                l_policy[[1]]$tot_cost_societal,
                                l_policy[[5]]$tot_cost_societal,
                                l_policy[[9]]$tot_cost_societal,
                                l_policy[[13]]$tot_cost_societal,
                                l_policy[[17]]$tot_cost_societal))
colnames(df_psa_cost_soc) <- c("cost_soc_base", "cost_soc_5",
                               "cost_soc_10", "cost_soc_15",
                               "cost_soc_20", "cost_soc_25")

df_psa_eff <- data.frame(cbind(base_df_pp$tot_qalys,
                                    l_policy[[1]]$tot_qalys,
                                    l_policy[[5]]$tot_qalys,
                                    l_policy[[9]]$tot_qalys,
                                    l_policy[[13]]$tot_qalys,
                                    l_policy[[17]]$tot_qalys))
colnames(df_psa_eff) <- c("eff_base", "eff_5", "eff_10",
                          "eff_15", "eff_20", "eff_25")

l_psa <- list(df_psa_cost_health = df_psa_cost_health,
              df_psa_cost_soc = df_psa_cost_soc,
              df_psa_eff = df_psa_eff)

v_scenarios <- c("Base_Case", "Policy_5_Increase", "Policy_10_Increase",
                 "Policy_15_Increase", "Policy_20_Increase",
                 "Policy_25_Increase")

psa_obj_health <- make_psa_obj(cost = l_psa$df_psa_cost_health,
                               effectiveness = l_psa$df_psa_eff,
                               strategies = v_scenarios)
psa_obj_soc <- make_psa_obj(cost = l_psa$df_psa_cost_soc,
                            effectiveness = l_psa$df_psa_eff,
                            strategies = v_scenarios)

plot(psa_obj_health)
plot(psa_obj_soc)

### CEAC/CEAF
ceac_obj_health <- ceac(wtp = seq(from = 0, to = 100000, by = 2000),
                        psa = psa_obj_health)
ceac_obj_soc <- ceac(wtp = seq(from = 0, to = 100000, by = 2000),
                     psa = psa_obj_soc)

## Using dampack
# ceaf_health <- plot(ceac_obj_health, frontier = TRUE, points = TRUE) +
#   scale_colour_discrete(labels = c("Status Quo", "5% Increase",
#                                    "10% Increase", "15% Increase",
#                                    "20% Increase", "25% Increase"))
# ggsave(filename = "results/PSA/ceaf_health.pdf",
#        plot = ceaf_health)

# ceaf_soc <- plot(ceac_obj_soc, frontier = TRUE, points = TRUE) +
#   scale_color_nejm() +
#   scale_colour_discrete(labels = c("Status Quo", "5% Increase",
#                                    "10% Increase", "15% Increase",
#                                    "20% Increase", "25% Increase"))
# ggsave(filename = "results/PSA/ceaf_soc.pdf",
#        plot = ceaf_soc)

## Custom
ceaf_health <- plot.ceac(ceac_obj_health)
ceaf_soc    <- plot.ceac(ceac_obj_soc)

# Edit and Combine plots
ceaf_health <- ceaf_health +
  labs(title = "Healthcare Sector Perspective", tag = "A") +
  theme(title = element_text(size = 12))

ceaf_soc <- ceaf_soc +
  labs(title = "Modified Healthcare Sector Perspective", tag = "B") +
  theme(title = element_text(size = 12))

ceaf_combined <- ggarrange(ceaf_health, ceaf_soc, nrow = 1,
                           common.legend = TRUE, legend = "bottom")
ggsave(filename = "results/PSA/ceaf_combined.png",
       plot = ceaf_combined, width = 10, height = 6)
ggsave(filename = "results/PSA/ceaf_combined.pdf",
       plot = ceaf_combined, width = 10, height = 6)
### Expected Loss Curve
el_health <- calc_exp_loss(wtp = seq(from = 0, to = 100000, by = 2000),
                           psa = psa_obj_health)
el_soc <- calc_exp_loss(wtp = seq(from = 0, to = 100000, by = 2000),
                        psa = psa_obj_soc)

## Dampack
# plot_el_health <- plot(el_health, n_x_ticks = 10, n_y_ticks = 2) +
#   scale_colour_discrete(labels = c("Status Quo", "5% Increase",
#                                    "10% Increase", "15% Increase",
#                                    "20% Increase", "25% Increase")) +
#   scale_y_continuous(breaks=c(2000, 20000, 60000), trans = "log") +
#   labs(title = "Healthcare Sector Perspective", tag = "A")
#
#
# plot_el_soc <- plot(el_soc, n_x_ticks = 10, n_y_ticks = 2) +
#   scale_colour_discrete(labels = c("Status Quo", "5% Increase",
#                                    "10% Increase", "15% Increase",
#                                    "20% Increase", "25% Increase")) +
#   scale_y_continuous(breaks=c(5000, 20000, 60000), trans = "log") +
#   labs(title = "Modified Healthcare Sector Perspective", tag = "B")

## Custom
plot_el_health <- plot.exp_loss(el_health, n_x_ticks = 10, n_y_ticks = 3) +
  scale_y_continuous(breaks=c(2000, 10000, 60000),
                     labels = labfun, trans = "log") +
  labs(title = "Healthcare Sector Perspective", tag = "A") +
  theme(title = element_text(size = 12))

plot_el_soc <- plot.exp_loss(el_soc, n_x_ticks = 10, n_y_ticks = 3) +
  scale_y_continuous(breaks=c(2000, 10000, 60000),
                     labels = labfun, trans = "log") +
  labs(title = "Modified Healthcare Sector Perspective", tag = "B") +
  theme(title = element_text(size = 12))

# Edit and combine plots
el_combined <- ggarrange(plot_el_health, plot_el_soc, nrow = 1,
                           common.legend = TRUE, legend = "bottom")
ggsave(filename = "results/PSA/el_combined.pdf",
       plot = el_combined, width = 10, height = 6)
ggsave(filename = "results/PSA/el_combined.png",
       plot = el_combined, width = 10, height = 6)
## Incremental Net Monetary Benefit (INMB)
wtp <- 100000
df_inmb_health_psa <- cbind(psa_obj_health$cost, psa_obj_health$effectiveness) %>%
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
  mutate(nmb_societal_all = (tot_qalys * wtp) - tot_cost_societal)
