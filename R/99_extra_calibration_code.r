n.top <- n
all_df_post <- lazy_dt(all_df) %>%
  filter(ddtx==1) %>%
  filter(id<=n.top | id==n + 1) %>%
  as_tibble()

# fit_gl_death <- survfit(Surv(time = months_to_gl_mort, gl_death) ~ id,
#                         data = all_df_post)  # Death after graft loss
df_post_id <- split(all_df_post, all_df_post$id)
l_fits <- list()
for (i in 1:length(df_post_id)) {
  l_fits[[i]] <- survfit(Surv(time = months_to_gl_mort, gl_death) ~ 1,
                          data = df_post_id[[i]])
  print(i)
}

split_all <- list()
for (i in 1:length(l_fits)) {
  split_all[[i]] <- cbind(l_fits[[i]][["time"]], l_fits[[i]][["surv"]])
}

v_trim <- rep(F, length(split_all))
for (i in 1:length(split_all)) {
  tmp_df <- as.data.frame(split_all[[i]])
  tmp_df <- mutate(tmp_df, 
                   trim = ifelse(V1 <= 5 & V2 < 0.25, T, 
                                 ifelse(V1 <= 10 & V2 < 0.18, T,
                                        ifelse(V1 >= 30 & V2 > 0.45, T,
                                               ifelse(V1 >= 45 & V2 > 0.3, T, F))))) 
  v_trim[i] <- any(tmp_df$trim)
}

# Test 
# split_points <- rep(F, length(fit_gl_death[["time"]]))
# for (i in 1:length(fit_gl_death[["time"]])) {
#   if (i != 1) {
#     if (fit_gl_death[["time"]][i] <= fit_gl_death[["time"]][i-1] | 
#         fit_gl_death[["time"]][i] == 0) {
#       split_points[i] <- T
#     }
#   }
# }
# split_points <- c(1, which(split_points))
# 
# # create split vector:
# split_code <- rep(1, length(fit_gl_death[["time"]]))
# for ( j in 1:length(split_points) ) {
#   
#   if (j!=length(split_points)) {
#     split_code[
#       split_points[j]:(split_points[j+1]-1)
#     ] <- j
#   } else {
#     split_code[
#       split_points[j]:length(fit_gl_death[["time"]])
#     ] <- j
#   }
#   
# }
# 
# split_times <- split(fit_gl_death[["time"]], split_code)
# split_surv  <- split(fit_gl_death[["surv"]], split_code)
# 
# split_all <- list()
# for (i in 1:length(split_surv)) {
#   split_all[[i]] <- rbind(split_times[[i]], split_surv[[i]])
# }
# 
# v_trim <- rep(F, length(split_all))
# for (i in 1:length(split_all)) {
#   if (ncol(split_all[[i]]) < 10) {
#     v_trim[i] <- T
#   } 
#   if (ncol(split_all[[i]]) >= 10) {
#     if (split_all[[i]][2,10] < 0.2) {
#       v_trim[i] <- T
#     }
#   } 
#   if (ncol(split_all[[i]]) >= 30){
#     if (split_all[[i]][2,30] > 0.5) {
#       v_trim[i] <- T
#     }
#   }
# }

v_trim_id <- which(v_trim %in% T)

all_df_post <- lazy_dt(all_df) %>%
  filter(ddtx==1) %>%
  filter(!(id %in% v_trim_id) | id==n + 1) %>%
  as_tibble()

v_pos_id <- unique(all_df_post$id)[-1]
v_good_id <- v_id[v_pos_id]

m_gof_trim <- lazy_dt(df_gof) %>%
  filter(id %in% v_good_id) %>%
  as_tibble()

write.csv(m_gof_trim, file = paste0(serv_path, "gof_recalibration_trim.csv"))

n.top <- length(unique(all_df_post$id)) - 1
top_group <- c(n.top, n + 1)
col_vector <- c(rep("gray", n.top), "red")

#KM curves
fit_gl <- survfit(Surv(time = months_to_gl, graft_loss) ~ id,
                  conf.int=.99,
                  data = all_df_post)  # Graft loss

fit_gs_death <- survfit(Surv(time = months_to_gs_mort, gs_death) ~ id,
                        conf.int=.99,
                        data = all_df_post)  # Death w/ function

fit_gl_death <- survfit(Surv(time = months_to_gl_mort, gl_death) ~ id,
                        conf.int=.99,
                        data = all_df_post)  # Death after graft loss

# Generate plots
gl_all <- surv_plot(fit = fit_gl, title_txt = "Graft Loss", n_top = n.top)
gs_death_all <- surv_plot(fit = fit_gs_death, title_txt = "Death w/ Function", n_top = n.top)
gl_death_all <- surv_plot(fit = fit_gl_death, title_txt = "Death after Graft Loss", n_top = n.top)

plots_post_trim <- grid.arrange(arrangeGrob(gl_all$plot, gs_death_all$plot,
                                            gl_death_all$plot +
                                              theme(legend.position="none"),
                                            nrow = 2),
                                nrow=2,heights=c(10, 1))
ggsave(filename = paste0(serv_path, "recalibration_post_trim.pdf"),
       plot = plots_post_trim, width = 7, height = 4)

# Proportions
df_mlogit_all_long <- summary_df_mlogit_long(data = all_df_post, n = n.top)
df_glmort_all_long <- summary_df_glmort_long(data = all_df_post, n = n.top)

mlogit_trim_plot <- mlogit_plot(data = df_mlogit_all_long, n = n.top)
glmort_trim_plot <- glmort_plot(data = df_glmort_all_long, n = n.top)

gl_trim <- list(mlogit_trim_plot, glmort_trim_plot)
plot_prop_trim <- gridExtra::grid.arrange(grobs=gl_trim,
                                    layout_matrix = lay)
ggsave(filename = paste0(serv_path, "recalibration_prop_trim_v2.pdf"),
       plot = plot_prop_trim, width = 7, height = 4)


# Test surv_plot mods
surv_plot <- function(fit, title_txt, n_top) {
  tmp_plot <- ggsurvplot(fit,
                         censor = FALSE,
                         conf.int = TRUE,
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
    scale_fill_manual(values = c(rep("gray", n_top), "red")) +
    guides(color = guide_none())
  
  
  return(tmp_plot)
}

surv_plot_gl_death <- function(fit, title_txt, n_top) {
  tmp_plot <- ggsurvplot(fit,
                         censor = FALSE,
                         conf.int = TRUE,
                         title = title_txt,
                         short.panel.labs = TRUE,
                         legend.title = "",
                         xlim = c(0, 42),
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
    scale_fill_manual(values = c(rep("gray", n_top), "red")) +
    guides(color = guide_none()) + 
    scale_x_continuous(labels = c(0,10,20,30,40))
  
  
  return(tmp_plot)
}

gl_all <- surv_plot(fit = fit_gl, title_txt = "Graft Loss", n_top = n.top)
p_build <- ggplot_build(gl_all$plot)
p_build$data[[1]] <- p_build$data[[1]] %>% 
  mutate(fill = case_when(colour == "gray" ~ "gray",
                          colour == "red" ~ "red"))
p_build$data[[3]] <- p_build$data[[3]] %>% 
  mutate(alpha = case_when(fill == "gray" ~ 0,
                           fill == "red" ~ 0.3))
gl_all <- ggplot_gtable(p_build)

gs_death_all <- surv_plot(fit = fit_gs_death, title_txt = "Death w/ Function", n_top = n.top)
p_build <- ggplot_build(gs_death_all$plot)
p_build$data[[1]] <- p_build$data[[1]] %>% 
  mutate(fill = case_when(colour == "gray" ~ "gray",
                          colour == "red" ~ "red"))
p_build$data[[3]] <- p_build$data[[3]] %>% 
  mutate(alpha = case_when(fill == "gray" ~ 0,
                           fill == "red" ~ 0.3))
gs_death_all <- ggplot_gtable(p_build)

gl_death_all <- surv_plot_gl_death(fit = fit_gl_death, title_txt = "Death after Graft Loss", n_top = n.top)
p_build <- ggplot_build(gl_death_all$plot)
p_build$data[[1]] <- p_build$data[[1]] %>% 
  mutate(fill = case_when(colour == "gray" ~ "gray",
                          colour == "red" ~ "red"))
p_build$data[[3]] <- p_build$data[[3]] %>% 
  mutate(alpha = case_when(fill == "gray" ~ 0,
                           fill == "red" ~ 0.3))
gl_death_all <- ggplot_gtable(p_build)

grid.arrange(arrangeGrob(gl_all, 
                         gs_death_all,
                         gl_death_all, nrow = 2))





# Export coefficients
m_ddtx_trim      <- m_ddtx_coef %>% filter(id %in% v_good_id) %>% arrange(factor(id, levels = v_good_id))
m_ldtx_trim      <- m_ldtx_coef %>% filter(id %in% v_good_id) %>% arrange(factor(id, levels = v_good_id))
m_mort_trim      <- m_mort_coef %>% filter(id %in% v_good_id) %>% arrange(factor(id, levels = v_good_id))
m_remove_trim    <- m_remove_coef %>% filter(id %in% v_good_id) %>% arrange(factor(id, levels = v_good_id))
m_mlogit_trim    <- m_mlogit_coef %>% filter(id %in% v_good_id) %>% arrange(factor(id, levels = v_good_id))
m_gs_mort_trim   <- m_gs_mort_coef %>% filter(id %in% v_good_id) %>% arrange(factor(id, levels = v_good_id))
m_gl_trim        <- m_gl_coef %>% filter(id %in% v_good_id) %>% arrange(factor(id, levels = v_good_id))
m_dial_mort_trim <- m_dial_mort_coef %>% filter(id %in% v_good_id) %>% arrange(factor(id, levels = v_good_id))
m_gl_mort_trim   <- m_gl_mort_coef %>% filter(id %in% v_good_id) %>% arrange(factor(id, levels = v_good_id))

m_ddtx_coef_trim        <- m_ddtx_trim   %>% dplyr::select(id:V27) %>% as.matrix()
m_ldtx_coef_trim        <- m_ldtx_trim   %>% dplyr::select(id:V27) %>% as.matrix()
m_mort_coef_trim        <- m_mort_trim   %>% dplyr::select(id:V27) %>% as.matrix()
m_remove_coef_trim      <- m_remove_trim %>% dplyr::select(id:V27) %>% as.matrix()
m_mlogit_coef_trim      <- m_mlogit_trim %>% as.matrix()
m_gs_mort_coef_trim     <- m_gs_mort_trim %>% dplyr::select(V1:V18, id) %>% as.matrix()
m_gl_coef_trim          <- m_gl_trim %>% dplyr::select(V1:V19, id) %>% as.matrix()
m_dial_mort_coef_trim   <- m_dial_mort_trim %>% dplyr::select(V1:V20, id) %>% as.matrix()
m_gl_mort_coef_trim     <- m_gl_mort_trim %>% dplyr::select(V1:V20, id) %>% as.matrix()

m_ddtx_shape_trim       <- m_ddtx_trim      %>% dplyr::select(id, V28) %>% as.matrix()
m_ldtx_shape_trim       <- m_ldtx_trim      %>% dplyr::select(id, V28) %>% as.matrix()
m_mort_shape_trim       <- m_mort_trim      %>% dplyr::select(id, V28) %>% as.matrix()
m_remove_shape_trim     <- m_remove_trim    %>% dplyr::select(id, V28) %>% as.matrix()
m_gs_mort_shape_trim    <- m_gs_mort_trim   %>% dplyr::select(id, V19) %>% as.matrix()
m_gl_shape_trim         <- m_gl_trim        %>% dplyr::select(id, V20) %>% as.matrix()
m_dial_mort_shape_trim  <- m_dial_mort_trim %>% dplyr::select(id, V21) %>% as.matrix()

# Export code
fwrite(m_ddtx_coef_trim, file = paste0(serv_path, "ddtx_coef_trim.csv"))
fwrite(m_ldtx_coef_trim, file = paste0(serv_path, "ldtx_coef_trim.csv"))
fwrite(m_mort_coef_trim, file = paste0(serv_path, "mort_coef_trim.csv"))
fwrite(m_remove_coef_trim, file = paste0(serv_path, "remove_coef_trim.csv"))
fwrite(m_mlogit_coef_trim, file = paste0(serv_path, "mlogit_coef_trim.csv"))
fwrite(m_gs_mort_coef_trim, file = paste0(serv_path, "gs_mort_coef_trim.csv"))
fwrite(m_gl_coef_trim, file = paste0(serv_path, "gl_coef_trim.csv"))
fwrite(m_dial_mort_coef_trim, file = paste0(serv_path, "dial_mort_coef_trim.csv"))
fwrite(m_gl_mort_coef_trim, file = paste0(serv_path, "gl_mort_coef_trim.csv"))
fwrite(m_ddtx_shape_trim, file = paste0(serv_path, "ddtx_shape_trim.csv"))
fwrite(m_ldtx_shape_trim, file = paste0(serv_path, "ldtx_shape_trim.csv"))
fwrite(m_mort_shape_trim, file = paste0(serv_path, "mort_shape_trim.csv"))
fwrite(m_remove_shape_trim, file = paste0(serv_path, "remove_shape_trim.csv"))
fwrite(m_gs_mort_shape_trim, file = paste0(serv_path, "gs_mort_shape_trim.csv"))
fwrite(m_gl_shape_trim, file = paste0(serv_path, "gl_shape_trim.csv"))
fwrite(m_dial_mort_shape_trim, file = paste0(serv_path, "dial_mort_shape_trim.csv"))

