# The following is all code to produce a figure about the distribution of donor
# kidney quality for SMDM 2023 poster and presentation

samp_sq <- as.data.frame(rand_kdpi_cold(num = 100000))
samp_p  <- as.data.frame(rand_kdpi_cold(num = 100000, france = TRUE, pct_shift = 0.05))

library(ggplot2)
library(ggsci)
library(gridExtra)

samp_sq$id <- "Status Quo"
samp_p$id <- "Policy"
df_samp <- rbind(samp_sq, samp_p)
df_samp$id <- factor(df_samp$id, levels = c("Status Quo", "Policy"))

plot <- ggplot(data = df_samp, mapping = aes(x = kdpi)) +
  geom_histogram(binwidth = 1, aes(y = ..density..), fill = "#0072B5FF") +
  facet_grid(id~.) +
  theme_bw() +
  scale_fill_nejm() +
  xlab("Kidney Quality: KDPI") +
  ylab("Density") +
  theme(
    legend.position = "none",
    strip.text = element_text(size = 13),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    axis.title.x = element_text(size = 16, face = "bold"),
    axis.title.y = element_text(size = 16, face = "bold")
  )
plot
ggsave("results/dist.pdf", plot, width = 8, height = 7)



## Compare shift distributions
samp_05 <- as.data.frame(rand_kdpi_cold(num = 100000, france = TRUE, pct_shift = 0.05))
samp_20 <- as.data.frame(rand_kdpi_cold(num = 100000, france = TRUE, pct_shift = 0.20))

samp_05$id <- "5%"
samp_20$id <- "20%"
df_samp_shift <- rbind(samp_05, samp_20)
df_samp_shift$id <- factor(df_samp_shift$id, levels = c("5%", "20%"))
# fill = "#0072B5FF"
plot <- ggplot(data = df_samp_shift, mapping = aes(x = kdpi, fill = id)) +
  geom_density(alpha = 0.5) +
  # facet_grid(id~.) +
  theme_bw() +
  scale_fill_nejm() +
  xlab("Kidney Quality: KDPI") +
  ylab("Density") +
  theme(
    legend.position = "right",
    strip.text = element_text(size = 13),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    axis.title.x = element_text(size = 16, face = "bold"),
    axis.title.y = element_text(size = 16, face = "bold")
  )
plot

library(dplyr)
df_samp_shift %>% group_by(id) %>%
  # calculate densities for each group over same range; store in list column
  summarise(d = list(density(kdpi, from = min(.$kdpi), to = max(.$kdpi)))) %>%
  # make a new data.frame from two density objects
  do(data.frame(x = .$d[[1]]$x,    # grab one set of x values (which are the same)
                y = .$d[[1]]$y - .$d[[2]]$y)) %>%    # and subtract the y values
  ggplot(aes(x, y)) +    # now plot
  geom_line()
ggsave("results/dist.pdf", plot, width = 8, height = 7)
