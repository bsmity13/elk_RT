############################################X
#--------------Elk Risk Taking--------------X
#---------------Brian J. Smith--------------X
#---------------Visualize iSSA--------------X
############################################X

# Set options ----
# For full printing of tibble columns
options(pillar.width = 1000)

# Load packages ----
library(tidyverse)
library(lubridate)
library(amt)
library(ragg)
library(RColorBrewer)
library(lubridate)
library(patchwork)

# Custom functions
source("99_fun.R")

# Create directory ----
# Repeating this across scripts in case they are not all run
dir.create("fig", showWarnings = FALSE)

# Options ----
theme_set(theme_bw())

# Load data ----
# iSSA data
issa_dat <- readRDS("../out/final_nested_data.rds")
# Fitted iSSA
dat <- readRDS("../out/fitted_issa.rds")
# Sightability-corrected Population density
# from Metz et al. monograph (appendix B)
pop <- read.csv("../elk_RT_data/Metz_etal_ssm_estimates.csv") %>% 
  # Keep just the population size inside of YNP
  filter(where == "inside") %>% 
  select(winter, elk_lwr = lwr, elk_mean = mean, elk_upr = upr)
# Predator density
wolf <- read.csv("../elk_RT_data/wolf_density.csv") %>% 
  select(winter, wolf = Density)
cougar <- read.csv("../elk_RT_data/cougar_density.csv") %>% 
  select(winter, cougar = Density)
# Scaling data.frame
scale_df <- read.csv("../out/scaling_data.csv")
row.names(scale_df) <- scale_df$term

# Correlation between food and risk (available) ----
dat2 <- issa_dat %>% 
  select(covs) %>% 
  unnest(cols = covs) %>% 
  filter(!case_)
cor(dat2$bio_end, dat2$wolf_risk_end, method = "spearman")
cor(dat2$bio_end, dat2$cougar_risk_end, method = "spearman")

# Correlation between forage and risk ----
# Hex-bin plot

hb_wolf <- ggplot(dat2, aes(x = bio_end, y = wolf_risk_end)) +
  geom_hex(binwidth = c(1, 0.35)) +
  scale_fill_viridis_c(name = "# Steps") +
  xlab("Forage Biomass (SD)") +
  ylab("Wolf Risk (SD)") +
  theme(legend.position = "none")

hb_cougar <- ggplot(dat2, aes(x = bio_end, y = cougar_risk_end)) +
  geom_hex(binwidth = c(1, 0.8)) +
  scale_fill_viridis_c(name = "# Steps") +
  xlab("Forage Biomass (SD)") +
  ylab("Cougar Risk (SD)") +
  theme(legend.position = "none")

hb_plot <- hb_wolf + hb_cougar +
  plot_layout(ncol = 2, nrow = 1) +
  plot_annotation(tag_levels = "A", tag_prefix = "(", tag_suffix = ")")

ggsave("fig/forage_risk_hexbin.tif", plot = hb_plot, device = agg_tiff, 
       width = 6.5, height = 6.5, units = "in", dpi = 300, compression = "lzw")


# Bootstrap ----
# Resample beta coefficients to account for uncertainty
set.seed(20221214)
boot_beta <- bootstrap_betas(data = dat, col = "issf", iter = 2000)

# Separate ID and winter
boot_beta <- split_ID_winter_season(boot_beta)

# Save
saveRDS(boot_beta, "../out/bootstrapped_iSSA_coefs.rds")

# If you want to load
# boot_beta <- readRDS("../out/bootstrapped_iSSA_coefs.rds")


# Habitat coefficients ----
# Population mean coefficient across all winters for each bootstrap iteration
beta_means <- boot_beta %>% 
  group_by(term, iter, season) %>% 
  summarize(mean_est = mean(estimate))

(coef_plot_early <- beta_means %>% 
    filter(season == "early") %>% 
    # Habitat selection coefficients only
    filter(!str_detect(term, fixed("cos(ta_)")),
           !str_detect(term, fixed("log(sl_)"))) %>% 
    ggplot(aes(x = term, y = mean_est, group = term)) +
    geom_violin() +
    geom_hline(yintercept = 0, color = "red", linetype = "dashed") +
    xlab(NULL) +
    ylab("Estimate") +
    ggtitle("Early Winter") +
    coord_cartesian(ylim = c(-1.5, 1.5)) +
    theme(axis.text.x = element_blank()) +
    NULL)

(coef_plot_late <- beta_means %>% 
    filter(season == "late") %>% 
    # Habitat selection coefficients only
    filter(!str_detect(term, fixed("cos(ta_)")),
           !str_detect(term, fixed("log(sl_)"))) %>% 
    ggplot(aes(x = term, y = mean_est, group = term)) +
    geom_violin() +
    geom_hline(yintercept = 0, color = "red", linetype = "dashed") +
    xlab("Coefficient") +
    ylab("Estimate") +
    ggtitle("Late Winter") +
    coord_cartesian(ylim = c(-1.5, 1.5)) +
    theme(axis.text.x = element_text(angle = 30, hjust = 1)) +
    NULL)

coef_plot <- coef_plot_early + coef_plot_late +
  plot_layout(nrow = 2, ncol = 1, heights = c(0.4, 0.6)) +
  plot_annotation(tag_levels = "A", tag_prefix = "(", tag_suffix = ")")

ggsave("fig/issa_coefs.tif", plot = coef_plot, device = agg_tiff, 
       width = 7, height = 7, units = "in", dpi = 300, compression = "lzw")

# Fig. 3A ----
# biomass selection by time of day

# ... ... setup x1 and x2 ----
# Predict for early winter using December 1
t_early <- seq(ymd_hms("2021-12-01 00:00:00", tz = "US/Mountain"), 
               ymd_hms("2021-12-01 23:59:00", tz = "US/Mountain"),
               by = "15 min")

# Predict for late winter using March 15
t_late <- seq(ymd_hms("2022-03-15 00:00:00", tz = "US/Mountain"), 
              ymd_hms("2022-03-15 23:59:00", tz = "US/Mountain"),
              by = "15 min")

bio_early_x1 <- create_pred_dat(scale_df = scale_df,
                                start_time = t_early,
                                # One SD above the mean
                                bio_end = scale_df["bio_end", "mean"] +
                                  scale_df["bio_end", "sd"]) %>% 
  scale_dat(scale_df)

bio_late_x1 <- create_pred_dat(scale_df = scale_df,
                               start_time = t_late,
                               # One SD above the mean
                               bio_end = scale_df["bio_end", "mean"] +
                                 scale_df["bio_end", "sd"]) %>% 
  scale_dat(scale_df)

bio_early_x2 <- create_pred_dat(scale_df = scale_df,
                                start_time = t_early,
                                # The mean
                                bio_end = scale_df["bio_end", "mean"]) %>% 
  scale_dat(scale_df)

bio_late_x2 <- create_pred_dat(scale_df = scale_df,
                               start_time = t_late,
                               # The mean
                               bio_end = scale_df["bio_end", "mean"]) %>% 
  scale_dat(scale_df)


bio_x1 <- rbind(bio_early_x1, bio_late_x1)
bio_x2 <- rbind(bio_early_x2, bio_late_x2)

# ... ... calculate RSS for bootstrap iterations ----
system.time({ # 3090 s = 51.5 min
  set.seed(123456)
  lr_bio <- boot_rss(model_df = dat,
                     niter = 2000,
                     x1 = bio_x1,
                     x2 = bio_x2)
})

# ... ... population mean ----
bio_pop_mean <- lr_bio %>% 
  group_by(iter, season, start_time) %>% 
  summarize(log_rss = mean(log_rss, na.rm = TRUE))

# ... ... bootstrap mean ----
bio_boot_mean <- bio_pop_mean %>% 
  group_by(season, start_time) %>% 
  summarize(log_rss = mean(log_rss, na.rm = TRUE))

# ... ... plot ----
# Early winter
bio_plot_early <- bio_pop_mean %>% 
  filter(season == "early") %>% 
  ggplot(mapping = aes(x = start_time, y = log_rss)) +
  # Left-side night
  geom_rect(aes(xmin = as.POSIXct("2021-12-01 00:00:00", tz = "US/Mountain"), 
                xmax = as.POSIXct("2021-12-01 07:41:00", tz = "US/Mountain"),  
                ymin = -Inf, 
                ymax = Inf),
            fill = "gray90", alpha = 0.1, color = NA) +
  # Right-side night
  geom_rect(aes(xmin = as.POSIXct("2021-12-01 16:43:00", tz = "US/Mountain"), 
                xmax = as.POSIXct("2021-12-02 00:00:00", tz = "US/Mountain"),
                ymin = -Inf, 
                ymax = Inf),
            fill = "gray90", alpha = 0.1, color = NA) +
  geom_line(aes(group = iter), color = "black", alpha = 0.02) +
  geom_line(data = filter(bio_boot_mean, season == "early"), 
            color = "#bbffbb", linewidth = 1) +
  geom_hline(yintercept = 0, color = "red", linetype = "dashed") +
  xlab(NULL) +
  ylab("log-RSS Biomass") +
  scale_x_datetime(date_labels = "%H:%M") +
  coord_cartesian(ylim = c(-0.5, 0.6), expand = FALSE) +
  # ggtitle("Early Winter") +
  theme(axis.text.x = element_blank())

# Late winter
bio_plot_late <- bio_pop_mean %>% 
  filter(season == "late") %>% 
  ggplot(mapping = aes(x = start_time, y = log_rss)) +
  # Left-side night
  geom_rect(aes(xmin = as.POSIXct("2022-03-15 00:00:00", tz = "US/Mountain"), 
                xmax = as.POSIXct("2022-03-15 07:35:00", tz = "US/Mountain"), 
                ymin = -Inf, 
                ymax = Inf),
            fill = "gray90", alpha = 0.1, color = NA) +
  # Right-side night
  geom_rect(aes(xmin = as.POSIXct("2022-03-15 19:29:00", tz = "US/Mountain"), 
                xmax = as.POSIXct("2022-03-16 00:00:00", tz = "US/Mountain"), 
                ymin = -Inf, 
                ymax = Inf),
            fill = "gray90", alpha = 0.1, color = NA) +
  geom_line(aes(group = iter), color = "black", alpha = 0.02) +
  geom_line(data = filter(bio_boot_mean, season == "late"), 
            color = "#bbffbb", linewidth = 1) +
  geom_hline(yintercept = 0, color = "red", linetype = "dashed") +
  xlab("Step Start Time") +
  ylab("log-RSS Biomass") +
  scale_x_datetime(date_labels = "%H:%M") +
  coord_cartesian(ylim = c(-0.5, 0.6), expand = FALSE) + 
  # ggtitle("Late Winter")
  NULL

# Combine
bio_plot <- bio_plot_early + bio_plot_late +
  plot_layout(nrow = 2, ncol = 1, heights = c(0.4, 0.6)) +
  plot_annotation(tag_levels = "A", tag_prefix = "(", tag_suffix = ")")

# ggsave("fig/biomass_by_time.tif", plot = bio_plot, device = agg_tiff,
#        width = 9, height = 7, units = "in", compression = "lzw")

# Fig. 3B ----
# Risk selection by time of day
# ... ... setup x1 and x2 ----
# Note that predator activities are only estimated on the hour
# Predict for early winter using December 1
t_early2 <- seq(ymd_hms("2021-12-01 00:00:00", tz = "US/Mountain"), 
                ymd_hms("2021-12-02 00:00:00", tz = "US/Mountain"),
                by = "1 hour")

# Predict for late winter using March 15
t_late2 <- seq(ymd_hms("2022-03-15 00:00:00", tz = "US/Mountain"), 
               ymd_hms("2022-03-16 00:00:00", tz = "US/Mountain"),
               by = "1 hour")

# Wolf, early
wolf_risk_early_x1 <- create_pred_dat(scale_df = scale_df,
                                      start_time = t_early2,
                                      # One SD above the mean
                                      wolf_risk_end = scale_df["wolf_risk_end", "mean"] +
                                        scale_df["wolf_risk_end", "sd"]) %>% 
  scale_dat(scale_df)  %>% 
  mutate(predator = "wolf")

wolf_risk_early_x2 <- create_pred_dat(scale_df = scale_df,
                                      start_time = t_early2,
                                      # The mean
                                      wolf_risk_end = scale_df["wolf_risk_end", "mean"]) %>% 
  scale_dat(scale_df) %>% 
  mutate(predator = "wolf")

# Wolf, late
wolf_risk_late_x1 <- create_pred_dat(scale_df = scale_df,
                                     start_time = t_late2,
                                     # One SD above the mean
                                     wolf_risk_end = scale_df["wolf_risk_end", "mean"] +
                                       scale_df["wolf_risk_end", "sd"]) %>% 
  scale_dat(scale_df) %>% 
  mutate(predator = "wolf")

wolf_risk_late_x2 <- create_pred_dat(scale_df = scale_df,
                                     start_time = t_late2,
                                     # The mean
                                     wolf_risk_end = scale_df["wolf_risk_end", "mean"]) %>% 
  scale_dat(scale_df) %>% 
  mutate(predator = "wolf")

# Cougar, early
cougar_risk_early_x1 <- create_pred_dat(scale_df = scale_df,
                                        start_time = t_early2,
                                        # One SD above the mean
                                        cougar_risk_end = scale_df["cougar_risk_end", "mean"] +
                                          scale_df["cougar_risk_end", "sd"]) %>% 
  scale_dat(scale_df) %>% 
  mutate(predator = "cougar")

cougar_risk_early_x2 <- create_pred_dat(scale_df = scale_df,
                                        start_time = t_early2,
                                        # The mean
                                        cougar_risk_end = scale_df["cougar_risk_end", "mean"]) %>% 
  scale_dat(scale_df) %>% 
  mutate(predator = "cougar")

# Cougar, late
cougar_risk_late_x1 <- create_pred_dat(scale_df = scale_df,
                                       start_time = t_late2,
                                       # One SD above the mean
                                       cougar_risk_end = scale_df["cougar_risk_end", "mean"] +
                                         scale_df["cougar_risk_end", "sd"]) %>% 
  scale_dat(scale_df) %>% 
  mutate(predator = "cougar")

cougar_risk_late_x2 <- create_pred_dat(scale_df = scale_df,
                                       start_time = t_late2,
                                       # The mean
                                       cougar_risk_end = scale_df["cougar_risk_end", "mean"]) %>% 
  scale_dat(scale_df) %>% 
  mutate(predator = "cougar")


# wolf_risk_x1 <- rbind(wolf_risk_early_x1,
#                       wolf_risk_late_x1)
# 
# wolf_risk_x2 <- rbind(wolf_risk_early_x2,
#                       wolf_risk_late_x2)
# 
# cougar_risk_x1 <- rbind(cougar_risk_early_x1,
#                       cougar_risk_late_x1)
# 
# cougar_risk_x2 <- rbind(cougar_risk_early_x2,
#                       cougar_risk_late_x2)

risk_x1 <- rbind(wolf_risk_early_x1,
                 wolf_risk_late_x1,
                 cougar_risk_early_x1,
                 cougar_risk_late_x1)

risk_x2 <- rbind(wolf_risk_early_x2,
                 wolf_risk_late_x2,
                 cougar_risk_early_x2,
                 cougar_risk_late_x2)

# ... ... calculate RSS for bootstrap iterations ----
system.time({
  set.seed(123456)
  lr_risk <- boot_rss(model_df = dat,
                      niter = 2000,
                      x1 = risk_x1,
                      x2 = risk_x2)
})
# 
# system.time({
#   set.seed(123456)
#   lr_risk_wolf <- boot_rss(model_df = dat,
#                            niter = 30,
#                            x1 = wolf_risk_x1,
#                            x2 = wolf_risk_x2)
# 
#   lr_risk_cougar <- boot_rss(model_df = dat,
#                              niter = 30,
#                              x1 = cougar_risk_x1,
#                              x2 = cougar_risk_x2)
# })
# 
# lr_risk <- rbind(lr_risk_wolf, lr_risk_cougar)

# ... ... population mean ----
risk_pop_mean <- lr_risk %>% 
  group_by(iter, season, predator, start_time) %>% 
  summarize(log_rss = mean(log_rss, na.rm = TRUE)) %>% 
  mutate(pred_iter = paste0(predator, iter))

# ... ... bootstrap mean ----
risk_boot_mean <- risk_pop_mean %>% 
  group_by(season, predator, start_time) %>% 
  summarize(log_rss = mean(log_rss, na.rm = TRUE))

# ... ... plot ----
# Early winter
risk_plot_early_wolf <- risk_pop_mean %>% 
  filter(season == "early", predator == "wolf") %>% 
  ggplot(mapping = aes(x = start_time, y = log_rss)) +
  # Left-side night
  geom_rect(aes(xmin = as.POSIXct("2021-12-01 00:00:00", tz = "US/Mountain"), 
                xmax = as.POSIXct("2021-12-01 07:41:00", tz = "US/Mountain"),  
                ymin = -Inf, 
                ymax = Inf),
            fill = "gray90", alpha = 0.1, color = NA) +
  # Right-side night
  geom_rect(aes(xmin = as.POSIXct("2021-12-01 16:43:00", tz = "US/Mountain"), 
                xmax = as.POSIXct("2021-12-01 23:59:59", tz = "US/Mountain"),
                ymin = -Inf, 
                ymax = Inf),
            fill = "gray90", alpha = 0.1, color = NA) +
  geom_line(aes(group = iter), alpha = 0.02) +
  geom_line(data = filter(risk_boot_mean, 
                          season == "early",
                          predator == "wolf"),
            color = "#993333",
            linewidth = 1) +
  geom_hline(yintercept = 0, color = "red", linetype = "dashed") +
  xlab(NULL) +
  ylab("log-RSS Risk") +
  scale_x_datetime(date_labels = "%H:%M") +
  coord_cartesian(ylim = c(-0.5, 0.6), expand = FALSE) +
  # ggtitle("Wolf, Early Winter") +
  theme(axis.text.x = element_blank())

# Late winter
risk_plot_late_wolf <- risk_pop_mean %>% 
  filter(season == "late", predator == "wolf") %>% 
  ggplot(mapping = aes(x = start_time, y = log_rss)) +
  # Left-side night
  geom_rect(aes(xmin = as.POSIXct("2022-03-15 00:00:00", tz = "US/Mountain"), 
                xmax = as.POSIXct("2022-03-15 07:35:00", tz = "US/Mountain"), 
                ymin = -Inf, 
                ymax = Inf),
            fill = "gray90", alpha = 0.1, color = NA) +
  # Right-side night
  geom_rect(aes(xmin = as.POSIXct("2022-03-15 19:29:00", tz = "US/Mountain"), 
                xmax = as.POSIXct("2022-03-15 23:59:59", tz = "US/Mountain"), 
                ymin = -Inf, 
                ymax = Inf),
            fill = "gray90", alpha = 0.1, color = NA) +
  geom_line(aes(group = iter), alpha = 0.02) +
  geom_line(data = filter(risk_boot_mean, 
                          season == "late",
                          predator == "wolf"),
            color = "#993333",
            linewidth = 1) +
  geom_hline(yintercept = 0, color = "red", linetype = "dashed") +
  xlab("Step Start Time") +
  ylab("log-RSS Risk") +
  scale_x_datetime(date_labels = "%H:%M") +
  coord_cartesian(ylim = c(-0.5, 0.6), expand = FALSE) +
  # ggtitle("Wolf, Late Winter") +
  NULL

# Early winter
risk_plot_early_cougar <- risk_pop_mean %>% 
  filter(season == "early", predator == "cougar") %>% 
  ggplot(mapping = aes(x = start_time, y = log_rss)) +
  # Left-side night
  geom_rect(aes(xmin = as.POSIXct("2021-12-01 00:00:00", tz = "US/Mountain"), 
                xmax = as.POSIXct("2021-12-01 07:41:00", tz = "US/Mountain"),  
                ymin = -Inf, 
                ymax = Inf),
            fill = "gray90", alpha = 0.1, color = NA) +
  # Right-side night
  geom_rect(aes(xmin = as.POSIXct("2021-12-01 16:43:00", tz = "US/Mountain"), 
                xmax = as.POSIXct("2021-12-01 23:59:59", tz = "US/Mountain"),
                ymin = -Inf, 
                ymax = Inf),
            fill = "gray90", alpha = 0.1, color = NA) +
  geom_line(aes(group = iter), alpha = 0.02) +
  geom_line(data = filter(risk_boot_mean, 
                          season == "early",
                          predator == "cougar"),
            color = "#993333",
            linewidth = 1) +
  geom_hline(yintercept = 0, color = "red", linetype = "dashed") +
  xlab(NULL) +
  ylab("log-RSS Risk") +
  scale_x_datetime(date_labels = "%H:%M") +
  coord_cartesian(expand = FALSE) +
  # coord_cartesian(ylim = c(-0.4, 0.4), expand = FALSE) +
  # ggtitle("Cougar, Early Winter") +
  theme(axis.text.x = element_blank())

# Late winter
risk_plot_late_cougar <- risk_pop_mean %>% 
  filter(season == "late", predator == "cougar") %>% 
  ggplot(mapping = aes(x = start_time, y = log_rss)) +
  # Left-side night
  geom_rect(aes(xmin = as.POSIXct("2022-03-15 00:00:00", tz = "US/Mountain"), 
                xmax = as.POSIXct("2022-03-15 07:35:00", tz = "US/Mountain"), 
                ymin = -Inf, 
                ymax = Inf),
            fill = "gray90", alpha = 0.1, color = NA) +
  # Right-side night
  geom_rect(aes(xmin = as.POSIXct("2022-03-15 19:29:00", tz = "US/Mountain"), 
                xmax = as.POSIXct("2022-03-15 23:59:59", tz = "US/Mountain"), 
                ymin = -Inf, 
                ymax = Inf),
            fill = "gray90", alpha = 0.1, color = NA) +
  geom_line(aes(group = iter), alpha = 0.02) +
  geom_line(data = filter(risk_boot_mean, 
                          season == "late",
                          predator == "cougar"),
            color = "#993333",
            linewidth = 1) +
  geom_hline(yintercept = 0, color = "red", linetype = "dashed") +
  xlab("Step Start Time") +
  ylab("log-RSS Risk") +
  scale_x_datetime(date_labels = "%H:%M") +
  coord_cartesian(expand = FALSE) +
  # coord_cartesian(ylim = c(-0.4, 0.4), expand = FALSE) +
  # ggtitle("Cougar, Late Winter")
  NULL

# Combine
risk_plot <- risk_plot_early_wolf +
  risk_plot_early_cougar + 
  risk_plot_late_wolf  + 
  risk_plot_late_cougar +
  plot_layout(nrow = 2, ncol = 2) +
  plot_annotation(tag_levels = "A", tag_prefix = "(", tag_suffix = ")")

# ggsave("fig/risk_by_time.tif", plot = risk_plot, device = agg_tiff, 
#        width = 9, height = 7, units = "in", compression = "lzw")

# Fig. 3C ----
# Selection for terrain
terr_dat <- beta_means %>% 
  filter(term %in% c("open_end", "rough_end"))

terr_plot_early <- terr_dat %>% 
  filter(season == "early") %>% 
  ggplot(aes(x = term, y = mean_est, group = term)) +
  geom_violin() +
  geom_hline(yintercept = 0, color = "red", linetype = "dashed") +
  xlab(NULL) +
  ylab("log-RSS Terrain") +
  # ggtitle("Early Winter") +
  coord_cartesian(ylim = c(-0.5, 0.6)) +
  theme(axis.text.x = element_blank())

terr_plot_late <- terr_dat %>% 
  filter(season == "late") %>% 
  ggplot(aes(x = term, y = mean_est, group = term)) +
  geom_violin() +
  geom_hline(yintercept = 0, color = "red", linetype = "dashed") +
  xlab("Terrain") +
  ylab("log-RSS Terrain") +
  # ggtitle("Late Winter") +
  scale_x_discrete(breaks = c("open_end", "rough_end"),
                   labels = c("Openness", "Roughness")) +
  coord_cartesian(ylim = c(-0.5, 0.6)) +
  theme()

# Combine
terr_plot <- terr_plot_early +
  terr_plot_late +
  plot_layout(nrow = 2, ncol = 1) +
  plot_annotation(tag_levels = "A", tag_prefix = "(", tag_suffix = ")")

ggsave("fig/terrain_by_time.tif", plot = terr_plot, device = agg_tiff, 
       width = 6.5, height = 4, units = "in", compression = "lzw")

# Full Figure 3 ----
## Row labels
# Row 1
row1_title <- ggplot() +
  geom_text(aes(x = 0, y = 0, label = "Early Winter"),
            size = 7, angle = 90) +
  theme_void()

# Row 2
row2_title <- ggplot() +
  geom_text(aes(x = 0, y = 0, label = "Late Winter"),
            size = 7, angle = 90) +
  theme_void()

row_labels <- wrap_elements(full = row1_title + 
                              row2_title +
                              plot_layout(ncol = 1, nrow = 2))

col1 <- wrap_elements(
  full = 
    bio_plot_early + bio_plot_late +
    plot_layout(ncol = 1, nrow = 2,
                heights = c(0.45, 0.55), tag_level = "keep") +
    plot_annotation(title = "Biomass",
                    tag_levels = list(c("A", "B")),
                    tag_prefix = "(",
                    tag_suffix = ")",
                    theme = theme(
                      plot.title = element_text(size = 18, 
                                                hjust = 0.7))) & 
    theme(plot.tag = element_text(face = "bold"))
)

col2 <- wrap_elements(
  full = 
    risk_plot_early_wolf + risk_plot_late_wolf +
    plot_layout(ncol = 1, nrow = 2,
                heights = c(0.45, 0.55), tag_level = "keep") +
    plot_annotation(title = "Wolf Risk",
                    tag_levels = list(c("C", "D")),
                    tag_prefix = "(",
                    tag_suffix = ")",
                    theme = theme(
                      plot.title = element_text(size = 18, 
                                                hjust = 0.7))) & 
    theme(plot.tag = element_text(face = "bold"))
)

col3 <- wrap_elements(
  full = 
    risk_plot_early_cougar + risk_plot_late_cougar +
    plot_layout(ncol = 1, nrow = 2,
                heights = c(0.45, 0.55), tag_level = "keep") +
    plot_annotation(title = "Cougar Risk",
                    tag_levels = list(c("E", "F")),
                    tag_prefix = "(",
                    tag_suffix = ")",
                    theme = theme(
                      plot.title = element_text(size = 18, 
                                                hjust = 0.7))) & 
    theme(plot.tag = element_text(face = "bold"))
)

col4 <- wrap_elements(
  full = 
    terr_plot_early + terr_plot_late +
    plot_layout(ncol = 1, nrow = 2,
                heights = c(0.45, 0.55), tag_level = "keep") +
    plot_annotation(title = "Terrain",
                    tag_levels = list(c("G", "H")),
                    tag_prefix = "(",
                    tag_suffix = ")",
                    theme = theme(
                      plot.title = element_text(size = 18, 
                                                hjust = 0.7))) & 
    theme(plot.tag = element_text(face = "bold"))
)


fig3 <- row_labels + col1 + col2 + col3 + col4 +
  plot_layout(nrow = 1, ncol = 5,
              widths = c(0.04, 0.24, 0.24, 0.24, 0.24))

ggsave("fig/fig3_issa.tif", plot = fig3, device = agg_tiff, 
       width = 12, height = 6, units = "in", dpi = 600, compression = "lzw")

# Movement coefficients ----

# ... step-length by ta ----

# Predict for early winter using December 1
ta_dat_early <- create_pred_dat(scale_df, 
                                ta_ = seq(-pi, pi, length.out = 200),
                                start_time = lubridate::ymd_hms(
                                  "2021-12-01 07:00:00", tz = "US/Mountain")) %>% 
  scale_dat(scale_df)

# Predict for late winter using March 15
ta_dat_late <- create_pred_dat(scale_df, 
                               ta_ = seq(-pi, pi, length.out = 200),
                               start_time = lubridate::ymd_hms(
                                 "2022-03-15 07:00:00", tz = "US/Mountain")) %>% 
  scale_dat(scale_df)

# Predict turn angle for all models
ta_pred_early <- predict_sl(dat, niter = 2000, newdata = ta_dat_early)
ta_pred_late <- predict_sl(dat, niter = 2000, newdata = ta_dat_late)

# ... ... population mean ----
ta_pop_mean_early <- ta_pred_early %>% 
  group_by(iter, ta_) %>% 
  summarize(mean_sl = mean(mean_sl, na.rm = TRUE))

ta_pop_mean_late <- ta_pred_late %>% 
  group_by(iter, ta_) %>% 
  summarize(mean_sl = mean(mean_sl, na.rm = TRUE))

# ... ... bootstrap mean ----
ta_boot_mean_early <- ta_pred_early %>% 
  group_by(ta_) %>% 
  summarize(mean_sl = mean(mean_sl, na.rm = TRUE))


ta_boot_mean_late <- ta_pred_late %>% 
  group_by(ta_) %>% 
  summarize(mean_sl = mean(mean_sl, na.rm = TRUE))

# ... ... plot ----
ta_plot_early <- ta_pop_mean_early %>% 
  ggplot(mapping = aes(x = ta_, y = mean_sl)) +
  geom_line(aes(group = iter), color = "black", alpha = 0.02) +
  geom_line(data = ta_boot_mean_early, color = "orange", linewidth = 1) +
  scale_x_continuous(name = "Turn Angle (radians)",
                     breaks = c(-pi, -pi/2, 0, pi/2, pi),
                     labels = expression(-pi, -pi/2, 0, pi/2, pi)) +
  ylab("Mean Step Length (m)") +
  coord_cartesian(ylim = c(0, 300))


ta_plot_late <- ta_pop_mean_late %>% 
  ggplot(mapping = aes(x = ta_, y = mean_sl)) +
  geom_line(aes(group = iter), color = "black", alpha = 0.02) +
  geom_line(data = ta_boot_mean_late, color = "orange", linewidth = 1) +
  scale_x_continuous(name = "Turn Angle (radians)",
                     breaks = c(-pi, -pi/2, 0, pi/2, pi),
                     labels = expression(-pi, -pi/2, 0, pi/2, pi)) +
  ylab("Mean Step Length (m)") +
  coord_cartesian(ylim = c(0, 300))

# Combine
ta_plot <- ta_plot_early + ta_plot_late +
  plot_annotation(tag_levels = "A", tag_prefix = "(", tag_suffix = ")")

ggsave("fig/sl_by_ta.tif", plot = ta_plot, device = agg_tiff, 
       width = 8, height = 5, units = "in", dpi = 600, compression = "lzw")

# ... step-length by SWE ----

# SWE in early winter
dat2 %>% 
  filter(season == "early") %>% 
  pull(swe_start_orig) %>% 
  hist()
# Values from 0 to 80 are reasonable

# SWE in late winter
dat2 %>% 
  filter(season == "late") %>% 
  pull(swe_start_orig) %>% 
  hist()
# Values from 0 to 300 are reasonable


# Predict for early winter using December 1
swe_dat_early <- create_pred_dat(scale_df, 
                                 swe_start = seq(0, 80, length.out = 100),
                                 start_time = lubridate::ymd_hms(
                                   "2021-12-01 07:00:00", tz = "US/Mountain")) %>% 
  scale_dat(scale_df)

# Predict for late winter using March 15
swe_dat_late <- create_pred_dat(scale_df, 
                                swe_start = seq(0, 300, length.out = 200),
                                start_time = lubridate::ymd_hms(
                                  "2022-03-15 07:00:00", tz = "US/Mountain")) %>% 
  scale_dat(scale_df)

# Predict turn angle for all models
swe_pred_early <- predict_sl(dat, niter = 2000, newdata = swe_dat_early)
swe_pred_late <- predict_sl(dat, niter = 2000, newdata = swe_dat_late)

# ... ... population mean ----
swe_pop_mean_early <- swe_pred_early %>% 
  group_by(iter, swe_start_orig) %>% 
  summarize(mean_sl = mean(mean_sl, na.rm = TRUE))

swe_pop_mean_late <- swe_pred_late %>% 
  group_by(iter, swe_start_orig) %>% 
  summarize(mean_sl = mean(mean_sl, na.rm = TRUE))

# ... ... bootstrap mean ----
swe_boot_mean_early <- swe_pred_early %>% 
  group_by(swe_start_orig) %>% 
  summarize(mean_sl = mean(mean_sl, na.rm = TRUE))


swe_boot_mean_late <- swe_pred_late %>% 
  group_by(swe_start_orig) %>% 
  summarize(mean_sl = mean(mean_sl, na.rm = TRUE))

# ... ... plot ----
swe_plot_early <- swe_pop_mean_early %>% 
  ggplot(mapping = aes(x = swe_start_orig, y = mean_sl)) +
  geom_line(aes(group = iter), color = "black", alpha = 0.02) +
  geom_line(data = swe_boot_mean_early, color = "#80CCFF", linewidth = 1) +
  xlab(expression("Snow-Water Equivalent" ~ (kg/m^2))) +
  ylab("Mean Step Length (m)") +
  coord_cartesian(xlim = c(0, 300), 
                  ylim = c(0, 300))


swe_plot_late <- swe_pop_mean_late %>% 
  ggplot(mapping = aes(x = swe_start_orig, y = mean_sl)) +
  geom_line(aes(group = iter), color = "black", alpha = 0.02) +
  geom_line(data = swe_boot_mean_late, color = "#80CCFF", linewidth = 1) +
  xlab(expression("Snow-Water Equivalent" ~ (kg/m^2))) +
  ylab("Mean Step Length (m)") +
  coord_cartesian(xlim = c(0, 300), 
                  ylim = c(0, 300))

# Combine
swe_plot <- swe_plot_early + swe_plot_late +
  plot_annotation(tag_levels = "A", tag_prefix = "(", tag_suffix = ")")

ggsave("fig/sl_by_swe.tif", plot = swe_plot, device = agg_tiff, 
       width = 8, height = 5, units = "in", dpi = 600, compression = "lzw")

# ... mean SL by time of day ----
# Predict for early winter using December 1
tod_dat_early <- create_pred_dat(scale_df, 
                                 start_time = seq(
                                   from = lubridate::ymd_hms(
                                     "2021-12-01 00:00:00", tz = "US/Mountain"),
                                   to = lubridate::ymd_hms(
                                     "2021-12-02 00:00:00", tz = "US/Mountain"),
                                   by = "1 hour"
                                 )) %>% 
  scale_dat(scale_df) %>% 
  mutate(hour = 0:24)

# Predict for late winter using March 15
tod_dat_late <- create_pred_dat(scale_df, 
                                start_time = seq(
                                  from = lubridate::ymd_hms(
                                    "2022-03-15 00:00:00", tz = "US/Mountain"),
                                  to = lubridate::ymd_hms(
                                    "2022-03-16 00:00:00", tz = "US/Mountain"),
                                  by = "1 hour"
                                )) %>% 
  scale_dat(scale_df) %>% 
  mutate(hour = 0:24)

# Predict turn angle for all models
tod_pred_early <- predict_sl(dat, niter = 2000, newdata = tod_dat_early)
tod_pred_late <- predict_sl(dat, niter = 2000, newdata = tod_dat_late)

# ... ... population mean ----
tod_pop_mean_early <- tod_pred_early %>% 
  group_by(iter, hour) %>% 
  summarize(mean_sl = mean(mean_sl, na.rm = TRUE))

tod_pop_mean_late <- tod_pred_late %>% 
  group_by(iter, hour) %>% 
  summarize(mean_sl = mean(mean_sl, na.rm = TRUE))

# ... ... bootstrap mean ----
tod_boot_mean_early <- tod_pred_early %>% 
  group_by(hour) %>% 
  summarize(mean_sl = mean(mean_sl, na.rm = TRUE))


tod_boot_mean_late <- tod_pred_late %>% 
  group_by(hour) %>% 
  summarize(mean_sl = mean(mean_sl, na.rm = TRUE))

# ... ... plot ----
tod_plot_early <- tod_pop_mean_early %>% 
  ggplot(mapping = aes(x = hour, y = mean_sl)) +
  geom_line(aes(group = iter), color = "black", alpha = 0.02) +
  geom_line(data = tod_boot_mean_early, color = "hotpink", linewidth = 1) +
  scale_x_continuous(name = "Hour",
                     breaks = seq(0, 24, by = 4),
                     labels = c("00:00", "04:00", "08:00", "12:00",
                                "16:00", "20:00", "00:00")) +
  ylab("Mean Step Length (m)") +
  coord_cartesian(ylim = c(0, 300))


tod_plot_late <- tod_pop_mean_late %>% 
  ggplot(mapping = aes(x = hour, y = mean_sl)) +
  geom_line(aes(group = iter), color = "black", alpha = 0.02) +
  geom_line(data = tod_boot_mean_late, color = "hotpink", linewidth = 1) +
  scale_x_continuous(name = "Hour",
                     breaks = seq(0, 24, by = 4),
                     labels = c("00:00", "04:00", "08:00", "12:00",
                                "16:00", "20:00", "00:00")) +
  ylab("Mean Step Length (m)") +
  coord_cartesian(ylim = c(0, 300))

# Combine
tod_plot <- tod_plot_early + tod_plot_late +
  plot_annotation(tag_levels = "A", tag_prefix = "(", tag_suffix = ")")

ggsave("fig/sl_by_tod.tif", plot = tod_plot, device = agg_tiff, 
       width = 8, height = 5, units = "in", dpi = 600, compression = "lzw")

# ... one movement figure ----
mov_plot <- (ta_plot_early + ggtitle("Early Winter")) + 
  (ta_plot_late + ggtitle("Late Winter")) +
  swe_plot_early + swe_plot_late +
  tod_plot_early + tod_plot_late +
  plot_layout(ncol = 2, nrow = 3) +
  plot_annotation(tag_levels = "A",
                  tag_prefix = "(",
                  tag_suffix = ")") &
  theme(plot.title = element_text(hjust = 0.5, size = 20))

ggsave("fig/mov_fig.tif", plot = mov_plot, device = agg_tiff, 
       width = 7, height = 6, units = "in", dpi = 600, compression = "lzw")

# Correlations ----
b_wide <- b_df %>% 
  dplyr::select(ID, winter, season, n_steps, term, estimate) %>% 
  pivot_wider(names_from = term, values_from = estimate)

b_wide %>% 
  select(bio_end:`cougar_risk_end:cougar_act_end`) %>% 
  cor() %>% 
  image(col = hcl.colors(12, palette = "viridis"))

{agg_tiff("fig/issa_coefs_corr.tif", res = 600,
          width = 12, height = 10, units = "in", compression = "lzw")
  b_wide %>% 
    select(bio_end:`cougar_risk_end:cougar_act_end`) %>% 
    plot()
  dev.off()
  }
