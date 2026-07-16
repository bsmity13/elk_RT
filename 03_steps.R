############################################X
#--------------Elk Risk Taking--------------X
#---------------Brian J. Smith--------------X
#----------------Format Steps---------------X
############################################X

# Set options ----
# For full printing of tibble columns
options(pillar.width = 1000)

# Load packages ----
library(dplyr)
library(lubridate)
library(purrr)
library(ggplot2)
library(ragg)
library(amt)
library(sf)
library(terra)

# Source functions
source("99_fun.R")

# Load data ----
# Yellowstone boundary (rasterized)
ynp <- rast("../elk_RT_data/geo/ynp.tif")

# 1h GPS data
dat <- readRDS("../out/GPS_1h.rds")

# Format as steps ----
stp <- dat %>% 
  select(ID, winter, clean_season, burst_, season, x_, y_, t_) %>% 
  nest(data = burst_:t_) %>% 
  mutate(steps = map(data, function(d) {
    d %>% 
      steps_by_burst(zero_dir = "N",
                     clockwise = TRUE,
                     lonlat = FALSE,
                     keep_cols = "start") %>% 
      # Keep only those steps within a season
      filter(!is.na(season)) %>% 
      # Get rid of 0-length steps
      filter(sl_ != 0) %>% 
      # Flag NA turn angles
      mutate(flag = is.na(ta_)) %>% 
      return()
  }),
  n_steps = map_dbl(steps, function(s) {
    return(nrow(s))
  })) %>% 
  # Require at least 20 steps
  filter(n_steps >= 20) %>% 
  # Fit SL and TA distributions
  # Fit step length distributions
  mutate(gamma_dist = map(steps, function(s) {
    return(fit_distr(s$sl_, "gamma")) # Warnings here seem to be fine
  }),
  # Fit turn angle distribution (von Mises)
  ta_dist = map(steps, function(s) {
    return(fit_distr(s$ta_, "vonmises"))
  }),
  # Gamma shape parameter
  shape = map_dbl(gamma_dist, function(x) {
    return(x$params$shape)
  }),
  # Gamma scale parameter
  scale = map_dbl(gamma_dist, function(x) {
    return(x$params$scale)
  }),
  # von Mises concentration parameter
  kappa = map_dbl(ta_dist, function(x) {
    return(x$params$kappa)
  }),
  # Start and end times
  start = map_chr(steps, function(x) {
    return(format(min(x$t1_, na.rm = TRUE), "%Y-%m-%d %H:%M:%S"))
  }) %>% 
    ymd_hms(tz = "US/Mountain"),
  end = map_chr(steps, function(x) {
    return(format(max(x$t2_, na.rm = TRUE), "%Y-%m-%d %H:%M:%S"))
  }) %>% 
    ymd_hms(tz = "US/Mountain"))

# Population-level tentative distributions ----
## For step length, choose large values to sample enough long steps
# Shape
shp <- quantile(stp$shape, 0.90)
hist(stp$shape, breaks = 100)
abline(v = shp, col = "red")

# Scale
scl <- quantile(stp$scale, 0.90)
hist(stp$scale, breaks = 100)
abline(v = scl, col = "red")

# Average step length
unname(shp * scl)

## Don't need to exaggerate turn angle like step length -- more diffuse is 
# probably better
# Kappa
k <- median(stp$kappa)
hist(stp$kappa, breaks = 100)
abline(v = k, col = "red")

# ... tentative step-length distribution ----
tent_sl <- make_gamma_distr(shape = shp, scale = scl)

tent_sl_fig <- data.frame(x = seq(0, 5000, length.out = 300)) %>% 
  mutate(y = dgamma(x, shape = shp, scale = scl)) %>% 
  ggplot(aes(x = x, y = y)) +
  geom_line() +
  xlab("Step Length (m)") +
  ylab("Probability Density") +
  ggtitle("Tentative Step-length Distribution") +
  theme_bw() +
  NULL

ggsave("fig/tent_sl.tif", plot = tent_sl_fig, device = agg_tiff,
       width = 6, height = 4, units = "in", compression = "lzw")


# ... tentative turn-angle distribution ----
tent_ta <- make_vonmises_distr(kappa = k)

tent_ta_fig <- data.frame(x = seq(-pi, pi, length.out = 100)) %>% 
  mutate(y = suppressWarnings({
    circular::dvonmises(x, mu = 0, kappa = k)
  })) %>% 
  ggplot(aes(x = x, y = y)) +
  geom_line() +
  xlab("Turn Angle (radians)") +
  ylab("Probability Density") +
  scale_x_continuous(breaks = c(-pi, -pi/2, 0, pi/2, pi),
                     labels = expression(-pi, -pi/2, 0, pi/2, pi)) +
  coord_cartesian(ylim = c(0, 0.25)) +
  ggtitle("Tentative Turn-angle Distribution") +
  theme_bw() +
  NULL

ggsave("fig/tent_ta.tif", plot = tent_ta_fig, device = agg_tiff,
       width = 6, height = 4, units = "in", compression = "lzw")

# Remove steps outside YNP ----
# Removing steps that start or end outside of YNP
# Important to do this now because any elk near the boundary will have
# some steps in and some steps out, but we wanted to create the steps (above)
# on the continuous trajectory before subsampling.
stp2 <- stp %>% 
  mutate(steps = map(steps, function(x) {
    x %>% 
      amt::extract_covariates(ynp, where = "both") %>% 
      filter(ynp_start == 1 & ynp_end == 1) %>% 
      return()
  }),
  n_steps = map_dbl(steps, nrow)) %>% 
  # None of these chosen individuals winter outside of YNP, so this shouldn't
  # be relevant, but just in case
  filter(n_steps >= 100) %>% 
  mutate(
  start = map_vec(steps, function(x) min(x$t1_)),
  end = map_vec(steps, function(x) max(x$t2_))
  )

c("full_data" = nrow(stp),
  "inside_YNP" = nrow(stp2))

# Sample available steps ----
set.seed(20220726)

stp2 <- stp2 %>% 
  mutate(rand = map(steps, function(x) {
    x %>% 
      random_steps(n_control = 50,
                   sl_distr = tent_sl,
                   ta_distr = tent_ta) %>% 
      # Drop flagged steps (no TA)
      filter(!flag) %>% 
      select(-flag) %>% 
      # Combine step_id_ and burst_
      mutate(step_id_ = paste(burst_, step_id_, sep = "_")) %>% 
      return()
  }),
  tot_steps = map_dbl(rand, nrow))

# Check available step spatial distribution
# {
#   par(mar = rep(0, 4))
#   plot(nr$geometry)
#   xx <- stp %>% 
#     select(ID, winter, rand) %>% 
#     unnest(cols = rand)
#   points(xx$x2_, xx$y2_, pch = ".", col = "red")
#   }

# Save ----
# Don't save redundant data columns
res <- stp2 %>% 
  dplyr::select(ID, winter, clean_season, start, end, shape, scale, kappa, 
                n_steps, tot_steps, rand)

saveRDS(res, "../out/steps_df.rds")


if(!interactive()){
  message("\n========================================\nFinished at ", Sys.time())
} else {
  beepr::beep(11)
}