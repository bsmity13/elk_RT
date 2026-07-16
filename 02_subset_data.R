############################################X
#--------------Elk Risk Taking--------------X
#---------------Brian J. Smith--------------X
#----------------Data Subset----------------X
############################################X

# Making sure I have consistent, 1-h steps.
# Also visualizing the timeseries of data for all individuals.

dir.create("fig")

# Load packages ----
library(dplyr)
library(stringr)
library(amt)
library(lubridate)
library(ggplot2)
library(patchwork)
library(ragg)

# Source functions
source("99_fun.R")

# Save plots?
save_plots <- TRUE

# Load data ----
dat <- readRDS("../out/clean_GPS.rds")
fate <- readRDS("../out/fate_formatted.rds")

# Format IDs as factor ----
# Get all IDs
all_IDs <- sort(unique(c(dat$ID, fate$ID)))
# Pad with 0s for correct chronological sorting
all_IDs <- sort(str_pad(all_IDs, width = 4, side = "left", pad = "0"))
all_IDs <- factor(all_IDs)

dat <- dat %>% 
  mutate(ID = str_pad(ID, width = 4, side = "left", pad = "0")) %>% 
  mutate(ID = factor(ID, levels = levels(all_IDs)))

fate <- fate %>% 
  mutate(ID = str_pad(ID, width = 4, side = "left", pad = "0")) %>% 
  mutate(ID = factor(ID, levels = levels(all_IDs)))

# Existing fix rates ----
fr <- lapply(split(dat, as.character(dat$ID)), function(x){
  sr <- x %>% 
    make_track(x_, y_, t_, crs = 32612, all_cols = TRUE) %>% 
    summarize_sampling_rate(time_unit = "hour")
  summ <- x %>% 
    summarize(ID = unique(ID),
              start = min(t_),
              end = max(t_))
  res <- cbind(summ, sr)
  return(res)
}) %>% 
  bind_rows() %>% 
  mutate(fr = round(median))

table(fr$fr)

# Some potentially wrong start dates (re-used collars)
# First two digits of the ID (for all those > 1000) are the deployment year
# Although note they might have been collared in December of the previous year
fr <- fr %>% 
  mutate(collar_year = case_when(
    as.numeric(as.character(ID)) > 1000 ~ as.numeric(paste0(20, 
                                                            stringr::str_sub(ID, start = 1, end = 2))),
    TRUE ~ 2000)) %>% 
  mutate(collar_start = ymd_hm(paste0(collar_year - 1, 
                                      "-12-01 00:00"))) %>% 
  rowwise() %>% 
  mutate(start2 = max(c(start, collar_start))) %>% 
  as.data.frame()

st2 <- fr %>% 
  dplyr::select(ID, start2)

dat <- left_join(dat, st2, by = "ID") %>% 
  filter(t_ >= start2)

# Plot
fr_plot <- ggplot(fr, aes(x = start2, xend = end, 
                          y = factor(ID), yend = factor(ID),
                          color = factor(fr))) +
  geom_segment() +
  scale_color_viridis_d(name = "Fix Rate (h)") +
  xlab(NULL) +
  ylab("ID") +
  theme_bw() +
  theme(axis.text.y = element_text(size = 6)) +
  NULL

if(save_plots) {
  ggsave("fig/fixrate.tif", plot = fr_plot, width = 6, height = 12,
         units = "in", dpi = 500, compression = "lzw", device = agg_tiff)
}

# Subset to 1h in winter ----
# Kohl et al. 2019 used 1 Nov to 30 Apr
# Here, I only have those steps labeled with 'clean_season'
dat1 <- dat %>%
  mutate(id_clean_winter_season = paste(as.character(ID), 
                                        winter,
                                        clean_season,
                                        sep = "_")) %>% 
  split(.$id_clean_winter_season) %>% 
  lapply(function(id) {
    if(nrow(id) < 10) {
      return(NULL)
    } else {
      id %>% 
        make_track(x_, y_, t_, crs = 32612, all_cols = TRUE) %>% 
        # This function finds the starting time that results in
        # the maximum number of 1-hour fixes.
        subset_fixrate(rate = "1 hour") %>% 
        # Find the 1h steps
        track_resample(rate = hours(1), tolerance = minutes(5)) %>% 
        filter_min_n_burst(min_n = 3) %>% 
        return()
    }
  }) %>% 
  bind_rows()

# Plot
dat1_plot <- dat1 %>% 
  ggplot(aes(x = t_, y = ID)) +
  geom_point(size = 0.5) +
  xlab(NULL) +
  ylab("ID") +
  ggtitle("Winter Only, 1h Fixes") +
  theme_bw() +
  theme(axis.text.y = element_text(size = 6)) +
  NULL

if(save_plots) {
  ggsave("fig/fixes_winter_1h.tif", plot = dat1_plot, width = 6, height = 12,
         units = "in", dpi = 500, compression = "lzw", device = agg_tiff)
}

# How many ...
#   individuals?
length(unique(dat1$ID))
#   total points?
nrow(dat1)
#   average points per individual?
dat1 %>% 
  group_by(ID) %>% 
  tally() %>% 
  pull(n) %>% 
  mean()

#   average points per season?
dat1 %>% 
  group_by(id_clean_winter_season) %>% 
  tally() %>% 
  pull(n) %>% 
  mean()

# Add mortalities ----

# Instead, need to combine datasets
a <- dat1 %>% 
  select(ID, t_) %>% 
  mutate(type = "GPS")
b <- fate %>% 
  select(ID, t_ = last_loc) %>% 
  mutate(type = "Mortality")
dat1_mort <- rbind(a, b) %>% 
  ggplot(aes(x = t_, y = ID, color = type, shape = type)) +
  geom_point(size = 0.5) +
  xlab(NULL) +
  ylab("ID") +
  ggtitle("Winter Only, 1h Fixes") +
  scale_color_manual(breaks = c("GPS", "Mortality"), 
                     values = c("black", "red")) +
  
  scale_shape_manual(breaks = c("GPS", "Mortality"), 
                     values = c(16, 4)) +
  theme_bw() +
  theme(axis.text.y = element_text(size = 6)) +
  NULL

# Save data ----
# Revert ID back to unpadded character
dat1 <- dat1 %>% 
  mutate(ID = as.character(as.numeric(as.character(ID))))
saveRDS(dat1, "../out/GPS_1h.rds")


if(!interactive()){
  message("\n========================================\nFinished at ", Sys.time())
} else {
  beepr::beep(11)
}
