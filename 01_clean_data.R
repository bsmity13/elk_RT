############################################X
#--------------Elk Risk Taking--------------X
#---------------Brian J. Smith--------------X
#-------------GPS Data Cleaning-------------X
############################################X

# Load packages ----
library(amt)
library(dplyr)
library(tidyr)
library(purrr)
library(lubridate)

# Load data ----
# GPS data
raw <- read.csv("../elk_RT_data/raw_GPS.csv") %>% 
  mutate(ID = as.character(ID),
         dt = ymd_hms(dt, tz = "US/Mountain")) %>% 
  rename(x = utm_e, y = utm_n, t = dt) %>% 
  # Multiply number of GPS satellites by -1
  # so that lower value is always better precision
  # (whether DOP or Number of Sats)
  mutate(precision = case_when(
    precision_type == "GPS Satellite Count" ~ precision * -1,
    TRUE ~ precision
  )) %>% 
  # Some collars went to 5 or 10 minute fixes during aerial survey
  # Subset to the hour to get rid of these
  filter(minute(t) %in% c(58:59, 0:2))

# Fates
# Just need deaths that occurred within a season to filter out the remaining
# data
fate <- read.csv("../elk_RT_data/fate.csv") %>% 
  mutate(last_loc = as.POSIXct(last_loc)) %>% 
  mutate(ID = as.character(ID),
         year = year(last_loc),
         month = month(last_loc)) %>% 
  filter(month %in% c(12, 1:3)) %>% 
  mutate(winter = case_when(
    month > 11 ~ year + 1,
    TRUE ~ year
  )) %>% 
  mutate(cause2 = case_when(
    cause %in% c("Malnutrition", "Natural-Unknown", "Natural-Other") ~ "Natural",
    TRUE ~ cause
  )) %>% 
  select(ID, winter, last_loc, cause, cause2)

# Save fates ----
saveRDS(fate, "../out/fate_formatted.rds")

# Attach fates to GPS data ----
raw2 <- raw %>% 
  left_join(fate, by = c("ID", "winter")) %>% 
  # If last_loc is NA, set to final date of tracking
  group_by(ID) %>% 
  mutate(last_gps = max(t)) %>% 
  ungroup() %>% 
  mutate(last_loc = case_when(
    is.na(last_loc) ~ last_gps,
    TRUE ~ last_loc
  )) %>% 
  # Drop unneeded column
  select(-last_gps) %>% 
  # Keep only GPS locations before the last_loc
  filter(t <= last_loc)

# Separate into winter-seasons ----
# Separating into winters (named for that January) and seasons (early or late).
# Will keep a couple of days on either side of each season to facilitate
# calculations, but will ultimately only keep correct season dates:
#   Correct season dates:
#     Early Winter = Nov 15 -- Dec 15
#     Late Winter  = Mar 01 -- Mar 31
#   Season dates for cleaning:
#     Early Winter = Nov 12 -- Dec 18
#     Late Winter  = Feb 26 -- Apr 03

raw2 <- raw2 %>% 
  # Correct winter-seasons
  mutate(winter = case_when(
    month(t) < 6 ~ year(t),
    month(t) > 6 ~ year(t) + 1
  ),
  season = case_when(
    (month(t) == 11 & day(t) >= 15) |
      (month(t) == 12 & day(t) <= 15) ~ "early",
    month(t) == 3 ~ "late"
  ),
  winter_season = case_when(
    is.na(season) ~ NA_character_,
    TRUE ~ paste(winter, season, sep = "_")
    )) %>% 
  # Cleaning winter-seasons
  mutate(clean_season = case_when(
    (month(t) == 11 & day(t) >= 12) |
      (month(t) == 12 & day(t) <= 18) ~ "early",
    (month(t) == 3) |
      (month(t) == 2 & day(t) >= 26) |
      (month(t) == 4 & day(t) <= 3) ~ "late"
  ),
  clean_winter_season = paste(winter, clean_season, sep = "_"))

# Subset to just the clean_seasons
raw3 <- raw2 %>% 
  filter(!is.na(clean_season))

# Take a quick look
# View(head(raw3, 100))

## Data cleaning

# Locations without precision ----
table(raw3$precision_type, useNA = "always")
# How many NAs?
sum(is.na(raw3$precision))/nrow(raw3)
# None

# Assign unique row ID ----
raw3 <- raw3 %>% 
  arrange(ID, t, precision) %>% 
  mutate(row = 1:nrow(.)) %>% 
  dplyr::select(row, everything())

# Convert to nested data.frame ----
ndf <- raw3 %>% 
  nest(dat = -c(ID, clean_winter_season)) %>% 
  mutate(trk1 = map(dat, function(dd) {
    dd %>% 
      make_track(x, y, t, crs = 32612, all_cols = TRUE)
  }))

# List of rows to flag ----
flag <- list()

# Remove poor-quality fixes ----
# What do poor-quality fixes look like?
# Note: original data had collars that could have DOP or number of satellites
# as the measure of precision. This subset has only DOP.
if (interactive()){
  ndf$trk1 %>% 
    bind_rows() %>% 
    group_by(precision_type) %>% 
    summarize(min = min(precision, na.rm = TRUE),
              med = median(precision, na.rm = TRUE),
              mean = mean(precision, na.rm = TRUE),
              max = max(precision, na.rm = TRUE))
  
  # DOP
  ndf$trk1 %>% 
    bind_rows() %>% 
    filter(precision_type != "GPS Satellite Count") %>% 
    pull(precision) %>% 
    quantile(c(0.9, 0.95, 0.99))
  
  ndf$trk1 %>% 
    bind_rows() %>% 
    filter(precision_type != "GPS Satellite Count") %>% 
    filter(precision < 10) %>% 
    pull(precision) %>% 
    hist(main = "Frequency of DOPs")
}

# Flag any DOP > 6 and any Satellites > (-)5
ndf <- ndf %>% 
  mutate(trk2 = map(trk1, function(tt) {
    tt %>% 
      mutate(low_prec = case_when(
        precision_type == "GPS Satellite Count" ~ precision > -5,
        precision_type %in% c("DOP", "HDOP") ~ precision > 6
      ))
  }))

flag$low_prec <- ndf$trk2 %>% 
  bind_rows() %>% 
  filter(low_prec) %>% 
  pull(row)

# Remove low-quality duplicates ----
# This is slow
system.time({
  ndf <- ndf %>% 
    mutate(trk3 = map(trk2, function(tt) {
      tt %>% 
        filter(!low_prec) %>% 
        arrange(t_) %>% 
        flag_duplicates(gamma = 5, time_unit = "mins", DOP = "precision")
    }))
})


print(ndf$trk3[[1]], n = 5, width = 200)

flag$duplicate_ <- ndf$trk3 %>% 
  bind_rows() %>% 
  filter(duplicate_) %>% 
  pull(row)
  
# Remove unreasonably fast steps ----
# SDR
# Say an elk can go 50 km/hr for 5 minutes
calculate_sdr(50, time = minutes(5))
# Use 5e4 as the cutoff
sdr_max <- 5e4

# This is fast
system.time({
  ndf <- ndf %>% 
    mutate(trk4 = map(trk3, function(tt) {
      tt %>% 
        filter(!duplicate_) %>% 
        arrange(t_) %>% 
        flag_fast_steps(delta = sdr_max)
    }))
})

print(ndf$trk4[[1]], n = 5, width = 200)

flag$fast_step_ <- ndf$trk4 %>% 
  bind_rows() %>% 
  filter(fast_step_) %>% 
  pull(row)

length(flag$fast_step_)

# Remove unreasonably fast round trips ----
# This is medium
system.time({
  ndf <- ndf %>% 
    mutate(trk5 = map(trk4, function(tt) {
      tt %>% 
        filter(!fast_step_) %>% 
        arrange(t_) %>% 
        flag_roundtrips(delta = sdr_max, epsilon = 10)
    }))
})

print(ndf$trk5[[1]], n = 5, width = 200)

flag$fast_roundtrip_ <- ndf$trk5 %>% 
  bind_rows() %>% 
  filter(fast_roundtrip_) %>% 
  pull(row)

length(flag$fast_roundtrip_)

# Clean data ----
ndf <- ndf %>% 
  mutate(clean = map(trk5, function(tt) {
    tt %>% 
      filter(!fast_roundtrip_) %>% 
      arrange(t_)
  })) %>% 
  mutate(n_orig = map_dbl(trk1, nrow),
         n_clean = map_dbl(clean, nrow),
         n_diff = n_orig - n_clean,
         perc_diff = n_diff/n_orig)

clean_df <- ndf %>% 
  dplyr::select(ID, clean:perc_diff)

clean <- clean_df %>% 
  dplyr::select(ID, clean) %>% 
  unnest(cols = clean)

# Save ----
saveRDS(clean, "../out/clean_GPS.rds")
saveRDS(flag, "../out/flagged_rows.rds")

if(!interactive()){
  message("\n========================================\nFinished at ", Sys.time())
} else {
  beepr::beep(11)
}
