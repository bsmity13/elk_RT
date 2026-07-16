############################################X
#--------------Elk Risk Taking--------------X
#---------------Brian J. Smith--------------X
#----------------Import Data----------------X
############################################X

# Loads data from BJS' computer and saves for use in this project
# This script should not be run by users of the GitHub/Zenodo archive

# Load packages ----
library(dplyr)
library(lubridate)

# Individuals to keep ----
# Hourly data during early/late winter
inds <- read.csv("../elk_RT_data/individuals_to_use.csv") %>% 
  mutate(start = ymd_hms(start, tz = "US/Mountain"),
         end = ymd_hms(end, tz = "US/Mountain"))

# GPS ----
# Load combined GPS data
d <- "../../../Data/Elk Data/Elk GPS/Elk GPS Proc/dat"
gps_all <- readRDS(file.path(d, "combined.rds"))

# Subset to individuals and dates
# Early winter: November 15 to December 15
# Late winter: March 1 - March 31
gps <- lapply(1:nrow(inds), function(i) {
  r <- inds[i, ]
  res <- gps_all %>% 
    filter(ID == r$ID,
           # Keep a couple of days on either side to help with
           # data cleaning and creating steps with turn angles.
           dt >= (r$start - days(2)) & dt <= (r$end + days(2))) %>% 
    mutate(winter = r$winter, season = r$season)
}) %>% 
  bind_rows() %>% 
  # Manually convert timestamp to character before saving to CSV
  mutate(dt = format(dt, "%Y-%m-%d %H:%M:%S"))

# Save CSV of GPS data
write.csv(gps, "../elk_RT_data/raw_GPS.csv", row.names = FALSE)

# GPS metadata ----
dmeta <- "../../../Data/Elk Data"
meta_all <- read.csv(file.path(dmeta, "NR_Elk_Study_Capture_ID.csv"))

# Select columns
meta <- meta_all %>% 
  select(ID, cap_date = Capture.Date, 
         cementum_age = Cementum.Age,
         birth = BirthDate, birth_year = BirthYear) %>% 
  # Format columns
  mutate(ID = as.character(ID),
         cap_date = mdy_hms(cap_date, tz = "US/Mountain"),
         birth = mdy_hms(birth, tz = "US/Mountain")) %>% 
  # Keep only relevant individuals
  filter(ID %in% inds$ID)

# Check birth year field
all.equal(meta$birth_year, year(meta$birth))

# Save CSV of capture metadata
write.csv(meta, "../elk_RT_data/capture_metadata.csv", row.names = FALSE)

# Mortality ----
# Load mortality data
dmort <- "../../../Data/Elk Data/"
mort <- read.csv(file.path(dmort, "NR-Elk_Mortality_4-18-22.csv"))
unk <- read.csv(file.path(dmort, "NR-Elk_Unknown_Fate_4-18-22.csv"))

# Select columns
mort <- mort %>% 
  select(ID, last_loc = Date.Last.Location,
         mort_utm_e = UTM.Easting, mort_utm_n = UTM.Northing, 
         cause = Cause.of.Death..COD., cause_certainty = COD.Certainty) %>% 
  mutate(last_loc = mdy_hms(last_loc))

unk <- unk %>% 
  select(ID, last_loc = Date.Last.Location,
         mort_utm_e = UTM.Easting, mort_utm_n = UTM.Northing, 
         cause = Cause.of.Death..COD., cause_certainty = COD.Certainty) %>% 
  mutate(last_loc = mdy_hms(last_loc),
         cause = "Unknown")

fate <- rbind(mort, unk) %>% 
  mutate(cause = case_when(
    cause == "" ~ "Unknown",
    cause == "Harvest Wound" ~ "Harvest",
    TRUE ~ cause),
    cause_certainty = case_when(
      cause_certainty == "" ~ "Unknown",
      TRUE ~ cause_certainty)
  ) %>% 
  # Keep only relevant individuals
  filter(ID %in% inds$ID)

# Save CSV of fate
write.csv(fate, "../elk_RT_data/fate.csv", row.names = FALSE)
