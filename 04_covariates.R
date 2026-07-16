############################################X
#--------------Elk Risk Taking--------------X
#---------------Brian J. Smith--------------X
#-------------Attach Covariates-------------X
############################################X

# Set options ----
# For full printing of tibble columns
options(pillar.width = 1000)

# Load packages ----
library(dplyr)
library(lubridate)
library(purrr)
library(amt)
library(sf)
library(terra)

# Source functions
source("99_fun.R")

# Load data ----
dat <- readRDS("../out/steps_df.rds") %>% 
  # clean_season points outside of season have already been filtered
  rename(season = clean_season)

# Attach covariates ----
# Add elk ID
dat$rand <- mapply(add_id, x = dat$rand, id = dat$ID, SIMPLIFY = FALSE)
# Add winter
dat$rand <- mapply(add_winter, x = dat$rand, winter = dat$winter, 
                   SIMPLIFY = FALSE)

# Takes a while...
# Actual time was 16650 s = 16650/60/60 hr
system.time({
  dat <- dat %>%
    mutate(covs = map(rand, function(x) {
      x %>% 
        attach_snow(verbose = TRUE) %>%
        attach_RAP(verbose = TRUE) %>%
        attach_DEM(verbose = TRUE) %>%
        attach_hour(verbose = TRUE) %>%
        attach_risk(verbose = TRUE) %>% 
        # Take log(swe)
        mutate(log_swe_start = log(swe_start + 1),
               log_swe_end = log(swe_end + 1)) %>% 
        return()
    }))
})

# Save ----
dat <- dat %>% 
  dplyr::select(-rand)

saveRDS(dat, "../out/steps_w_covs.rds")


if(!interactive()){
  message("\n========================================\nFinished at ", Sys.time())
} else {
  beepr::beep(11)
}