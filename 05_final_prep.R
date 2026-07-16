############################################X
#--------------Elk Risk Taking--------------X
#---------------Brian J. Smith--------------X
#-----------------Final Prep----------------X
############################################X

# Set options ----
# For full printing of tibble columns
options(pillar.width = 1000)

# Load packages ----
library(dplyr)
library(tidyr)
library(amt)

# Custom functions
source("99_fun.R")

# Load data ----
dat <- readRDS("out/steps_w_covs.rds")

# Scale and center covariates ----
# Based on mean and SD of entire dataset
df <- dat %>% 
  select(covs) %>% 
  unnest(cols = covs)

scale_df <- df %>% 
  select(swe_start:rough_end, 
         cougar_act_start:log_swe_end,
         -period, -solar_time) %>% 
  pivot_longer(everything(), names_to = "term") %>% 
  group_by(term) %>% 
  summarize(mean = mean(value, na.rm = TRUE),
            sd = sd(value, na.rm = TRUE))

# Save
write.csv(scale_df, "../out/scaling_data.csv", row.names = FALSE)

dat <- dat %>% 
  mutate(covs = map(covs, function(x) {
    x %>% 
      scale_dat(scale_df = scale_df) %>% 
      return()
  }))

# Save final dataset ----
saveRDS(dat, "../out/final_nested_data.rds")


if(!interactive()){
  message("\n========================================\nFinished at ", Sys.time())
} else {
  beepr::beep(11)
}