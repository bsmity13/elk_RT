############################################X
#--------------Elk Risk Taking--------------X
#---------------Brian J. Smith--------------X
#===========================================X
#------------------Fit iSSA-----------------X
############################################X

# Set options ----
# For full printing of tibble columns
options(pillar.width = 1000)

# Load packages ----
library(dplyr)
library(tidyr)
library(amt)
library(broom)

# Custom functions
source("99_fun.R")

# Load data ----
dat <- readRDS("../out/final_nested_data.rds")

# Question:
# How many individuals monitored for multiple years?
multi <- dat %>% 
  select(ID, winter, season) %>% 
  group_by(ID) %>% 
  summarize(yrs = n_distinct(winter)) %>% 
  arrange(desc(yrs))

table(multi$yrs)

# Fit models ----
# Varying shape parameter by including intxn with log(sl_)
# (Would vary scale parameter by including intxn with sl_)
mods <- dat %>% 
  mutate(issf = map(covs, function(x) {
    x %>% 
      fit_issf(case_ ~ 
                 
                 ## Resources ##
                 # Biomass
                 bio_end + 
                 bio_end:sin(solar_time) + 
                 bio_end:cos(solar_time) +
                 
                 ## Risks ##
                 # Risk from wolves
                 wolf_risk_end + 
                 wolf_risk_end:wolf_act_end + 
                 # Risk from cougars
                 cougar_risk_end + 
                 cougar_risk_end:cougar_act_end +
                 
                 ## Conditions ##
                 # Openness and roughness (above and beyond the risk effect)
                 open_end + 
                 rough_end +
                 
                 ## Movement ##
                 # Correlated SL and TA
                 log(sl_) + 
                 cos(ta_) + 
                 log(sl_):cos(ta_) +
                 # Movement by snow depth
                 log(sl_):log_swe_start + 
                 # Movement by time of day
                 log(sl_):sin(solar_time) + 
                 log(sl_):cos(solar_time) +
                 
                 # Stratum
                 strata(step_id_),
               model = TRUE)
  }))

mods <- mods %>% 
  mutate(coefs = map(issf, function(x) {
    broom::tidy(x$model)
  }))

# Save ----
system.time({
  saveRDS(mods, "../out/fitted_issa.rds")
})

# Takes about 7 minutes to save

if(!interactive()){
  message("\n========================================\nFinished at ", Sys.time())
} else {
  beepr::beep(11)
}
