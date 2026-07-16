############################################X
#--------------Elk Risk Taking--------------X
#---------------Brian J. Smith--------------X
#-----------------Model Risk----------------X
############################################X

# Set options ----
# For full printing of tibble columns
options(pillar.width = 1000)

# Load packages ----
library(dplyr)
library(amt)
library(survival)
library(mgcv)
library(tidyr)

# Custom functions
source("99_fun.R")

# Load data ----
# Fitted iSSA
dat <- readRDS("../out/fitted_issa.rds")
# All covariates
cov_dat <- dat$covs %>%
  bind_rows()

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

# Predator activity
wolf_activity <- read.csv("../elk_RT_data/wolf_activity_table.csv")
cougar_activity <- read.csv("../elk_RT_data/cougar_activity_table.csv")

# Scaling data.frame
scale_df <- read.csv("../out/scaling_data.csv")
row.names(scale_df) <- scale_df$term

# Individual metadata
meta <- read.csv("../elk_RT_data/capture_metadata.csv") %>% 
  filter(!is.na(birth_year)) %>% 
  mutate(cap_date = ymd(cap_date)) %>% 
  # Keep only first capture
  group_by(ID) %>% 
  filter(cap_date == min(cap_date))

# Calculate risk exposure ----
# Average risk available to each individual in that winter-season
risk_exp <- dat %>% 
  select(ID, winter, season, covs) %>% 
  mutate(risk_exp = map(covs, function(x) {
    x %>% 
      # Keep available steps only
      filter(!case_) %>% 
      summarize(wolf_risk_exp = mean(wolf_risk_end_orig),
                cougar_risk_exp = mean(cougar_risk_end_orig)) %>% 
      return()
  })) %>% 
  select(-covs) %>% 
  unnest(cols = risk_exp)

# Mean and SD for scaling
risk_scale <- risk_exp %>% 
  pivot_longer(wolf_risk_exp:cougar_risk_exp,
               names_to = "term", values_to = "value") %>% 
  group_by(term) %>% 
  summarize(mean = mean(value),
            sd = sd(value))

# Scale and center
risk_exp <- scale_dat(risk_exp, risk_scale)

# Save
write.csv(risk_exp, "../elk_RT_data/risk_exposure.csv", row.names = FALSE)
write.csv(risk_scale, "../elk_RT_data/risk_exposure_scaling.csv", row.names = FALSE)

# Calculate time allocation ----

# ... risk from wolves ----
# High food, high risk
x1_wolf <- expand.grid(step_id_ = 2,
                       season = c("early", "late"),
                       bio_start = scale_df["bio_start", "mean"],
                       bio_end = quantile(cov_dat$bio_end_orig, 
                                          0.9),
                       wolf_risk_start = scale_df["wolf_risk_start", "mean"],
                       wolf_risk_end = quantile(cov_dat$wolf_risk_end_orig, 
                                                0.9),
                       cougar_risk_start = scale_df["cougar_risk_start", "mean"],
                       cougar_risk_end = scale_df["cougar_risk_end", "mean"],
                       open_end = scale_df["open_end", "mean"],
                       rough_end = scale_df["rough_end", "mean"],
                       x1_ = mean(cov_dat$x1_),
                       y1_ = mean(cov_dat$y1_),
                       t1_ = seq(from = as.POSIXct("2020-03-15 00:00:00", 
                                                   tz = "US/Mountain"),
                                 to = as.POSIXct("2020-03-15 23:00:00", 
                                                 tz = "US/Mountain"),
                                 by = "1 hour")) %>% 
  # Need this to attach the hour
  mutate(t2_ = t1_ + hours(1)) %>% 
  # These are meaningless, but needed for 'predict()' or 'scale_dat()'
  mutate(sl_ = mean(cov_dat$sl_),
         ta_ = 0,
         elev_end = scale_df["elev_end", "mean"],
         elev_start = scale_df["elev_start", "mean"],
         log_swe_end = 0,
         log_swe_start = 0,
         open_start = scale_df["open_start", "mean"],
         rough_start = scale_df["rough_start", "mean"],
         swe_end = 1,
         swe_start = 1) %>% 
  attach_hour() %>% 
  scale_dat(scale_df)

# Low food, low risk
x2_wolf <- expand.grid(step_id_ = 2,
                       season = c("early", "late"),
                       bio_start = scale_df["bio_start", "mean"],
                       bio_end = quantile(cov_dat$bio_end_orig, 
                                          0.1),
                       wolf_risk_start = scale_df["wolf_risk_start", "mean"],
                       wolf_risk_end = quantile(cov_dat$wolf_risk_end_orig, 
                                                0.1),
                       cougar_risk_start = scale_df["cougar_risk_start", "mean"],
                       cougar_risk_end = scale_df["cougar_risk_end", "mean"],
                       open_end = scale_df["open_end", "mean"],
                       rough_end = scale_df["rough_end", "mean"],
                       x1_ = mean(cov_dat$x1_),
                       y1_ = mean(cov_dat$y1_),
                       t1_ = seq(from = as.POSIXct("2020-03-15 00:00:00", 
                                                   tz = "US/Mountain"),
                                 to = as.POSIXct("2020-03-15 23:00:00", 
                                                 tz = "US/Mountain"),
                                 by = "1 hour")) %>% 
  # Need this to attach the hour
  mutate(t2_ = t1_ + hours(1)) %>% 
  # These are meaningless, but needed for 'predict()' or 'scale_dat()'
  mutate(sl_ = mean(cov_dat$sl_),
         ta_ = 0,
         elev_end = scale_df["elev_end", "mean"],
         elev_start = scale_df["elev_start", "mean"],
         log_swe_end = 0,
         log_swe_start = 0,
         open_start = scale_df["open_start", "mean"],
         rough_start = scale_df["rough_start", "mean"],
         swe_end = 1,
         swe_start = 1) %>% 
  attach_hour() %>% 
  scale_dat(scale_df)

# ... ... calculate log_rss ----
# 51 sec for 10 rows
51/10*nrow(dat)/60

system.time({ # 1379 s for 1k iter
  set.seed(123456)
  lr_wolf <- boot_rss(model_df = dat,
                      niter = 1000,
                      x1 = x1_wolf,
                      x2 = x2_wolf)
})


# ... risk from cougars ----
# High food, high risk
x1_cougar <- expand.grid(step_id_ = 2,
                         season = c("early", "late"),
                         bio_start = scale_df["bio_start", "mean"],
                         bio_end = quantile(cov_dat$bio_end_orig, 
                                            0.9),
                         cougar_risk_start = scale_df["cougar_risk_start", "mean"],
                         cougar_risk_end = quantile(cov_dat$cougar_risk_end_orig, 
                                                    0.9),
                         wolf_risk_start = scale_df["wolf_risk_start", "mean"],
                         wolf_risk_end = scale_df["wolf_risk_end", "mean"],
                         open_end = scale_df["open_end", "mean"],
                         rough_end = scale_df["rough_end", "mean"],
                         x1_ = mean(cov_dat$x1_),
                         y1_ = mean(cov_dat$y1_),
                         t1_ = seq(from = as.POSIXct("2020-03-15 00:00:00", 
                                                     tz = "US/Mountain"),
                                   to = as.POSIXct("2020-03-15 23:00:00", 
                                                   tz = "US/Mountain"),
                                   by = "1 hour")) %>% 
  # Need this to attach the hour
  mutate(t2_ = t1_ + hours(1)) %>% 
  # These are meaningless, but needed for 'predict()' or 'scale_dat()'
  mutate(sl_ = mean(cov_dat$sl_),
         ta_ = 0,
         elev_end = scale_df["elev_end", "mean"],
         elev_start = scale_df["elev_start", "mean"],
         log_swe_end = 0,
         log_swe_start = 0,
         open_start = scale_df["open_start", "mean"],
         rough_start = scale_df["rough_start", "mean"],
         swe_end = 1,
         swe_start = 1) %>% 
  attach_hour() %>% 
  scale_dat(scale_df)

# Low food, low risk
x2_cougar <- expand.grid(step_id_ = 2,
                         season = c("early", "late"),
                         bio_start = scale_df["bio_start", "mean"],
                         bio_end = quantile(cov_dat$bio_end_orig, 
                                            0.1),
                         cougar_risk_start = scale_df["cougar_risk_start", "mean"],
                         cougar_risk_end = quantile(cov_dat$cougar_risk_end_orig, 
                                                    0.1),
                         wolf_risk_start = scale_df["wolf_risk_start", "mean"],
                         wolf_risk_end = scale_df["wolf_risk_end", "mean"],
                         open_end = scale_df["open_end", "mean"],
                         rough_end = scale_df["rough_end", "mean"],
                         x1_ = mean(cov_dat$x1_),
                         y1_ = mean(cov_dat$y1_),
                         t1_ = seq(from = as.POSIXct("2020-03-15 00:00:00", 
                                                     tz = "US/Mountain"),
                                   to = as.POSIXct("2020-03-15 23:00:00", 
                                                   tz = "US/Mountain"),
                                   by = "1 hour")) %>% 
  # Need this to attach the hour
  mutate(t2_ = t1_ + hours(1)) %>% 
  # These are meaningless, but needed for 'predict()' or 'scale_dat()'
  mutate(sl_ = mean(cov_dat$sl_),
         ta_ = 0,
         elev_end = scale_df["elev_end", "mean"],
         elev_start = scale_df["elev_start", "mean"],
         log_swe_end = 0,
         log_swe_start = 0,
         open_start = scale_df["open_start", "mean"],
         rough_start = scale_df["rough_start", "mean"],
         swe_end = 1,
         swe_start = 1) %>% 
  attach_hour() %>% 
  scale_dat(scale_df)

# ... ... calculate log_rss ----
system.time({ # 1479 s for 1k iter
  set.seed(654321)
  lr_cougar <- boot_rss(model_df = dat,
                        niter = 1000,
                        x1 = x1_cougar,
                        x2 = x2_cougar)
})


# Save ----
saveRDS(lr_wolf, "../out/bootstrapped_log_rss_wolf.rds")
saveRDS(lr_cougar, "../out/bootstrapped_log_rss_cougar.rds")

# If you want to load
# lr_wolf <- readRDS("../out/bootstrapped_log_rss_wolf.rds")
# lr_cougar <- readRDS("../out/bootstrapped_log_rss_cougar.rds")

# Calculate risk metrics ----

# ... wolves ----
wolf_metrics_early <- lr_wolf %>% 
  filter(season == "early") %>% 
  arrange(hour_end) %>% 
  group_by(ID, winter, season, iter) %>% 
  summarize(mean_risk = mean(log_rss),
            var_risk = var(log_rss),
            # The hour when wolves are most active
            max_risk = getElement(log_rss, 
                                  which(wolf_activity$index_early == 
                                          max(wolf_activity$index_early))),
            # The hour when wolves are least active
            min_risk = getElement(log_rss, 
                                  which(wolf_activity$index_early == 
                                          min(wolf_activity$index_early))),
            # Difference (min - max)
            #   Positive values are more risk when wolves are less active [expected]
            #   Negative values are more risk when wolves are more active [not exp]
            diff_risk = min_risk - max_risk)


wolf_metrics_late <- lr_wolf %>% 
  filter(season == "late") %>% 
  arrange(hour_end) %>% 
  group_by(ID, winter, season, iter) %>% 
  summarize(mean_risk = mean(log_rss),
            var_risk = var(log_rss),
            # The hour when wolves are most active
            max_risk = getElement(log_rss, 
                                  which(wolf_activity$index_late == 
                                          max(wolf_activity$index_late))),
            # The hour when wolves are least active
            min_risk = getElement(log_rss, 
                                  which(wolf_activity$index_late == 
                                          min(wolf_activity$index_late))),
            # Difference (min - max)
            #   Positive values are more risk when wolves are less active [expected]
            #   Negative values are more risk when wolves are more active [not exp]
            diff_risk = min_risk - max_risk)



# Combine
wolf_metrics <- rbind(wolf_metrics_early,
                      wolf_metrics_late) %>% 
  # Add ages
  left_join(select(meta, ID, birth), by = "ID") %>% 
  mutate(now = case_when(
    season == "early" ~ as.Date(paste0(winter - 1, "-12-01")),
    season == "late" ~ as.Date(paste0(winter, "-03-15"))
  )) %>% 
  mutate(age = as.numeric(difftime(now, birth, units = "days"))/365.25) %>% 
  # Add predator densities
  left_join(wolf, by = "winter") %>% 
  left_join(cougar, by = "winter") %>% 
  # Add elk density
  left_join(select(pop, winter, elk = elk_mean), by = "winter") %>% 
  # Add risk exposure
  left_join(risk_exp, 
            by = c("ID", "winter", "season"))

# ... cougars ----
cougar_metrics_early <- lr_cougar %>% 
  filter(season == "early") %>% 
  arrange(hour_end) %>% 
  group_by(ID, winter, season, iter) %>% 
  summarize(mean_risk = mean(log_rss),
            var_risk = var(log_rss),
            # The hour when cougars are most active
            max_risk = getElement(log_rss, 
                                  which(cougar_activity$index_early == 
                                          max(cougar_activity$index_early))),
            # The hour when cougars are least active
            min_risk = getElement(log_rss, 
                                  which(cougar_activity$index_early == 
                                          min(cougar_activity$index_early))),
            # Difference (min - max)
            #   Positive values are more risk when cougars are less active [expected]
            #   Negative values are more risk when cougars are more active [not exp]
            diff_risk = min_risk - max_risk)


cougar_metrics_late <- lr_cougar %>% 
  filter(season == "late") %>% 
  arrange(hour_end) %>% 
  group_by(ID, winter, season, iter) %>% 
  summarize(mean_risk = mean(log_rss),
            var_risk = var(log_rss),
            # The hour when cougars are most active
            max_risk = getElement(log_rss, 
                                  which(cougar_activity$index_late == 
                                          max(cougar_activity$index_late))),
            # The hour when cougars are least active
            min_risk = getElement(log_rss, 
                                  which(cougar_activity$index_late == 
                                          min(cougar_activity$index_late))),
            # Difference (min - max)
            #   Positive values are more risk when cougars are less active [expected]
            #   Negative values are more risk when cougars are more active [not exp]
            diff_risk = min_risk - max_risk)

# Combine
cougar_metrics <- rbind(cougar_metrics_early,
                        cougar_metrics_late) %>% 
  # Add ages
  left_join(select(meta, ID, birth), by = "ID") %>% 
  mutate(now = case_when(
    season == "early" ~ as.Date(paste0(winter - 1, "-12-01")),
    season == "late" ~ as.Date(paste0(winter, "-03-15"))
  )) %>% 
  mutate(age = as.numeric(difftime(now, birth, units = "days"))/365.25) %>% 
  # Add predator densities
  left_join(wolf, by = "winter") %>% 
  left_join(cougar, by = "winter") %>% 
  # Add elk density
  left_join(select(pop, winter, elk = elk_mean), by = "winter") %>% 
  # Add risk exposure
  left_join(risk_exp, 
            by = c("ID", "winter", "season"))


# ... ... just to check ----
## EARLY SEASON
# 'max_risk' should be same as:
mxr_early <- lr_wolf %>% 
  filter(season == "early", hour_end == 8) %>% 
  arrange(ID, winter, season, iter)
all.equal(wolf_metrics_early$max_risk, mxr_early$log_rss)

# 'min_risk' should be same as:
mnr_early <- lr_wolf %>% 
  filter(season == "early", hour_end == 5) %>% 
  arrange(ID, winter, season, iter)
all.equal(wolf_metrics_early$min_risk, mnr_early$log_rss)

# 'max_risk' should be same as:
mxr_early <- lr_cougar %>% 
  filter(season == "early", hour_end == 17) %>% 
  arrange(ID, winter, season, iter)
all.equal(cougar_metrics_early$max_risk, mxr_early$log_rss)

# 'min_risk' should be same as:
mnr_early <- lr_cougar %>% 
  filter(season == "early", hour_end == 5) %>% 
  arrange(ID, winter, season, iter)
all.equal(cougar_metrics_early$min_risk, mnr_early$log_rss)

## LATE SEASON
# 'max_risk' should be same as:
mxr_late <- lr_wolf %>% 
  filter(season == "late", hour_end == 7) %>% 
  arrange(ID, winter, season, iter)
all.equal(wolf_metrics_late$max_risk, mxr_late$log_rss)

# 'min_risk' should be same as:
mnr_late <- lr_wolf %>% 
  filter(season == "late", hour_end == 23) %>% 
  arrange(ID, winter, season, iter)
all.equal(wolf_metrics_late$min_risk, mnr_late$log_rss)

# 'max_risk' should be same as:
mxr_late <- lr_cougar %>% 
  filter(season == "late", hour_end == 20) %>% 
  arrange(ID, winter, season, iter)
all.equal(cougar_metrics_late$max_risk, mxr_late$log_rss)

# 'min_risk' should be same as:
mnr_late <- lr_cougar %>% 
  filter(season == "late", hour_end == 13) %>% 
  arrange(ID, winter, season, iter)
all.equal(cougar_metrics_late$min_risk, mnr_late$log_rss)

# Save ----
saveRDS(wolf_metrics, "../out/bootstrapped_wolf_risk_metrics.rds")
saveRDS(cougar_metrics, "../out/bootstrapped_cougar_risk_metrics.rds")

# If you want to load
# wolf_metrics <- readRDS("../out/bootstrapped_wolf_risk_metrics.rds")
# cougar_metrics <- readRDS("../out/bootstrapped_cougar_risk_metrics.rds")

# Model ----
system.time({ # 256172 s (3 days) for 2k iter; XX s for 1k iter
  mods <- lapply(1:max(wolf_metrics$iter), function (i) {
    cat("Iteration", i, "of", max(wolf_metrics$iter), "          \r")
    
    # ... wolves ----
    # Subset data
    wolf_ss <- wolf_metrics[which(wolf_metrics$iter == i), ]
    wolf_ss$ID <- factor(wolf_ss$ID)
    wolf_ss$winter <- factor(wolf_ss$winter)
    
    # Daily mean
    wolf_mean <- gam(mean_risk ~ 
                       s(age, bs = "cr", k = 20) + 
                       s(wolf, bs = "cr", k = 4) +
                       s(elk, bs = "cr", k = 4) +
                       s(wolf_risk_exp, bs = "cr", k = 20) +
                       s(winter, bs = "re") +
                       s(ID, bs = "re") +
                       season,
                     knots = list(
                       age = seq(1.5, 21.5, length.out = 20),
                       wolf = seq(2.7, 4.7, length.out = 4),
                       elk = seq(4900, 7300, length.out = 4)
                     ),
                     family = gaussian(), data = wolf_ss,
                     method = "REML")
    
    # Riskiest time
    wolf_max <- gam(max_risk ~ 
                      s(age, bs = "cr", k = 20) + 
                      s(wolf, bs = "cr", k = 4) +
                      s(elk, bs = "cr", k = 4) +
                      s(wolf_risk_exp, bs = "cr", k = 20) +
                      s(winter, bs = "re") +
                      s(ID, bs = "re") +
                      season,
                    knots = list(
                      age = seq(1.5, 21.5, length.out = 20),
                      wolf = seq(2.7, 4.7, length.out = 4),
                      elk = seq(4900, 7300, length.out = 4)
                    ),
                    family = gaussian(), data = wolf_ss,
                    method = "REML")
    
    # Safest time
    wolf_min <- gam(min_risk ~ 
                      s(age, bs = "cr", k = 20) + 
                      s(wolf, bs = "cr", k = 4) +
                      s(elk, bs = "cr", k = 4) +
                      s(wolf_risk_exp, bs = "cr", k = 20) +
                      s(winter, bs = "re") +
                      s(ID, bs = "re") +
                      season,
                    knots = list(
                      age = seq(1.5, 21.5, length.out = 20),
                      wolf = seq(2.7, 4.7, length.out = 4),
                      elk = seq(4900, 7300, length.out = 4)
                    ),
                    family = gaussian(), data = wolf_ss,
                    method = "REML")
    
    # ... cougars ----
    # Subset data
    cougar_ss <- cougar_metrics[which(cougar_metrics$iter == i), ]
    cougar_ss$ID <- factor(cougar_ss$ID)
    cougar_ss$winter <- factor(cougar_ss$winter)
    
    # Daily mean
    cougar_mean <- gam(mean_risk ~ 
                         s(age, bs = "cr", k = 20) + 
                         # Not enough unique values to support a smooth
                         cougar +
                         s(elk, bs = "cr", k = 4) +
                         s(cougar_risk_exp, bs = "cr", k = 20) +
                         s(winter, bs = "re") +
                         s(ID, bs = "re") +
                         season,
                       knots = list(
                         age = seq(1.5, 21.5, length.out = 20),
                         elk = seq(4900, 7300, length.out = 4)
                       ),
                       family = gaussian(), data = cougar_ss,
                       method = "REML")
    
    # Riskiest time
    cougar_max <- gam(max_risk ~ 
                        s(age, bs = "cr", k = 20) + 
                        # Not enough unique values to support a smooth
                        cougar +
                        s(elk, bs = "cr", k = 4) +
                        s(cougar_risk_exp, bs = "cr", k = 20) +
                        s(winter, bs = "re") +
                        s(ID, bs = "re") +
                        season,
                      knots = list(
                        age = seq(1.5, 21.5, length.out = 20),
                        elk = seq(4900, 7300, length.out = 4)
                      ),
                      family = gaussian(), data = cougar_ss,
                      method = "REML")
    
    # Safest time
    cougar_min <- gam(min_risk ~ 
                        s(age, bs = "cr", k = 20) + 
                        # Not enough unique values to support a smooth
                        cougar +
                        s(elk, bs = "cr", k = 4) +
                        s(cougar_risk_exp, bs = "cr", k = 20) +
                        s(winter, bs = "re") +
                        s(ID, bs = "re") +
                        season,
                      knots = list(
                        age = seq(1.5, 21.5, length.out = 20),
                        elk = seq(4900, 7300, length.out = 4)
                      ),
                      family = gaussian(), data = cougar_ss,
                      method = "REML")
    
    # ... list to return ----
    ll <- list(
      wolf_mean = wolf_mean,
      wolf_max = wolf_max,
      wolf_min = wolf_min,
      cougar_mean = cougar_mean,
      cougar_max = cougar_max,
      cougar_min = cougar_min
    )
    return(ll)
  })
})

# Save
saveRDS(mods, "../out/bootstrapped_risk_models.rds")

# If you want to load
# mods <- readRDS("../out/bootstrapped_risk_models.rds")

beepr::beep(11)

# Predict ----
system.time({ # 4226 s
  preds <- lapply(1:length(mods), function(i) {
    cat("Iteration", i, "of", length(mods), "        \r")
    # Get models for this iteration
    m <- mods[[i]]
    
    # ... p-values ----
    p_vals <- do.call(rbind, 
                      lapply(1:length(m), function(j) {
                        ## Get model
                        mm <- m[[j]]
                        
                        ## Smooth terms
                        # Smooth terms table
                        st <- summary(mm)$s.table
                        # Format data.frame
                        sdf <- data.frame(term = rownames(st),
                                          p = unname(st[, "p-value"]))
                        
                        ## Parametric terms
                        # Parametric terms table
                        pt <- summary(mm)$p.table
                        # Format data.frame
                        pdf <- data.frame(term = rownames(pt),
                                          p = unname(pt[, "Pr(>|t|)"]))
                        
                        ## Combine
                        res <- rbind(pdf, sdf)
                        res$model <- names(m)[j]
                        return(res)
                      })
    )
    p_vals$iter <- i
    
    # ... age prediction ----
    age_df <- data.frame(
      age = seq(4, 20, length.out = 100),
      wolf = mean(wolf_metrics$wolf),
      cougar = mean(cougar_metrics$cougar),
      elk = mean(wolf_metrics$elk),
      wolf_risk_exp = mean(wolf_metrics$wolf_risk_exp),
      cougar_risk_exp = mean(cougar_metrics$cougar_risk_exp),
      winter = "2018",
      season = "early",
      ID = wolf_metrics$ID[1]
    )
    
    age_preds <- do.call(rbind, 
                         lapply(1:length(m), function(j) {
                           mm <- m[[j]]
                           res <- age_df
                           res$log_rss <- predict(mm, 
                                                  newdata = age_df, 
                                                  type = "response",
                                                  exclude = c("s(ID)", "s(winter)"))
                           res$model <- names(m)[j]
                           return(res)
                         })
    )
    age_preds$iter <- i
    
    # ... wolf density prediction ----
    wolf_df <- data.frame(
      age = mean(wolf_metrics$age, na.rm = TRUE),
      wolf = seq(min(wolf_metrics$wolf), 
                 max(wolf_metrics$wolf), 
                 length.out = 50),
      cougar = mean(cougar_metrics$cougar),
      elk = mean(wolf_metrics$elk),
      wolf_risk_exp = mean(wolf_metrics$wolf_risk_exp),
      cougar_risk_exp = mean(cougar_metrics$cougar_risk_exp),
      winter = "2018",
      season = "early",
      ID = wolf_metrics$ID[1]
    )
    
    wolf_preds <- do.call(rbind, 
                          lapply(1:length(m), function(j) {
                            mm <- m[[j]]
                            res <- wolf_df
                            res$log_rss <- predict(mm, 
                                                   newdata = wolf_df, 
                                                   type = "response",
                                                   exclude = c("s(ID)", "s(winter)"))
                            res$model <- names(m)[j]
                            return(res)
                          })
    )
    wolf_preds$iter <- i
    
    # ... cougar density prediction ----
    cougar_df <- data.frame(
      age = mean(wolf_metrics$age, na.rm = TRUE),
      wolf = mean(wolf_metrics$wolf),
      cougar = seq(min(cougar_metrics$cougar), 
                   max(cougar_metrics$cougar), 
                   length.out = 50),
      elk = mean(wolf_metrics$elk),
      wolf_risk_exp = mean(wolf_metrics$wolf_risk_exp),
      cougar_risk_exp = mean(cougar_metrics$cougar_risk_exp),
      winter = "2018",
      season = "early",
      ID = wolf_metrics$ID[1]
    )
    
    cougar_preds <- do.call(rbind, 
                            lapply(1:length(m), function(j) {
                              mm <- m[[j]]
                              res <- cougar_df
                              res$log_rss <- predict(mm, 
                                                     newdata = cougar_df, 
                                                     type = "response",
                                                     exclude = c("s(ID)", "s(winter)"))
                              res$model <- names(m)[j]
                              return(res)
                            })
    )
    cougar_preds$iter <- i
    
    # ... elk density prediction ----
    elk_df <- data.frame(
      age = mean(wolf_metrics$age, na.rm = TRUE),
      wolf = mean(wolf_metrics$wolf),
      cougar = mean(cougar_metrics$cougar),
      elk = seq(min(wolf_metrics$elk), 
                max(wolf_metrics$elk), 
                length.out = 100),
      wolf_risk_exp = mean(wolf_metrics$wolf_risk_exp),
      cougar_risk_exp = mean(cougar_metrics$cougar_risk_exp),
      winter = "2018",
      season = "early",
      ID = wolf_metrics$ID[1]
    )
    
    elk_preds <- do.call(rbind, 
                         lapply(1:length(m), function(j) {
                           mm <- m[[j]]
                           res <- elk_df
                           res$log_rss <- predict(mm, 
                                                  newdata = elk_df, 
                                                  type = "response",
                                                  exclude = c("s(ID)", "s(winter)"))
                           res$model <- names(m)[j]
                           return(res)
                         })
    )
    elk_preds$iter <- i
    
    # ... winter prediction ----
    winter_df <- data.frame(
      age = mean(wolf_metrics$age, na.rm = TRUE),
      wolf = mean(wolf_metrics$wolf),
      cougar = mean(cougar_metrics$cougar),
      elk = mean(wolf_metrics$elk),
      winter = sort(unique(wolf_metrics$winter)),
      wolf_risk_exp = mean(wolf_metrics$wolf_risk_exp),
      cougar_risk_exp = mean(cougar_metrics$cougar_risk_exp),
      season = "early",
      ID = wolf_metrics$ID[1]
    )
    
    winter_preds <- do.call(rbind, 
                            lapply(1:length(m), function(j) {
                              mm <- m[[j]]
                              res <- winter_df
                              res$log_rss <- predict(mm, 
                                                     newdata = winter_df, 
                                                     type = "response",
                                                     exclude = c("s(ID)"))
                              res$model <- names(m)[j]
                              return(res)
                            })
    )
    winter_preds$iter <- i
    
    # ... elk ID prediction ----
    id_df <- data.frame(
      age = mean(wolf_metrics$age, na.rm = TRUE),
      wolf = mean(wolf_metrics$wolf),
      cougar = mean(cougar_metrics$cougar),
      elk = mean(wolf_metrics$elk),
      wolf_risk_exp = mean(wolf_metrics$wolf_risk_exp),
      cougar_risk_exp = mean(cougar_metrics$cougar_risk_exp),
      winter = "2018",
      season = "early",
      ID = sort(unique(wolf_metrics$ID))
    )
    
    id_preds <- do.call(rbind, 
                        lapply(1:length(m), function(j) {
                          mm <- m[[j]]
                          res <- id_df
                          res$log_rss <- predict(mm, 
                                                 newdata = id_df, 
                                                 type = "response",
                                                 exclude = c("s(winter)"))
                          res$model <- names(m)[j]
                          return(res)
                        })
    )
    id_preds$iter <- i
    
    # ... season prediction ----
    seas_df <- data.frame(
      age = mean(wolf_metrics$age, na.rm = TRUE),
      wolf = mean(wolf_metrics$wolf),
      cougar = mean(cougar_metrics$cougar),
      elk = mean(wolf_metrics$elk),
      wolf_risk_exp = mean(wolf_metrics$wolf_risk_exp),
      cougar_risk_exp = mean(cougar_metrics$cougar_risk_exp),
      winter = "2018",
      season = c("early", "late"),
      ID = wolf_metrics$ID[1]
    )
    
    seas_preds <- do.call(rbind, 
                          lapply(1:length(m), function(j) {
                            mm <- m[[j]]
                            res <- seas_df
                            res$log_rss <- predict(mm, 
                                                   newdata = seas_df, 
                                                   type = "response",
                                                   exclude = c("s(ID)", "s(winter)"))
                            res$model <- names(m)[j]
                            return(res)
                          })
    )
    seas_preds$iter <- i
    
    # ... wolf exposure prediction ----
    wolf_exp_df <- data.frame(
      age = mean(wolf_metrics$age, na.rm = TRUE),
      wolf = mean(wolf_metrics$wolf),
      cougar = mean(cougar_metrics$cougar),
      elk = mean(wolf_metrics$elk),
      wolf_risk_exp = seq(min(wolf_metrics$wolf_risk_exp), 
                          max(wolf_metrics$wolf_risk_exp), 
                          length.out = 100),
      cougar_risk_exp = mean(cougar_metrics$cougar_risk_exp),
      winter = "2018",
      season = "early",
      ID = wolf_metrics$ID[1]
    )
    
    wolf_exp_preds <- do.call(rbind, 
                              lapply(1:length(m), function(j) {
                                mm <- m[[j]]
                                res <- wolf_exp_df
                                res$log_rss <- predict(mm, 
                                                       newdata = wolf_exp_df, 
                                                       type = "response",
                                                       exclude = c("s(ID)", "s(winter)"))
                                res$model <- names(m)[j]
                                return(res)
                              })
    )
    wolf_exp_preds$iter <- i
    
    # ... cougar exposure prediction ----
    cougar_exp_df <- data.frame(
      age = mean(wolf_metrics$age, na.rm = TRUE),
      wolf = mean(wolf_metrics$wolf),
      cougar = mean(cougar_metrics$cougar),
      elk = mean(wolf_metrics$elk),
      wolf_risk_exp = mean(wolf_metrics$wolf_risk_exp),
      cougar_risk_exp = seq(min(cougar_metrics$cougar_risk_exp), 
                            max(cougar_metrics$cougar_risk_exp), 
                            length.out = 100),
      winter = "2018",
      season = "early",
      ID = wolf_metrics$ID[1]
    )
    
    cougar_exp_preds <- do.call(rbind, 
                              lapply(1:length(m), function(j) {
                                mm <- m[[j]]
                                res <- cougar_exp_df
                                res$log_rss <- predict(mm, 
                                                       newdata = cougar_exp_df, 
                                                       type = "response",
                                                       exclude = c("s(ID)", "s(winter)"))
                                res$model <- names(m)[j]
                                return(res)
                              })
    )
    cougar_exp_preds$iter <- i
    
    
    # List to return
    XX <- list(
      p = p_vals,
      age = age_preds,
      wolf = wolf_preds,
      cougar = cougar_preds,
      elk = elk_preds,
      winter = winter_preds,
      ID = id_preds,
      season = seas_preds,
      wolf_risk_exp = wolf_exp_preds,
      cougar_risk_exp = cougar_exp_preds
    )
    
    return(XX)
  })
})

# Save
saveRDS(preds, "../out/bootstrapped_model_predictions.rds")

# If you want to load
# preds <- readRDS("../out/bootstrapped_model_predictions.rds")

beepr::beep(11)
