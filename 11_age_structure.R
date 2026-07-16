############################################X
#--------------Elk Risk Taking--------------X
#---------------Brian J. Smith--------------X
#----------Population Age Structure---------X
############################################X

# Set options ----
# For full printing of tibble columns
options(pillar.width = 1000)

# Colors
# 1995
col1995 <- "#285CA0"
# 2009
col2009 <- "#FDC086"

my_colors <- c(col1995, col2009)

# Create directories ----
dir.create("../out/table", showWarnings = FALSE, recursive = TRUE)
dir.create("fig", showWarnings = FALSE, recursive = TRUE)

# Load packages ----
library(dplyr)
library(amt)
library(survival)
library(mgcv)
library(tidyr)
library(ggplot2)
library(patchwork)
library(ragg)

# Custom functions
source("99_fun.R")

# Load data ----
# Age reconstruction of Hoy et al. (2019)
agestr <- read.csv("../elk_RT_data/Hoy_etal_female_elk_age_reconstruction.csv")
# Bootstrapped risk models
mods <- readRDS("../out/bootstrapped_risk_models.rds")

# Predator metrics
wolf_metrics <- readRDS("../out/bootstrapped_wolf_risk_metrics.rds")
cougar_metrics <- readRDS("../out/bootstrapped_cougar_risk_metrics.rds")

# GPS collared age distribution ----
wolf_gps_age <- wolf_metrics %>% 
  filter(iter == 1) %>% 
  select(ID, winter, season, age)

cougar_gps_age <- cougar_metrics %>% 
  filter(iter == 1) %>% 
  select(ID, winter, season, age)

# Make sure wolf and cougar data.frames have the same exact information
identical(wolf_gps_age, cougar_gps_age)

age_hist <- ggplot(wolf_gps_age, aes(x = age)) +
  # facet_wrap(~ winter) +
  geom_histogram(fill = "gray90", color = "black",
                 breaks = seq(1.5, 20.5, by = 1)) +
  scale_x_continuous(name = "Age (y)", breaks = seq(2, 20, by = 3)) +
  ylab("Number of Elk-year-seasons") +
  theme_bw()

ggsave("fig/age_dist.tif", plot = age_hist,
       device = agg_tiff, width = 7, height = 4, units = "in",
       dpi = 300, compression = "lzw")

# Process age structure ----
agestr2 <- as.matrix(agestr[, 2:21])
row.names(agestr2) <- agestr[, 1]

# Note that this transposes the matrix
age_prop <- apply(agestr2, 1, function(x){
  return(x/sum(x))
  })
colSums(age_prop)

# Average age by year
age_mean <- apply(age_prop, 2, function(w){
  return(weighted.mean(1:20, w = w))
})

age_median <- apply(agestr2, 1, function(xx){
  # Essentially want a "weighted median" sort of function
  # Repeat each age for the number of individuals, then take median
  age_list <- lapply(1:20, function(a) {
    return(rep(a, times = xx[[a]]))
  })
  ages <- do.call(c, age_list)
  return(median(ages))
})

# Tidy
age_tidy <- agestr %>% 
  rename(year = age) %>% 
  pivot_longer(cols = -year,
               names_to = "age",
               names_transform = readr::parse_number,
               values_to = "n")

# Mean again (double-check)
age_tidy %>% 
  group_by(year) %>% 
  summarize(mean_age = weighted.mean(age, w = n))

# Age structure figure ----
age_struct_fig <- age_tidy %>% 
  filter(year %in% c(1995, 2009)) %>% 
  group_by(year) %>% 
  mutate(prop = n/(sum(n))) %>% 
  mutate(struct = case_when(
    year == 1995 ~ 4,
    year == 2009 ~ 10
  )) %>% 
  mutate(struct = factor(struct, levels = c(4, 10),
                         labels = c("4",
                                    "10"))) %>% 
  ggplot(aes(x = age, y = prop, fill = struct)) +
  geom_bar(stat = "identity", position = position_dodge(0.75),
           key_glyph = "rect") +
  xlab("Elk Age (y)") +
  ylab("Proportion of Population") +
  scale_fill_manual(name = "Median Age", 
                    breaks = c("4", "10"), 
                    values = my_colors) +
  theme_bw() +
  theme(legend.position = "bottom")

# Predict ----
system.time({ # 77 sec
  age_pred <- lapply(1:length(mods), function(i) {
    # For every bootstrap iteration, calculate the weighted average trait
    # for 1995 (young elk population) and 2009 (old elk population)
    
    # Status
    cat("Iteration", i, "of", length(mods), "        \r")
    
    # Get models for this iteration
    # There are 6 models in here: the 3 temporal summaries 
    # (riskiest time, safest time, daily mean) x the 2 predators
    m <- mods[[i]]
    
    # Setup age prediction
    age_df <- data.frame(
      age = 1:20,
      wolf = mean(wolf_metrics$wolf),
      cougar = mean(cougar_metrics$cougar),
      elk = mean(wolf_metrics$elk),
      wolf_risk_exp = mean(wolf_metrics$wolf_risk_exp),
      cougar_risk_exp = mean(cougar_metrics$cougar_risk_exp),
      winter = "2018",
      season = "early",
      ID = wolf_metrics$ID[1]
    )
    
    # Make predictions for each of the 6 models
    age_preds <- lapply(1:length(m), function(j) {
                           mm <- m[[j]]
                           res <- age_df
                           # Predict log-RSS
                           res$log_rss <- predict(mm, 
                                                  newdata = age_df, 
                                                  type = "response",
                                                  exclude = c("s(ID)", "s(winter)"))
                           # Attach model name
                           res$model <- names(m)[j]
                           
                           # Attach population proportion for each age class
                           res$w1995 <- age_prop[, "1995"] 
                           res$w2009 <- age_prop[, "2009"] 
                           
                           return(res)
                         }) %>% 
      # Combine the 6 models into a single data.frame
      bind_rows()
    
    # Add the bootstrap iteration number
    age_preds$iter <- i
    return(age_preds)
  }) %>% 
    # Combine all bootstrap iterations into a single data.frame
    bind_rows()
})

# Save
saveRDS(age_pred, "../out/bootstrapped_age_structure_results.rds")

# If you want to load
# age_pred <- readRDS("../out/bootstrapped_age_structure_results.rds")

# Weighted means ----
wmean <- age_pred %>% 
  group_by(model, iter) %>% 
  summarize(log_trait_1995 = weighted.mean(log_rss, w = w1995),
            log_trait_2009 = weighted.mean(log_rss, w = w2009)) %>% 
  ungroup() %>% 
  pivot_longer(cols = log_trait_1995:log_trait_2009,
               names_to = "year",
               names_prefix = "log_trait_",
               values_to = "log_trait") %>% 
  mutate(trait = exp(log_trait))

# Summarize mean and 95% CI
wmean_summ <- wmean %>% 
  group_by(model, year) %>% 
  summarize(log_mean = mean(log_trait),
            log_lwr = quantile(log_trait, 0.025),
            log_upr = quantile(log_trait, 0.975),
            mean = mean(trait),
            lwr = quantile(trait, 0.025),
            upr = quantile(trait, 0.975)) %>% 
  mutate(across(log_mean:upr, \(x) round(x, digits = 2)))

write.csv(wmean_summ, "../out/table/pop_mean_table.csv", row.names = FALSE)

# Summarize % decrease from 1995 to 2009
wmean %>% 
  select(-log_trait) %>% 
  pivot_wider(names_from = year,
              values_from = trait,
              names_prefix = "trait_") %>% 
  mutate(perc_diff = (trait_1995 - trait_2009)/trait_1995) %>% 
  group_by(model) %>% 
  summarize(mean_diff = mean(perc_diff),
            lwr = quantile(perc_diff, 0.05),
            upr = quantile(perc_diff, 0.95)) %>% 
  mutate(across(mean_diff:upr, \(x) round(x, digits = 2) * 100))

# Figure
wolf_mean_fig <- wmean %>% 
  filter(model == "wolf_mean") %>% 
  mutate(struct = case_when(
    year == 1995 ~ 4,
    year == 2009 ~ 10
  )) %>% 
  mutate(struct = factor(struct, levels = c(4, 10),
                         labels = c("4",
                                    "10"))) %>% 
  ggplot(aes(x = struct, y = trait, fill = struct)) +
  geom_boxplot() +
  xlab("Median Elk Age") +
  ylab("Population Risk Taking (log-RSS)") +
  scale_fill_manual(name = "Year", 
                    breaks = c("4",
                               "10"), 
                    values = my_colors) +
  theme_bw() +
  theme(legend.position = "none")

# Combined figure ----
comb <- age_struct_fig + wolf_mean_fig +
  plot_layout(ncol = 1, nrow = 2) +
  plot_annotation(tag_levels = "A",
                  tag_prefix = "(",
                  tag_suffix = ")")

ggsave("fig/risk/risk_age_structure.tif", plot = comb,
       device = agg_tiff, width = 3, height = 7, units = "in",
       dpi = 300, compression = "lzw")
