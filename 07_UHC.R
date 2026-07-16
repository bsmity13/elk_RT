############################################X
#--------------Elk Risk Taking--------------X
#---------------Brian J. Smith--------------X
#---------------iSSA UHC Plots--------------X
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
library(survMisc)

# Custom functions
source("99_fun.R")

# Create necessary directories ----
dir.create("fig/uhc_plots", showWarnings = FALSE)
dir.create("fig/uhc_plots/by_cov", showWarnings = FALSE)
dir.create("fig/uhc_plots/by_time", showWarnings = FALSE)

# Options ----
theme_set(theme_bw())

# Load data ----
# Data and fitted iSSA
mod <- readRDS("../out/fitted_issa.rds")

# Summary stats ----
# ... concordance ----
concord <- unlist(lapply(mod$issf, function(m) {
  survival::concordance(m$model)$concordance
}))

hist(concord)
range(concord)
mean(concord)
quantile(concord, c(0.025, 0.975))
median(concord)

# ... measure of explained randomness (MER) ----
# Concordance
mer <- unlist(lapply(mod$issf, function(m) {
  survMisc::rsq(m$model, sigD = 4)$mer
}))

hist(mer)
range(mer)
quantile(mer, c(0.025, 0.975))
mean(mer)
median(mer)

# UHC plots ----
# UHC plots typically use out-of-sample testing data, but here
# I am using in-sample data.
lapply(1:nrow(mod), function(i) {
  try({
    # Report status
    cat("Row", i, "of", nrow(mod), "...               \n")
    # Get elk, winter, season
    e <- mod$ID[i]
    w <- mod$winter[i]
    s <- mod$season[i]
    # Get the dataset
    d <- mod$covs[[i]] %>% 
      # Keep only relevant columns
      select(x1_:ta_, t1_:dt_, case_, step_id_, swe_start:log_swe_end, 
             -t1_round, -period) %>% 
      # Drop unwanted 'start' columns (protect log_swe_start)
      rename(log_swe_st = log_swe_start) %>%
      select(-ends_with("_start")) %>%
      rename(log_swe_start = log_swe_st) %>% 
      # Drop unwanted time columns
      select(-hour_end) %>% 
      # Make sure there are no NAs
      drop_na()
    # Get the model
    m <- mod$issf[[i]]
    # Prep UHC plots
    uhc <- amt:::prep_uhc.fit_clogit(object = m, 
                                     test_dat = d, 
                                     n_samp = 200,
                                     verbose = FALSE)
    # Convert to confidence envelopes
    conf <- uhc %>% 
      as.data.frame() %>% 
      conf_envelope(levels = c(0.9, 0.95, 1)) %>% 
      # Split into list elements
      split(.$var)
    
    # Plot each list element
    lapply(conf, function(x) {
      p <- x |>
        ggplot(aes(x = x)) +
        geom_ribbon(aes(ymin = CI100_lwr, ymax = CI100_upr), 
                    color = NA, fill = "gray50") +
        geom_ribbon(aes(ymin = CI95_lwr, ymax = CI95_upr), 
                    color = NA, fill = "gray70") +
        geom_ribbon(aes(ymin = CI90_lwr, ymax = CI90_upr), 
                    color = NA, fill = "gray90") +
        geom_line(aes(y = A, color = "A", linetype = "A")) +
        geom_line(aes(y = U, color = "U", linetype = "U")) +
        scale_color_manual(name = "Distribution",
                           breaks = c("U", "A"),
                           labels = c("Used", "Avail"),
                           values = c("black", "red")) +
        scale_linetype_manual(name = "Distribution",
                              breaks = c("U", "A"),
                              labels = c("Used", "Avail"),
                              values = c("solid", "dashed")
        ) +
        xlab(x$var[1]) +
        ylab("Density") +
        ggtitle(label = paste("Elk", e, w, s, "winter", 
                              paste0("(", x$var[1], ")")),
                subtitle = paste(formula(m$model))[3]) +
        theme_bw() +
        theme(plot.subtitle = element_text(size = 6))
      
      # Save
      # Saving these in two places so they're easy to flip through in either
      # By winter and season first
      dir1 <- paste0("fig/uhc_plots/by_time/", w, "/", s, "/", x$var[1])
      dir.create(dir1, recursive = TRUE)
      fn1 <- paste0(dir1, "/elk", e, "_", w, "_", s, "_", x$var[1], ".tif")
      ggsave(fn1, plot = p, device = agg_tiff, width = 12, height = 6, 
             units = "in", dpi = 150, compression = "lzw")
      # By covariate
      dir2 <- paste0("fig/uhc_plots/by_cov/", x$var[1])
      dir.create(dir2, recursive = TRUE)
      fn2 <- paste0(dir2, "/elk", e, "_", w, "_", s, "_", x$var[1], ".tif")
      ggsave(fn2, plot = p, device = agg_tiff, width = 12, height = 6, 
             units = "in", dpi = 150, compression = "lzw")
    })
  })
})
