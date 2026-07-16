############################################X
#--------------Elk Risk Taking--------------X
#---------------Brian J. Smith--------------X
#---------------Visualize Risk--------------X
############################################X

# Set options ----
# For full printing of tibble columns
options(pillar.width = 1000)

# Load packages ----
library(dplyr)
library(amt)
library(survival)
library(mgcv)
library(ggplot2)
library(ragg)
library(RColorBrewer)
library(patchwork)
library(geomtextpath)

# Custom functions
source("99_fun.R")

# Create directory ----
# Repeating this across scripts in case they are not all run
dir.create("fig", showWarnings = FALSE)
dir.create("fig/risk", showWarnings = FALSE, recursive = TRUE)

# Options ----
theme_set(theme_bw())

# Load data ----
preds <- readRDS("../out/bootstrapped_model_predictions.rds")
lr_wolf <- readRDS("../out/bootstrapped_log_rss_wolf.rds")

# "Significance" ----
# 'mgcv' provides approximate p-values, but I find these rather hard to interpret
# Instead, quantify proportion of bootstrap iterations with an effect.

# ... age ----
age_diff <- lapply(preds, iter_diff, "age", mono = FALSE) %>% 
  bind_rows(.id = "iter")

diff_ci(age_diff, CI = 0.95)
diff_prop(age_diff)

# ... elk density ----
elk_diff <- lapply(preds, iter_diff, "elk", mono = FALSE) %>% 
  bind_rows(.id = "iter")

diff_ci(elk_diff, CI = 0.95)
diff_prop(elk_diff)

# ... wolf density ----
wolf_diff <- lapply(preds, iter_diff, "wolf", mono = FALSE) %>% 
  bind_rows(.id = "iter")

diff_ci(wolf_diff, CI = 0.95)
diff_prop(wolf_diff)

# ... cougar density ----
cougar_diff <- lapply(preds, iter_diff, "cougar", mono = FALSE) %>% 
  bind_rows(.id = "iter")

diff_ci(cougar_diff, CI = 0.95)
diff_prop(cougar_diff)

# ... wolf exposure ----
wolf_exp_diff <- lapply(preds, iter_diff, "wolf_risk_exp", mono = FALSE) %>% 
  bind_rows(.id = "iter")

diff_ci(wolf_exp_diff, CI = 0.95)
diff_prop(wolf_exp_diff)

# ... cougar exposure ----
cougar_exp_diff <- lapply(preds, iter_diff, "cougar_risk_exp", mono = FALSE) %>% 
  bind_rows(.id = "iter")

diff_ci(cougar_exp_diff, CI = 0.95)
diff_prop(cougar_exp_diff)

# ... season  ----
season_diff <- lapply(preds, iter_diff, "season", mono = FALSE) %>% 
  bind_rows(.id = "iter")

diff_ci(season_diff, CI = 0.95)
diff_prop(season_diff)


# Plot ----
# ... p-values ----
# * approximate p-values returned by 'mgcv'
p_dat <- lapply(preds, function(x) {
  return(x$p)
}) %>% 
  bind_rows()

p_plot <- ggplot(p_dat, aes(x = model, y = p)) +
  facet_wrap(~ term, scales = "free") +
  geom_boxplot() +
  xlab(NULL) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave("fig/risk/risk_taking_p-values.tiff", plot = p_plot, device = agg_tiff,
       width = 8, height = 10, units = "in", dpi = 300, compression = "lzw")

# ... age ----
age_dat <- lapply(preds, function(x) {
  return(x$age)
}) %>% 
  bind_rows()

age_plot <- risk_plots(age_dat, var = age, xlab = "Age (y)",
                       # These limits different than those in appendix
                       wolf_ylim = c(-2, 3),
                       cougar_ylim = c(-10, 15))

ggsave("fig/risk/risk_taking_age.tiff", plot = age_plot, device = agg_tiff,
       width = 8, height = 8, units = "in", dpi = 600, compression = "lzw")

# Age with least risk taking wrt cougars:
age_dat %>% 
  filter(model == "cougar_max") %>% 
  group_by(age) %>% 
  summarize(rt = mean(log_rss),
            se = sd(log_rss)) %>% 
  filter(rt == min(rt))

# Mean risk taking across all ages:
age_dat %>% 
  group_by(model) %>% 
  summarize(mean_rt = mean(log_rss),
            se = sd(log_rss))

# ... predator ----
wolf_dat <- lapply(preds, function(x) {
  return(x$wolf)
}) %>% 
  bind_rows()

cougar_dat <- lapply(preds, function(x) {
  return(x$cougar)
}) %>% 
  bind_rows()

pred_dat <- rbind(
  wolf_dat %>% 
    filter(grepl("wolf", model)),
  cougar_dat %>% 
    filter(grepl("cougar", model))
) %>% 
  mutate(predator = case_when(
    grepl("wolf", model) ~ wolf,
    grepl("cougar", model) ~ cougar
  ))

pred_plot <- risk_plots(pred_dat, var = predator, xlab = "Predator Density",
                        wolf_ylim = c(-3, 3),
                        cougar_ylim = c(-30, 15))

ggsave("fig/risk/risk_taking_predator.tiff", plot = pred_plot, device = agg_tiff,
       width = 6.5, height = 6.5, units = "in", dpi = 300, compression = "lzw",
       scale = 1.1)

# ... elk ----
elk_dat <- lapply(preds, function(x) {
  return(x$elk)
}) %>% 
  bind_rows()

elk_plot <- risk_plots(elk_dat, var = elk, xlab = "Elk Abundance",
                       wolf_ylim = c(-3, 3),
                       cougar_ylim = c(-30, 15))

ggsave("fig/risk/risk_taking_elk.tiff", plot = elk_plot, device = agg_tiff,
       width = 6.5, height = 6.5, units = "in", dpi = 300, compression = "lzw")

# ... predator risk exposure ----
# ... predator ----
wolf_exp_dat <- lapply(preds, function(x) {
  return(x$wolf_risk_exp)
}) %>% 
  bind_rows()

cougar_exp_dat <- lapply(preds, function(x) {
  return(x$cougar_risk_exp)
}) %>% 
  bind_rows()

risk_exp_dat <- rbind(
  wolf_exp_dat %>% 
    filter(grepl("wolf", model)),
  cougar_exp_dat %>% 
    filter(grepl("cougar", model))
) %>% 
  mutate(risk_exp = case_when(
    grepl("wolf", model) ~ wolf_risk_exp,
    grepl("cougar", model) ~ cougar_risk_exp
  ))

risk_exp_plot <- risk_plots(risk_exp_dat, var = risk_exp, xlab = "Predator Risk Exposure",
                        wolf_ylim = c(-3, 3),
                        cougar_ylim = c(-30, 15))

ggsave("fig/risk/risk_taking_risk_exp.tiff", plot = risk_exp_plot, device = agg_tiff,
       width = 6.5, height = 6.5, units = "in", dpi = 300, compression = "lzw",
       scale = 1.1)

# ... winter ----
winter_dat <- lapply(preds, function(x) {
  return(x$winter)
}) %>% 
  bind_rows()

winter_plot <- risk_plots_discrete(winter_dat, var = winter, 
                                   xlab = NULL,
                                   wolf_ylim = c(-3, 3),
                                   cougar_ylim = c(-30, 15))

ggsave("fig/risk/risk_taking_winter.tiff", plot = winter_plot, device = agg_tiff,
       width = 6.5, height = 6.5, units = "in", dpi = 300, compression = "lzw",
       scale = 1.2)

# ... season ----
season_dat <- lapply(preds, function(x) {
  return(x$season)
}) %>% 
  bind_rows()

season_plot <- risk_plots_discrete(season_dat, var = season, 
                                   xlab = NULL,
                                   wolf_ylim = c(-3, 3),
                                   cougar_ylim = c(-30, 15))

ggsave("fig/risk/risk_taking_season.tiff", plot = season_plot, device = agg_tiff,
       width = 6.5, height = 6.5, units = "in", dpi = 300, compression = "lzw")

# ... id ----
id_dat <- lapply(preds, function(x) {
  return(x$ID)
}) %>% 
  bind_rows()

id_plot <- id_dat %>% 
  ggplot(aes(x = model, y = log_rss)) +
  geom_boxplot() +
  xlab(NULL) +
  ggtitle("Individual ID Random Effect") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave("fig/risk/risk_taking_id.tiff", plot = id_plot, device = agg_tiff,
       width = 8, height = 10, units = "in", dpi = 300, compression = "lzw")

beepr::beep()
