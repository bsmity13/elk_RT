############################################X
#--------------Elk Risk Taking--------------X
#---------------Brian J. Smith--------------X
#--------------Graphical Model--------------X
############################################X

# Load packages ----
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggimage)
library(ragg)
library(patchwork)
library(magick)
library(cowplot)

# Create directory ----
# Repeating this across scripts in case they are not all run
dir.create("fig", showWarnings = FALSE)

# Silhouettes ----
cougar <- image_read_svg("../PhyloPic/cougar_Christian_Osorio.svg") %>% 
  image_ggplot()
wolf <- image_read_svg("../PhyloPic/wolf_Margot_Michaud.svg") %>% 
  image_ggplot()

# Colors ----
colors <- c("growth" = "#004880",
            "predation" = "#fc4f52",
            "fitness" = "#8310d5")

# Functions ----

# ... growth curve ----
# Hyperbola of the form y = (A/(x-c)) + d
# Rotated 45 degrees with respect to a normal parabola
# Has asymptotes at x = c and y = d

hyperbola <- function(x, A, c, d) {
  y <- (A/(x-c)) + d
  return(y)
}

# ... mortality curve ----
# Exponential decay
decay <- function(x, a, k) {
  return(a*exp(-1 * k * x))
}

# Note: fitness is growth minus mortality.


# Panel A ----
# Prime-aged prey; coursing predator
df_A <- data.frame(trait = seq(0, 1, length.out = 100),
                    growth = hyperbola(x = seq(0, 1, length.out = 100),
                                       A = 1/20,
                                       c = 1.01, d = 1.05),
                    predation = decay(
                      seq(0, 1, length.out = 100),
                      0.05, 2
                    )) %>% 
  mutate(fitness = growth - predation) %>% 
  pivot_longer(growth:fitness,
               names_to = "rate",
               values_to = "value") %>% 
  mutate(rate = factor(rate, levels = c("growth", "predation", "fitness")))

max_fit_A <- df_A %>% 
  group_by(rate) %>% 
  filter(value == max(value)) %>% 
  ungroup() %>% 
  filter(rate == "fitness")


sub_A <- ggplot(df_A, aes(x = trait, y = value, color = rate)) +
  geom_line(linewidth = 1.5, linetype = "solid", show.legend = FALSE) +
  geom_vline(xintercept = max_fit_A$trait, linetype = "dashed",
             color = "goldenrod") +
  geom_point(aes(x = trait, y = value),
             data = max_fit_A, inherit.aes = FALSE,
             shape = 21, color = "black", fill = "goldenrod",
             size = 4) +
  scale_x_continuous(name = "Defensive Trait", 
                     breaks = c(0.1, 0.9),
                     labels = c("Low", "High")) +
  ylab("Rate") +
  ggtitle("Prime-aged Elk") +
  scale_color_manual(name = "Rate", values = colors) +
  coord_cartesian(ylim = c(0, 1)) +
  theme_classic()
sub_A

# Panel B ----
# Old prey; coursing predator
df_B <- data.frame(trait = seq(0, 1, length.out = 100),
                   growth =  hyperbola(x = seq(0, 1, length.out = 100),
                                       A = 1/5,
                                       c = 1.01, d = 1.05),
                   predation = decay(
                     seq(0, 1, length.out = 100),
                     0.5, 2
                   )) %>% 
  mutate(fitness = growth - predation) %>% 
  pivot_longer(growth:fitness,
               names_to = "rate",
               values_to = "value") %>% 
  mutate(rate = factor(rate, levels = c("growth", "predation", "fitness")))

max_fit_B <- df_B %>% 
  group_by(rate) %>% 
  filter(value == max(value)) %>% 
  ungroup() %>% 
  filter(rate == "fitness")


sub_B <- ggplot(df_B, aes(x = trait, y = value, color = rate)) +
  geom_line(linewidth = 1.5, linetype = "solid", show.legend = FALSE) +
  geom_vline(xintercept = max_fit_B$trait, linetype = "dashed",
             color = "goldenrod") +
  geom_point(aes(x = trait, y = value),
             data = max_fit_B, inherit.aes = FALSE,
             shape = 21, color = "black", fill = "goldenrod",
             size = 4) +
  scale_x_continuous(name = "Defensive Trait", 
                     breaks = c(0.1, 0.9),
                     labels = c("Low", "High")) +
  ylab("Rate") +
  ggtitle("Senescent Elk") +
  scale_color_manual(name = "Rate", values = colors) +
  coord_cartesian(ylim = c(0, 1)) +
  theme_classic() +
  theme(legend.position = "bottom")
sub_B

# Panel C ----
# Prime-aged prey; ambush predator
df_C <- data.frame(trait = seq(0, 1, length.out = 100),
                    growth =  hyperbola(x = seq(0, 1, length.out = 100),
                                        A = 1/20,
                                        c = 1.01, d = 1.05),
                    predation = decay(
                      seq(0, 1, length.out = 100),
                      0.25, 4
                    )) %>% 
  mutate(fitness = growth - predation) %>% 
  pivot_longer(growth:fitness,
               names_to = "rate",
               values_to = "value") %>% 
  mutate(rate = factor(rate, levels = c("growth", "predation", "fitness")))

max_fit_C <- df_C %>% 
  group_by(rate) %>% 
  filter(value == max(value)) %>% 
  ungroup() %>% 
  filter(rate == "fitness")


sub_C <- ggplot(df_C, aes(x = trait, y = value, color = rate)) +
  geom_line(linewidth = 1.5) +
  geom_vline(xintercept = max_fit_C$trait, linetype = "dashed",
             color = "goldenrod") +
  geom_point(aes(x = trait, y = value),
             data = max_fit_C, inherit.aes = FALSE,
             shape = 21, color = "black", fill = "goldenrod",
             size = 4) +
  scale_x_continuous(name = "Defensive Trait", 
                     breaks = c(0.1, 0.9),
                     labels = c("Low", "High")) +
  ylab("Rate") +
  # ggtitle("Prime-aged Prey; Ambush Predator") +
  scale_color_manual(name = "Rate", values = colors) +
  coord_cartesian(ylim = c(0, 1)) +
  theme_classic() +
  theme(legend.position = "bottom")
sub_C

# Panel D ----
# Old prey; ambush predator
df_D <- data.frame(trait = seq(0, 1, length.out = 100),
                    growth = hyperbola(x = seq(0, 1, length.out = 100),
                                       A = 1/5,
                                       c = 1.01, d = 1.05),
                    predation = decay(
                      seq(0, 1, length.out = 100),
                      0.75, 4
                    )) %>% 
  mutate(fitness = growth - predation) %>% 
  pivot_longer(growth:fitness,
               names_to = "rate",
               values_to = "value") %>% 
  mutate(rate = factor(rate, levels = c("growth", "predation", "fitness")))

max_fit_D <- df_D %>% 
  group_by(rate) %>% 
  filter(value == max(value)) %>% 
  ungroup() %>% 
  filter(rate == "fitness")


sub_D <- ggplot(df_D, aes(x = trait, y = value, color = rate)) +
  geom_line(linewidth = 1.5) +
  geom_vline(xintercept = max_fit_D$trait, linetype = "dashed",
             color = "goldenrod") +
  geom_point(aes(x = trait, y = value),
             data = max_fit_D, inherit.aes = FALSE,
             shape = 21, color = "black", fill = "goldenrod",
             size = 4) +
  scale_x_continuous(name = "Defensive Trait", 
                     breaks = c(0.1, 0.9),
                     labels = c("Low", "High")) +
  ylab("Rate") +
  # ggtitle("Old Prey; Ambush Predator") +
  scale_color_manual(name = "Rate", values = colors) +
  coord_cartesian(ylim = c(0, 1)) +
  theme_classic() +
  theme(legend.position = "bottom")
sub_D

# ... wolf & cougar ----
wolf_img <- ggplot() +
  geom_image(
    data = tibble(trait = 0.8, value = 0.8),
    aes(x = trait, y = value, image = "../PhyloPic/wolf_Margot_Michaud.svg"),
    inherit.aes = FALSE,
    size = 1
  ) +
  theme_void() +
  # theme(plot.background = element_rect(fill = "hotpink")) +
  NULL

cougar_img <- ggplot() +
  geom_image(
    data = tibble(trait = 0.8, value = 0.8),
    aes(x = trait, y = value, image = "../PhyloPic/cougar_Christian_Osorio.svg"),
    inherit.aes = FALSE,
    size = 1
  ) +
  theme_void() +
  # theme(plot.background = element_rect(fill = "hotpink")) +
  NULL

# ... combine ----
panel_top <- sub_A + sub_C + sub_B + sub_D +
  plot_layout(ncol = 2, nrow = 2, byrow = FALSE, guides = "collect")

# Panel E ----
df_E <- expand.grid(predator = c("cougar",
                                 "wolf"),
                    age = c(7, 18)) %>% 
  mutate(trait = c(
    max_fit_C$trait,
    max_fit_A$trait,
    max_fit_D$trait,
    max_fit_B$trait
  )) %>% 
  mutate(risk_taking = 0.8 - trait)


panel_E <- ggplot(df_E, aes(x = age, y = risk_taking, group = predator)) +
  geom_line(linewidth = 1.5, color = "gray50") +
  geom_point(shape = 21, color = "black", fill = "goldenrod",
             size = 4) +
  # geom_vline(xintercept = 12, linetype = "dashed") +
  geom_image(
    data = tibble(age = 10, risk_taking = 0.30, predator = "cougar"),
    aes(image = "../PhyloPic/cougar_Christian_Osorio.svg"),
    size = 0.5
  ) +
  geom_image(
    data = tibble(age = 14, risk_taking = 0.7, predator = "wolf"),
    aes(image = "../PhyloPic/wolf_Margot_Michaud.svg"),
    size = 0.5
  ) +
  scale_x_continuous(name = "Prey Age",
                     breaks = c(7, 12, 20),
                     labels = c("Prime", "Senescent", "Maximum")) +
  coord_cartesian(ylim = c(0.1, 0.9),
                  xlim = c(5, 20)) +
  scale_y_continuous(name = "Risk Taking",
                     breaks = c(0.2, 0.8),
                     labels = c("Low", "High"),
                     sec.axis = dup_axis(name = "Defensive Trait",
                                         breaks = c(0.2, 0.8),
                                         labels = c("High", "Low"))) +
  # scale_linetype_discrete(name = "Predator") +
  theme_classic() +
  theme(legend.position = "bottom",
        margins = margin(l = 10, r = 20))


# With cowplot ----
no_leg <- theme(legend.position = "none",
                margins = margin(t = 5, r = 5, b = 5, l = 10))
blank <- ggplot() +
  theme_void() +
  theme(plot.background = element_rect(fill = "white"))

# Wolf Row
wolf_row <- plot_grid(sub_A + no_leg, blank, sub_B + no_leg,
                      nrow = 1, ncol = 3, rel_widths = c(0.45, 0.1, 0.45),
                      labels = c("(A)", "", "(B)"),
                      hjust = -0.2)
wolf_row_sil <- ggdraw(wolf_row) +
                         draw_plot(wolf_img, 0.3, 0.25, 0.4, 0.4)

# Cougar Row
cougar_row <- plot_grid(sub_C + no_leg, blank, sub_D + no_leg,
                      nrow = 1, ncol = 3, rel_widths = c(0.45, 0.1, 0.45),
                      labels = c("(C)", "", "(D)"),
                      hjust = -0.2)
cougar_row_sil <- ggdraw(cougar_row) +
                         draw_plot(cougar_img, 0.28, 0.28, 0.4, 0.4)

# Combined Row
comb_row <- plot_grid(panel_E, nrow = 1, ncol = 1, labels = "(E)",
                      hjust = -0.2)

# Legend
leg <- get_legend(sub_C)

fig <- plot_grid(wolf_row_sil, cougar_row_sil, comb_row, leg,
          ncol = 1,
          rel_heights = c(1, 1, 1, 0.3))

ggsave2("fig/graphical_model.tif", plot = fig, device = agg_tiff,
       width = 8, height = 9, units = "in", dpi = 300, compression = "lzw",
       background = "white")
