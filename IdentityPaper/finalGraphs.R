#### Script for producing final graphs:
#### We want to add greek letters through LaTeX

## Encoding: windows-1250
## Created:  2023-03-01 FrK
## Updated:  2023-03-02 FrK

## NOTES:
## Script does just graphs, nothing else!



# Head --------------------------------------------------------------------

# Clearing all
rm(list=ls())

# Packages:
library(tidyverse)
library(readr)
library(dplyr)
library(tibble)
library(forcats)
library(ggplot2)
library(rstatix)
library(stringr)
library(knitr)
library(kableExtra)
library(stargazer)
library(latex2exp)

#### Reading data in:
# source("SOM_loading.R")
load("SOM.RData")



# Producing graphs --------------------------------------------------------

# Fig.1 -------------------------------------------------------------------

# Data preparation
tcs = tc %>%
  mutate(Model = recode(Step, "Classical HK" = "DHK", "No Identity, Constant Boundary" = "RHK", "No Identity, Normal Boundary" = "VB",
                       "Constant SPIRO, Normal Boundary" = "VBI", "Normal SPIRO, Normal Boundary" = "VBVI"),
         N = recode(N, `100` = "Even (N=100)", `101` = "Odd (N=101)") %>% factor(),
         Conformity = recode(Conformity, `0.2` = r"($\mu_{\alpha} = 0.2$)", `0.8` = r"($\mu_{\alpha} = 0.8$)")) %>%
  group_by(Conformity, N, Boundary, Model) %>%
  summarise(across(diversity:ESBG, list(mean = mean, sd = sd))) %>%
  ungroup()



# Graph itself
tcs %>%
  ggplot() +
  aes(y = ESBG_mean, x = Boundary, group = Model, col = Model) +
  facet_grid(TeX(Conformity, output = "character") ~ N, scales = "fixed",
             labeller = labeller(.rows = label_parsed, .cols = label_value)) +
  # facet_grid(Conformity ~ N, scales = "fixed",
  #            labeller = labeller(.rows = label_value, .cols = label_value)) +
  geom_line(linewidth = 1, alpha = 0.4) +
  geom_point(alpha = 0.7, size = 1.5) +
  scale_color_viridis_d() +
  scale_x_continuous(breaks = seq(0.1, 0.3, 0.02)) +
  scale_y_continuous(breaks = seq(0, 0.4, 0.05)) +
  labs(x = TeX(r"( $\mu_{\epsilon}$ )"), y = "Mean ESBG") +
  theme_light(6) +
  theme( strip.background = element_rect(fill = "black", colour = NA),
         strip.text = element_text(colour = "white", size = 9),
         legend.position = "bottom",
         legend.text = element_text(size = 9),
         legend.title = element_text(size = 9),
         legend.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "pt"),
         legend.box.margin = margin(t = -7, r = 0, b = -3, l = 0, unit = "pt"),
         plot.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "pt"),
         axis.title = element_text(size = 9))

# Saving pics
ggsave("SelectedPics/figure1.png", units = "mm", height = 73, width = 122)
ggsave("SelectedPics/figure1.svg", units = "mm", height = 73, width = 122)



# Fig.2 -------------------------------------------------------------------

# Data preparation
tcs = tc %>% filter(Step == "Normal SPIRO, Normal Boundary") %>%
  mutate(Boundary_STD = recode(Boundary_STD, `0` = r"($\sigma_{\epsilon} = 0$)", `0.05` = r"($\sigma_{\epsilon} = 0.05$)",
                               `0.1` = r"($\sigma_{\epsilon} = 0.10$)", `0.15` = r"($\sigma_{\epsilon} = 0.15$)")) %>%
  group_by(SPIRO_Mean, Boundary_STD, Boundary) %>%
  summarise(across(diversity:ESBG, list(mean = mean, sd = sd))) %>%
  ungroup()

# Graph itself
tcs %>%
  ggplot() +
  aes(y = ESBG_mean, x = Boundary, group = SPIRO_Mean, col = SPIRO_Mean) +
  facet_wrap(~ TeX(Boundary_STD, output = "character"), scales = "fixed",
             labeller = label_parsed) +
  geom_line(linewidth = 1, alpha = 0.4) +
  geom_point(alpha = 0.7, size = 1.5) +
  scale_color_viridis_d() +
  scale_x_continuous(breaks = seq(0.1, 0.3, 0.02)) +
  scale_y_continuous(breaks = seq(0, 0.5, 0.05)) +
  labs(x = TeX(r"( $\mu_{\epsilon}$ )"), y = "Mean ESBG", col = TeX(r"( $\mu_{SPIRO}$ )")) +
  guides(colour = guide_legend(nrow = 1)) +
  theme_light(6) +
  theme( strip.background = element_rect(fill = "black", colour = NA),
         strip.text = element_text(colour = "white", size = 9),
         legend.position = "bottom",
         legend.text = element_text(size = 9),
         legend.title = element_text(size = 9),
         legend.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "pt"),
         legend.box.margin = margin(t = -7, r = 0, b = -3, l = 0, unit = "pt"),
         plot.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "pt"),
         axis.title = element_text(size = 9))

# Saving
ggsave("SelectedPics/figure2.png", units = "mm", height = 73, width = 122)
ggsave("SelectedPics/figure2.svg", units = "mm", height = 73, width = 122)



# Fig.3 -------------------------------------------------------------------

# Data preparation
tcs = tc %>% filter(Step == "Normal SPIRO, Normal Boundary", Boundary %in% c(0.1, 0.15, 0.2, 0.25, 0.3)) %>%
  mutate(Boundary_STD = recode(Boundary_STD, `0` = r"($\sigma_{\epsilon} = 0$)", `0.05` = r"($\sigma_{\epsilon} = 0.05$)",
                               `0.1` = r"($\sigma_{\epsilon} = 0.10$)", `0.15` = r"($\sigma_{\epsilon} = 0.15$)")) %>%
  group_by(SPIRO_Mean, Boundary_STD, Boundary) %>%
  summarise(across(diversity:ESBG, list(mean = mean, sd = sd))) %>%
  ungroup() %>%
  mutate(Boundary = factor(Boundary))

# Graph itself
tcs %>%
  ggplot() +
  aes(y = ESBG_mean, group = Boundary, x = as.numeric(as.character(SPIRO_Mean)), col = Boundary) +
  facet_wrap(~ TeX(Boundary_STD, output = "character"), scales = "fixed",
             labeller = label_parsed) +
  geom_line(linewidth = 1, alpha = 0.4) +
  geom_point(alpha = 0.7, size = 1.5) +
  scale_color_viridis_d() +
  scale_x_continuous(breaks = seq(0.25, 0.85, 0.12)) +
  scale_y_continuous(breaks = seq(0, 0.5, 0.05)) +
  labs(col = TeX(r"( $\mu_{\epsilon}$ )"), y = "Mean ESBG", x = TeX(r"( $\mu_{SPIRO}$ )")) +
  guides(colour = guide_legend(nrow = 1)) +
  theme_light(6) +
  theme( strip.background = element_rect(fill = "black", colour = NA),
         strip.text = element_text(colour = "white", size = 9),
         legend.position = "bottom",
         legend.text = element_text(size = 9),
         legend.title = element_text(size = 9),
         legend.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "pt"),
         legend.box.margin = margin(t = -7, r = 0, b = -3, l = 0, unit = "pt"),
         plot.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "pt"),
         axis.title = element_text(size = 9))

# Saving
ggsave("SelectedPics/figure3.png", units = "mm", height = 73, width = 122)
ggsave("SelectedPics/figure3.svg", units = "mm", height = 73, width = 122)

