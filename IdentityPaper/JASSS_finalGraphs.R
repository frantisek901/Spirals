#### Script for producing final graphs and tables for JASSS paper

## Encoding: windows-1250
## Created:  2023-03-01 FrK
## Updated:  2023-08-04 FrK

## NOTES:
## Script does just graphs and tables, nothing else!



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
library(ggnewscale)
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
  group_by(#Conformity, N,
           Boundary, Model) %>%
  summarise(across(diversity:ESBG, list(mean = mean, sd = sd))) %>%
  ungroup()



# Graph itself
tcs %>%
  ggplot() +
  aes(y = ESBG_mean, x = Boundary, group = Model, col = Model) +
  # facet_grid(TeX(Conformity, output = "character") ~ N, scales = "fixed",
  #            labeller = labeller(.rows = label_parsed, .cols = label_value)) +
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
tcs = tc %>% filter(Step == "Constant SPIRO, Normal Boundary") %>%
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
ggsave("SelectedPics/figure2.png", units = "mm", height = 146, width = 122)
ggsave("SelectedPics/figure2.svg", units = "mm", height = 73, width = 122)



# Fig.3 -------------------------------------------------------------------

# Data preparation
tcs = tc %>% filter(Step == "Constant SPIRO, Normal Boundary", Boundary %in% c(0.1, 0.15, 0.2, 0.25, 0.3)) %>%
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
ggsave("SelectedPics/figure3.png", units = "mm", height = 146, width = 122)
ggsave("SelectedPics/figure3.svg", units = "mm", height = 73, width = 122)



# Tables and regressions ------------------------------------------------------------------

tc1 = tc %>% filter(Step %in% c("Classical HK", "No Identity, Constant Boundary")#,
                    #Boundary <= .25, Boundary >= .15
                    ) %>%
  select(HK_distribution, SPIRO_Mean:Conformity_STD, diversity:ESBG) %>%
  mutate(across(N:Conformity_STD, ~factor(.x)))

tc4 = tc %>% filter(Step == "Normal SPIRO, Normal Boundary", Boundary %in% c(0.1, 0.15, 0.2, 0.25, 0.3)) %>%
  select(HK_distribution, SPIRO_Mean:Conformity_STD, diversity:ESBG) %>%
  mutate(across(N:Conformity_STD, ~factor(.x)))



# Regression showing interaction effect of evenness  and non-random start --------

mF = lm(ESBG ~ HK_distribution * N , data = tc1)
stargazer(mF, type =  'text', digits = 3)



# Regressions showing just main effects in Step 4 -------------------------

mE = lm(ESBG ~ SPIRO_STD + SPIRO_Mean + Boundary_STD + Boundary + HK_distribution + N + Conformity_STD + Conformity, data = tc4)
mD = lm(ESBG ~ SPIRO_Mean + Boundary_STD + Boundary + HK_distribution + N + Conformity_STD + Conformity, data = tc4)
mC = lm(ESBG ~ SPIRO_STD + Boundary_STD + Boundary + HK_distribution + N + Conformity_STD + Conformity, data = tc4)
mB = lm(ESBG ~ SPIRO_STD + SPIRO_Mean + Boundary + HK_distribution + N + Conformity_STD + Conformity, data = tc4)
mA = lm(ESBG ~ SPIRO_STD + SPIRO_Mean + Boundary_STD + HK_distribution + N + Conformity_STD + Conformity, data = tc4)
stargazer(mE, mD, mC, mB, mA, type =  'text', omit.stat = c("ser", "f"), digits = 3)



# Regressions with 4-way interactions of main variables -------------------

mG = lm(ESBG ~ SPIRO_STD * SPIRO_Mean * Boundary_STD * Boundary + HK_distribution + N + Conformity_STD + Conformity, data = tc4)
mH = lm(ESBG ~ SPIRO_Mean * Boundary_STD * Boundary + HK_distribution + N + Conformity_STD + Conformity, data = tc4)
mI = lm(ESBG ~ SPIRO_STD * Boundary_STD * Boundary + HK_distribution + N + Conformity_STD + Conformity, data = tc4)
mJ = lm(ESBG ~ SPIRO_STD * SPIRO_Mean * Boundary + HK_distribution + N + Conformity_STD + Conformity, data = tc4)
mK = lm(ESBG ~ SPIRO_STD * SPIRO_Mean * Boundary_STD + HK_distribution + N + Conformity_STD + Conformity, data = tc4)
stargazer(mG, mH, mI, mJ, mK, type = 'text', omit = 21:805, omit.stat = c("ser", "f"), digits = 3)



# Step-wise regression -------------------

mL = lm(ESBG ~ HK_distribution + N + Conformity_STD + Conformity, data = tc4)
mM = lm(ESBG ~ Boundary + HK_distribution + N + Conformity_STD + Conformity, data = tc4)
mN = lm(ESBG ~ Boundary_STD + Boundary + HK_distribution + N + Conformity_STD + Conformity, data = tc4)
mO = lm(ESBG ~ SPIRO_Mean + Boundary_STD + Boundary + HK_distribution + N + Conformity_STD + Conformity, data = tc4)
mP = lm(ESBG ~ SPIRO_STD + SPIRO_Mean + Boundary_STD + Boundary + HK_distribution + N + Conformity_STD + Conformity, data = tc4)
stargazer(mL, mM, mN, mO, mP, mG, type = 'text', omit = 21:805, omit.stat = c("ser", "f"), digits = 3)
stargazer(mL, mM, mN, mO, mP, type = 'latex', omit = 21:805, omit.stat = c("ser", "f"), digits = 3)

mQ = lm(ESBG ~ SPIRO_Mean + Boundary + HK_distribution + N + Conformity_STD + Conformity, data = tc4)
stargazer(mL, mM, mN, mQ, mO, mP, type = 'text', omit = 21:805, omit.stat = c("ser", "f"), digits = 3)


# Table of model with only main effects for proceedings -------------------

coef(summary(mE)) %>%
  knitr::kable(digits = 3, format = 'text')


# md = lm(diversity ~ SPIRO_STD + SPIRO_Mean + Boundary_STD + Boundary + HK_distribution + Conformity_STD + Conformity + N, data = tc4)
# me = lm(extremness ~ SPIRO_STD + SPIRO_Mean + Boundary_STD + Boundary + HK_distribution + Conformity_STD + Conformity + N, data = tc4)



# NEW Fig.2 -------------------------------------------------------------------

# Data preparation
tcs = tc %>% filter(Step == "Constant SPIRO, Normal Boundary") %>%
  mutate(Boundary_STD = recode(Boundary_STD, `0` = r"($\sigma_{\epsilon} = 0$)", `0.05` = r"($\sigma_{\epsilon} = 0.05$)",
                               `0.1` = r"($\sigma_{\epsilon} = 0.10$)", `0.15` = r"($\sigma_{\epsilon} = 0.15$)")) %>%
  group_by(SPIRO_Mean, Boundary_STD, Boundary) %>%
  summarise(ESBG_mean = mean(ESBG)) %>%
  ungroup() %>%
  mutate(lab_color = if_else(round(ESBG_mean, 2) < .37, "white", "black"))


# Drawing Figure 2
tcs %>%
  ggplot() +
  aes(x = SPIRO_Mean, y = Boundary, col = ESBG_mean,
      label = round(ESBG_mean, 2)) +
  facet_wrap(~ TeX(Boundary_STD, output = "character"), scales = "fixed",
             labeller = label_parsed, nrow = 1) +
  geom_point(shape = 15, size = 5.75) +
  scale_color_viridis_c() +
  new_scale_color() +
  # geom_text(aes(color = lab_color), size = 2) +
  scale_color_identity() +
  labs(y = TeX(r"( $\mu_{\epsilon}$ )"), col = "Mean ESBG", x = TeX(r"( $\mu_{SPIRO}$ )")) +
  theme_light(6) +
  theme( strip.background = element_rect(fill = "black", colour = NA),
         strip.text = element_text(colour = "white", size = 12),
         legend.position = "bottom",
         legend.text = element_text(size = 9),
         legend.title = element_text(size = 9),
         legend.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "pt"),
         legend.box.margin = margin(t = -3, r = 0, b = 0, l = 0, unit = "pt"),
         plot.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "pt"),
         axis.title = element_text(size = 12))


# Saving
ggsave("SelectedPics/JASS_newFig2.png", units = "px", height = 1475, width = 1485)
ggsave("SelectedPics/JASS_newFig2.svg", units = "px", height = 770, width = 746)



# Data preparation
tcs = tc %>% filter(Step == "Constant SPIRO, Normal Boundary") %>%
  mutate(Boundary_STD = recode(Boundary_STD, `0` = r"($\sigma_{\epsilon} = 0$)", `0.05` = r"($\sigma_{\epsilon} = 0.05$)",
                               `0.1` = r"($\sigma_{\epsilon} = 0.10$)", `0.15` = r"($\sigma_{\epsilon} = 0.15$)")) %>%
  group_by(SPIRO_Mean, Boundary_STD, Boundary, Conformity_STD) %>%
  summarise(ESBG_mean = mean(ESBG)) %>%
  ungroup() %>%
  mutate(lab_color = if_else(round(ESBG_mean, 2) < .37, "white", "black"))


# Drawing Figure 2
tcs %>%
  ggplot() +
  aes(x = SPIRO_Mean, y = Boundary, col = ESBG_mean,
      label = round(ESBG_mean, 2)) +
  facet_wrap(Conformity_STD ~ TeX(Boundary_STD, output = "character"), scales = "fixed",
             labeller = label_parsed, nrow = 2) +
  geom_point(shape = 15, size = 5.75) +
  scale_color_viridis_c() +
  new_scale_color() +
  # geom_text(aes(color = lab_color), size = 2) +
  scale_color_identity() +
  labs(y = TeX(r"( $\mu_{\epsilon}$ )"), col = "Mean ESBG", x = TeX(r"( $\mu_{SPIRO}$ )")) +
  theme_light(6) +
  theme( strip.background = element_rect(fill = "black", colour = NA),
         strip.text = element_text(colour = "white", size = 12),
         legend.position = "bottom",
         legend.text = element_text(size = 9),
         legend.title = element_text(size = 9),
         legend.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "pt"),
         legend.box.margin = margin(t = -3, r = 0, b = 0, l = 0, unit = "pt"),
         plot.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "pt"),
         axis.title = element_text(size = 12))


# Saving
ggsave("SelectedPics/JASS_NEWnewFig2.png", units = "px", height = 3732, width = 1888)


# NEW Fig.3 -------------------------------------------------------------------

# Data preparation
tcs = tc %>% filter(Step == "Constant SPIRO, Normal Boundary") %>%
  mutate(Boundary_STD = recode(Boundary_STD, `0` = r"($\sigma_{\epsilon} = 0$)", `0.05` = r"($\sigma_{\epsilon} = 0.05$)",
                               `0.1` = r"($\sigma_{\epsilon} = 0.10$)", `0.15` = r"($\sigma_{\epsilon} = 0.15$)")) %>%
  group_by(SPIRO_Mean, Boundary_STD, Boundary) %>%
  summarise(ESBG_sd = sd(ESBG)) %>%
  ungroup()


# Drawing Figure 3
tcs %>%
  ggplot() +
  aes(x = SPIRO_Mean, y = Boundary, col = ESBG_sd,
      label = round(ESBG_sd, 2)) +
  facet_wrap(~ TeX(Boundary_STD, output = "character"), scales = "fixed",
             labeller = label_parsed, nrow = 1) +
  geom_point(shape = 15, size = 5.75) +
  # geom_text(color = "white", size = 2) +
  scale_color_viridis_c() +
  labs(y = TeX(r"( $\mu_{\epsilon}$ )"), col = "Standard deviation of ESBG", x = TeX(r"( $\mu_{SPIRO}$ )")) +
  theme_light(6) +
  theme( strip.background = element_rect(fill = "black", colour = NA),
         strip.text = element_text(colour = "white", size = 12),
         legend.position = "bottom",
         legend.text = element_text(size = 9),
         legend.title = element_text(size = 9),
         legend.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "pt"),
         legend.box.margin = margin(t = -3, r = 0, b = 0, l = 0, unit = "pt"),
         plot.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "pt"),
         axis.title = element_text(size = 12))


# Saving
ggsave("SelectedPics/JASS_newFig3.png", units = "px", height = 1475, width = 1485)
ggsave("SelectedPics/JASS_newFig3.svg", units = "px", height = 740, width = 750)



# NEW Fig.4 -------------------------------------------------------------------

# Data preparation
tcs = tc %>% filter(Step == "Constant SPIRO, Normal Boundary") %>%
  mutate(Boundary_STD = recode(Boundary_STD, `0` = r"($\sigma_{\epsilon} = 0$)", `0.05` = r"($\sigma_{\epsilon} = 0.05$)",
                               `0.1` = r"($\sigma_{\epsilon} = 0.10$)", `0.15` = r"($\sigma_{\epsilon} = 0.15$)")) %>%
  group_by(SPIRO_Mean, Boundary_STD, Boundary) %>%
  summarise(diversity_mean = mean(diversity)) %>%
  ungroup()


# Drawing Figure 4
tcs %>%
  ggplot() +
  aes(x = SPIRO_Mean, y = Boundary, col = diversity_mean,
      label = round(diversity_mean, 2)) +
  facet_wrap(~ TeX(Boundary_STD, output = "character"), scales = "fixed",
             labeller = label_parsed, nrow = 1) +
  geom_point(shape = 15, size = 5.75) +
  # geom_text(color = "white", size = 2) +
  scale_color_viridis_c() +
  labs(y = TeX(r"( $\mu_{\epsilon}$ )"), col = "Mean diversity", x = TeX(r"( $\mu_{SPIRO}$ )")) +
  theme_light(6) +
  theme( strip.background = element_rect(fill = "black", colour = NA),
         strip.text = element_text(colour = "white", size = 12),
         legend.position = "bottom",
         legend.text = element_text(size = 9),
         legend.title = element_text(size = 9),
         legend.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "pt"),
         legend.box.margin = margin(t = -3, r = 0, b = 0, l = 0, unit = "pt"),
         plot.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "pt"),
         axis.title = element_text(size = 12))


# Saving
ggsave("SelectedPics/JASS_newFig4.png", units = "px", height = 1475, width = 1485)
ggsave("SelectedPics/JASS_newFig4.svg", units = "px", height = 770, width = 746)



# NEW Fig.5 -------------------------------------------------------------------

# Data preparation
tcs = tc %>% filter(Step == "Constant SPIRO, Normal Boundary") %>%
  mutate(Boundary_STD = recode(Boundary_STD, `0` = r"($\sigma_{\epsilon} = 0$)", `0.05` = r"($\sigma_{\epsilon} = 0.05$)",
                               `0.1` = r"($\sigma_{\epsilon} = 0.10$)", `0.15` = r"($\sigma_{\epsilon} = 0.15$)")) %>%
  group_by(SPIRO_Mean, Boundary_STD, Boundary) %>%
  summarise(extremness_mean = mean(extremness)) %>%
  ungroup()


# Drawing Figure 5
tcs %>%
  ggplot() +
  aes(x = SPIRO_Mean, y = Boundary, col = extremness_mean,
      label = round(extremness_mean, 2)) +
  facet_wrap(~ TeX(Boundary_STD, output = "character"), scales = "fixed",
             labeller = label_parsed, nrow = 1) +
  geom_point(shape = 15, size = 5.75) +
  # geom_text(color = "white", size = 2) +
  scale_color_viridis_c() +
  labs(y = TeX(r"( $\mu_{\epsilon}$ )"), col = "Mean extremness", x = TeX(r"( $\mu_{SPIRO}$ )")) +
  theme_light(6) +
  theme( strip.background = element_rect(fill = "black", colour = NA),
         strip.text = element_text(colour = "white", size = 12),
         legend.position = "bottom",
         legend.text = element_text(size = 9),
         legend.title = element_text(size = 9),
         legend.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "pt"),
         legend.box.margin = margin(t = -3, r = 0, b = 0, l = 0, unit = "pt"),
         plot.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "pt"),
         axis.title = element_text(size = 12))


# Saving
ggsave("SelectedPics/JASS_newFig5.png", units = "px", height = 1475, width = 1485)
ggsave("SelectedPics/JASS_newFig5.svg", units = "px", height = 770, width = 746)




# NEW Fig.6 -------------------------------------------------------------------

# Data preparation
tcs = tc %>% filter(Step == "Constant SPIRO, Normal Boundary") %>%
  mutate(Boundary_STD = recode(Boundary_STD, `0` = r"($\sigma_{\epsilon} = 0$)", `0.05` = r"($\sigma_{\epsilon} = 0.05$)",
                               `0.1` = r"($\sigma_{\epsilon} = 0.10$)", `0.15` = r"($\sigma_{\epsilon} = 0.15$)")) %>%
  group_by(SPIRO_Mean, Boundary_STD, Boundary) %>%
  summarise(ESBG_Mean = mean(ESBG), ESBG_SD = sd(ESBG)) %>%
  ungroup() %>%
  pivot_longer(cols = starts_with("ESBG"), names_prefix = "ESBG_", names_to = "Measure", values_to = "ESBG")


# Drawing Figure 6
tcs %>%
  ggplot() +
  aes(y = SPIRO_Mean, x = Boundary, col = ESBG,
      label = round(ESBG, 2)) +
  facet_grid(TeX(Boundary_STD, output = "character") + Measure ~ ., scales = "fixed",
             labeller = label_parsed, switch = "y" ) +
  geom_point(shape = 15, size = 5.75) +
  # geom_text(color = "white", size = 2) +
  scale_color_viridis_c() +
  labs(x = TeX(r"( $\mu_{\epsilon}$ )"), col = "ESBG", y = TeX(r"( $\mu_{SPIRO}$ )")) +
  theme_light(6) +
  theme( strip.background = element_rect(fill = "black", colour = NA),
         strip.text = element_text(colour = "white", size = 12),
         strip.placement = "outside",
         legend.position = "bottom",
         legend.text = element_text(size = 9),
         legend.title = element_text(size = 9),
         legend.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "pt"),
         legend.box.margin = margin(t = -3, r = 0, b = 0, l = 0, unit = "pt"),
         plot.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "pt"),
         axis.title = element_text(size = 12))


# Saving
ggsave("SelectedPics/JASS_newFig6.png", units = "px", height = 2950, width = 1400)
ggsave("SelectedPics/JASS_newFig6.svg", units = "px", height = 1050, width = 500)






# TO-DO! ------------------------------------------------------------------

# 1) Do pair-wise steps comparison step-by-step
# 2) confirm that subset of Step N+1 is equal to Step N

x = tc %>% filter(Step == "No Identity, Normal Boundary")

dx = x %>%
  group_by(Boundary, Boundary_STD, N, Conformity, Conformity_STD, HK_distribution) %>%
  summarise(n = n(), sd_ESBG = sd(ESBG), mean_ESBG = mean(ESBG)) %>% ungroup()

dy = dx %>%
  pivot_wider(id_cols = 1:5, names_from = HK_distribution, values_from = c(sd_ESBG, mean_ESBG))


dx %>%
  ggplot() +
  aes(x = sd_ESBG, fill = HK_distribution, col = HK_distribution) +
  # geom_histogram(alpha = .5) +
  geom_density(alpha = .5)

dx %>%
  ggplot() +
  aes(x = mean_ESBG, fill = HK_distribution, col = HK_distribution) +
  # geom_histogram(alpha = .5) +
  geom_density(alpha = .5)
