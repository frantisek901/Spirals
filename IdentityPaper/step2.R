#### Script for processing data from Identity paper experiments
#### Now we follow with sensitivity analysis of HK model with heterogeneous parameters


## Encoding: windows-1250
## Created:  2022-11-15 FrK
## Edited:   2023-02-03 FrK

## Notes:
##
##


# Head --------------------------------------------------------------------

# Clearing all
rm(list=ls())

# Packages
library(readr)
library(dplyr)
library(tibble)
library(forcats)
library(ggplot2)
library(rstatix)



# Loading and processing data ------------------------------------------------------------

# Creating object 'raw' (tibble): Loading....
raw = read_csv("ClassicalHK_heterogenousParameters_RS01-05.csv", skip = 6) %>%
  add_row(read_csv("ClassicalHK_heterogenousParameters_RS06-10.csv", skip = 6))
for (i in seq(11, 56, 5)) {
  raw = raw %>%
    add_row(read_csv(paste0("ClassicalHK_heterogenousParameters_RS", i, "-", i + 4, ".csv"), skip = 6))
}

# Transforming 'raw' to clean 'ts'
ts = raw %>%
  # Selecting and renaming...
  select(HK_distribution = 4, Present_opinion = 5,
         2, Use_identity = 12,
         N = 3, Boundary = 7, 8, Conformity = 10, 11,
         39, diversity = 40, extremness = 41, ESBG = 42,
         boundary_mean = 43, boundary_sd = 44,
         conformity_mean = 45, conformity_sd = 46) %>%

  # Processing.
  mutate(
    across(.cols = c(1:2, 4), factor),
    even_N = ((N %%2) == 0))

# For Joining script we need these data stored under different name:
ts20 = ts



# Checking set factors and their real mean and SD ------------------------------------------------------

hist(ts$boundary_mean)
bm = group_by(ts, Boundary) %>% filter(Boundary_STD > 0) %>% get_summary_stats(boundary_mean)
bm
ts %>% filter(Boundary_STD > 0) %>%
  ggplot() +
  aes(x = boundary_mean) +
  facet_wrap(vars(Boundary), scales = "free_x", nrow = 3) +
  geom_histogram() +
  theme_light()


hist(ts$boundary_sd)
bs = group_by(ts, Boundary_STD) %>% filter(Boundary_STD > 0) %>% get_summary_stats(boundary_sd)
bs
ts %>% filter(Boundary_STD > 0) %>%
  ggplot() +
  aes(x = boundary_sd) +
  facet_wrap(vars(Boundary_STD), scales = "free_x", nrow = 2) +
  geom_histogram(binwidth = 0.001) +
  theme_light()


hist(ts$conformity_mean)
cm = group_by(ts, Conformity) %>% filter(Conformity_STD > 0) %>% get_summary_stats(conformity_mean)
cm
ts %>% filter(Conformity_STD > 0) %>%
  ggplot() +
  aes(x = conformity_mean) +
  facet_wrap(vars(Conformity), scales = "free_x", nrow = 2) +
  geom_histogram() +
  theme_light()


hist(ts$conformity_sd)
cs = group_by(ts, Conformity_STD) %>% filter(Conformity_STD > 0) %>% get_summary_stats(conformity_sd)
cs
ts %>% filter(Conformity_STD > 0) %>%
  ggplot() +
  aes(x = conformity_sd) +
  facet_wrap(vars(Conformity_STD), scales = "free_x", nrow = 2) +
  geom_histogram() +
  theme_light()

# OK, my conclusion is, that for mapping purposes the real means and SDs are well
# mapped to parameters values, but for regression models it would be safer
# to use real measured values.

# Summary statistics ------------------------------------------------------

# Basic
ts %>% get_summary_stats(ESBG, extremness, diversity)

# According all factors:
for (i in names(ts)[c(1, 5:9, 18)]) {
  print(i)
  ts %>%
    group_by(eval(str2lang(i))) %>%
    get_summary_stats(ESBG, extremness, diversity) %>%
    print()
}



# Regression --------------------------------------------------------------

# Regression just on classical conditions of HK
m1 = lm(ESBG ~ Boundary + Conformity + even_N + HK_distribution, ts)
m2 = lm(extremness ~ Boundary + Conformity + even_N + HK_distribution, ts)
m3 = lm(diversity ~ Boundary + Conformity + even_N + HK_distribution, ts)
stargazer::stargazer(m1, m2, m3, type = "text")

m21 = lm(ESBG ~ factor(Boundary) + factor(Conformity) + even_N + HK_distribution, ts)
m22 = lm(extremness ~ factor(Boundary) + factor(Conformity) + even_N + HK_distribution, ts)
m23 = lm(diversity ~ factor(Boundary) + factor(Conformity) + even_N + HK_distribution, ts)
stargazer::stargazer(m21, m22, m23, type = "text")

m11 = lm(ESBG ~ Boundary * Conformity * even_N * HK_distribution, ts)
m12 = lm(extremness ~ Boundary * Conformity * even_N * HK_distribution, ts)
m13 = lm(diversity ~ Boundary * Conformity * even_N * HK_distribution, ts)
stargazer::stargazer(m11, m12, m13, type = "text")

stargazer::stargazer(m1, m11, m21, type = "text")
stargazer::stargazer(m2, m12, m22, type = "text")
stargazer::stargazer(m3, m13, m23, type = "text")

m31 = lm(ESBG ~ Boundary + Boundary_STD + Conformity + Conformity_STD + even_N + HK_distribution, ts)
m32 = lm(extremness ~ Boundary + Boundary_STD + Conformity + Conformity_STD + even_N + HK_distribution, ts)
m33 = lm(diversity ~ Boundary + Boundary_STD + Conformity + Conformity_STD + even_N + HK_distribution, ts)
stargazer::stargazer(m31, m32, m33, type = "text")

m41 = lm(ESBG ~ Boundary * Boundary_STD + Conformity * Conformity_STD + even_N + HK_distribution, ts)
m42 = lm(extremness ~ Boundary * Boundary_STD + Conformity * Conformity_STD + even_N + HK_distribution, ts)
m43 = lm(diversity ~ Boundary * Boundary_STD + Conformity * Conformity_STD + even_N + HK_distribution, ts)
stargazer::stargazer(m41, m42, m43, type = "text")

stargazer::stargazer(m1, m31, m41, type = "text", omit.stat = c("f", "ser"))
stargazer::stargazer(m2, m32, m42, type = "text", omit.stat = c("f", "ser"))
stargazer::stargazer(m3, m33, m43, type = "text", omit.stat = c("f", "ser"))

m51 = lm(ESBG ~ factor(Boundary) * factor(Boundary_STD) + factor(Conformity) * factor(Conformity_STD) + even_N + HK_distribution, ts)
m52 = lm(extremness ~ factor(Boundary) * factor(Boundary_STD) + factor(Conformity) * factor(Conformity_STD) + even_N + HK_distribution, ts)
m53 = lm(diversity ~ factor(Boundary) * factor(Boundary_STD) + factor(Conformity) * factor(Conformity_STD) + even_N + HK_distribution, ts)
stargazer::stargazer(m51, m52, m53, type = "text")

m61 = lm(ESBG ~ factor(Conformity) * factor(Boundary) + factor(Conformity_STD) * factor(Boundary_STD) + even_N + HK_distribution, ts)
m62 = lm(extremness ~ factor(Conformity) * factor(Boundary) + factor(Conformity_STD) * factor(Boundary_STD) + even_N + HK_distribution, ts)
m63 = lm(diversity ~ factor(Conformity) * factor(Boundary) + factor(Conformity_STD) * factor(Boundary_STD) + even_N + HK_distribution, ts)
stargazer::stargazer(m61, m62, m63, type = "text")

# Just short glance that 'almost full' model address 'only' cca 8pctp more:
stargazer::stargazer(m1, m31, m41, m51, type = "text", omit.stat = c("f", "ser"))
stargazer::stargazer(m2, m32, m42, m52, type = "text", omit.stat = c("f", "ser"))
stargazer::stargazer(m3, m33, m43, m53, type = "text", omit.stat = c("f", "ser"))

m71 = lm(ESBG ~ factor(Boundary) * factor(Boundary_STD), ts)
m72 = lm(extremness ~ factor(Boundary) * factor(Boundary_STD), ts)
m73 = lm(diversity ~ factor(Boundary) * factor(Boundary_STD), ts)
stargazer::stargazer(m71, m72, m73, type = "text")

m81 = lm(ESBG ~ (Boundary) * (Boundary_STD), ts)
m82 = lm(extremness ~ (Boundary) * (Boundary_STD), ts)
m83 = lm(diversity ~ (Boundary) * (Boundary_STD), ts)
stargazer::stargazer(m81, m82, m83, type = "text")

stargazer::stargazer(m41, m71, m81, type = "text", omit.stat = c("f", "ser"))
stargazer::stargazer(m42, m72, m82, type = "text", omit.stat = c("f", "ser"))
stargazer::stargazer(m43, m73, m83, type = "text", omit.stat = c("f", "ser"))

stargazer::stargazer(m1, m31, m41, m81, type = "text", omit.stat = c("f", "ser"))
stargazer::stargazer(m2, m32, m42, m82, type = "text", omit.stat = c("f", "ser"))
stargazer::stargazer(m3, m33, m43, m83, type = "text", omit.stat = c("f", "ser"))


## Regressions using real values, not set parameters
# Regression just on classical conditions of HK
r1 = lm(ESBG ~ boundary_mean + conformity_mean + even_N + HK_distribution, ts)
r2 = lm(extremness ~ boundary_mean + conformity_mean + even_N + HK_distribution, ts)
r3 = lm(diversity ~ boundary_mean + conformity_mean + even_N + HK_distribution, ts)
stargazer::stargazer(r1, r2, r3, type = "text")

r31 = lm(ESBG ~ boundary_mean + boundary_sd + conformity_mean + conformity_sd + even_N + HK_distribution, ts)
r32 = lm(extremness ~ boundary_mean + boundary_sd + conformity_mean + conformity_sd + even_N + HK_distribution, ts)
r33 = lm(diversity ~ boundary_mean + boundary_sd + conformity_mean + conformity_sd + even_N + HK_distribution, ts)
stargazer::stargazer(r31, r32, r33, type = "text")

r41 = lm(ESBG ~ boundary_mean * boundary_sd + conformity_mean * conformity_sd + even_N + HK_distribution, ts)
r42 = lm(extremness ~ boundary_mean * boundary_sd + conformity_mean * conformity_sd + even_N + HK_distribution, ts)
r43 = lm(diversity ~ boundary_mean * boundary_sd + conformity_mean * conformity_sd + even_N + HK_distribution, ts)
stargazer::stargazer(r41, r42, r43, type = "text")

r81 = lm(ESBG ~ boundary_mean * boundary_sd, ts)
r82 = lm(extremness ~ boundary_mean * boundary_sd, ts)
r83 = lm(diversity ~ boundary_mean * boundary_sd, ts)
stargazer::stargazer(r81, r82, r83, type = "text")

stargazer::stargazer(r1, r31, r41, r81, type = "text", omit.stat = c("f", "ser"))
stargazer::stargazer(r2, r32, r42, r82, type = "text", omit.stat = c("f", "ser"))
stargazer::stargazer(r3, r33, r43, r83, type = "text", omit.stat = c("f", "ser"))

# Comparison of both measurement approaches on best models:
stargazer::stargazer(r41, r81, m41, m81, type = "text", omit.stat = c("f", "ser"))
stargazer::stargazer(r42, r82, m42, m82, type = "text", omit.stat = c("f", "ser"))
stargazer::stargazer(r43, r83, m43, m83, type = "text", omit.stat = c("f", "ser"))

# Conclusion: Surprising! Despite exact measurement of key factors,
# the regression statistics are only slightly better. And still,
# big part of variability is not addressed. Why? Is it cos' of wrong function?

m91 = nls(formula = ESBG ~ Boundary * Boundary_STD,  data = ts)
m92 = glm(extremness ~ (Boundary) * (Boundary_STD), ts)
m93 = glm(diversity ~ (Boundary) * (Boundary_STD), ts)
stargazer::stargazer(m91, m82, m83, type = "text")


# Maps ------------------------------------------------------------------

# Heat map function:
heat_map = function(.data = tm, .var = "ESBG_sd", .y = "Boundary_STD", .x = "Boundary", .title = "") {
  .data %>%
    ggplot() +
    aes(x = eval(str2lang(.x)), y = eval(str2lang(.y)), col = eval(str2lang(.var))) +
    geom_point(shape = 15, size = 15) +
    scale_color_viridis_c() +
    #scale_x_continuous(breaks = seq(0.05, 0.35, 0.05)) +
    # scale_y_continuous(breaks = seq(0.1, 1, 0.1)) +
    labs(x = .x, y = .y, col = .var, title = .title) +
    theme_light()
}


# Data preparation on Boundary:
tm = ts %>%
  group_by(Boundary, Boundary_STD) %>%
  summarise(across(.cols = diversity:ESBG,
                 list(mean = mean, sd = sd),
                 .names = "{.col}_{.fn}")) %>%
  ungroup() %>%
  mutate(Boundary_STD = factor(Boundary_STD))

.height = 7.5
.width = 29.5

# Drawing:
heat_map(.var = "ESBG_sd", .title = "Heat Map 2nd Step") %>%
  ggsave("Pics/s02map02.png", plot = ., units = "cm", height = .height, width = .width)
heat_map(.var = "ESBG_mean", .title = "Heat Map 2nd Step") %>%
  ggsave("Pics/s02map01.png", plot = ., units = "cm", height = .height, width = .width)
heat_map(.var = "extremness_sd", .title = "Heat Map 2nd Step") %>%
  ggsave("Pics/s02map04.png", plot = ., units = "cm", height = .height, width = .width)
heat_map(.var = "extremness_mean", .title = "Heat Map 2nd Step") %>%
  ggsave("Pics/s02map03.png", plot = ., units = "cm", height = .height, width = .width)
heat_map(.var = "diversity_sd", .title = "Heat Map 2nd Step") %>%
  ggsave("Pics/s02map06.png", plot = ., units = "cm", height = .height, width = .width)
heat_map(.var = "diversity_mean", .title = "Heat Map 2nd Step") %>%
  ggsave("Pics/s02map05.png", plot = ., units = "cm", height = .height, width = .width)



# Data preparation on conformity:
tm = ts %>%
  group_by(Conformity, Conformity_STD) %>%
  summarise(across(.cols = diversity:ESBG,
                   list(mean = mean, sd = sd),
                   .names = "{.col}_{.fn}")) %>%
  ungroup() %>%
  mutate(Conformity_STD = factor(Conformity_STD))

.height = 6.5
.width = 15.5

# Drawing:
heat_map(.var = "ESBG_sd", .y = "Conformity_STD", .x = "Conformity", .title = "Heat Map 2nd Step") %>%
  ggsave("Pics/s02map12.png", plot = ., units = "cm", height = .height, width = .width)
heat_map(.var = "ESBG_mean", .y = "Conformity_STD", .x = "Conformity", .title = "Heat Map 2nd Step") %>%
  ggsave("Pics/s02map11.png", plot = ., units = "cm", height = .height, width = .width)
heat_map(.var = "extremness_sd", .y = "Conformity_STD", .x = "Conformity", .title = "Heat Map 2nd Step") %>%
  ggsave("Pics/s02map14.png", plot = ., units = "cm", height = .height, width = .width)
heat_map(.var = "extremness_mean", .y = "Conformity_STD", .x = "Conformity", .title = "Heat Map 2nd Step") %>%
  ggsave("Pics/s02map13.png", plot = ., units = "cm", height = .height, width = .width)
heat_map(.var = "diversity_sd", .y = "Conformity_STD", .x = "Conformity", .title = "Heat Map 2nd Step") %>%
  ggsave("Pics/s02map16.png", plot = ., units = "cm", height = .height, width = .width)
heat_map(.var = "diversity_mean", .y = "Conformity_STD", .x = "Conformity", .title = "Heat Map 2nd Step") %>%
  ggsave("Pics/s02map15.png", plot = ., units = "cm", height = .height, width = .width)
# OK, Conformity is not that much important



heat_map_facets = function(.data = tm, .var = "ESBG_sd", .y = "Boundary_STD", .x = "Boundary",
                           .y.facet = "Conformity", .x.facet = "Conformity_STD", .title = "") {
  .data %>%
    ggplot() +
    aes(x = eval(str2lang(.x)), y = eval(str2lang(.y)), col = eval(str2lang(.var))) +
    facet_grid(rows = vars(Conformity),
               cols = vars(Conformity_STD),
               labeller = "label_both") +
    geom_point(shape = 15, size = 6) +
    scale_color_viridis_c() +
    #scale_x_continuous(breaks = seq(0.05, 0.35, 0.05)) +
    # scale_y_continuous(breaks = seq(0.1, 1, 0.1)) +
    labs(x = .x, y = .y, col = .var, title = .title) +
    theme_light()
}

# Data preparation on boundary and conformity:
tm = ts %>%
  group_by(Boundary, Boundary_STD, Conformity, Conformity_STD) %>%
  summarise(across(.cols = diversity:ESBG,
                   list(mean = mean, sd = sd),
                   .names = "{.col}_{.fn}")) %>%
  ungroup() %>%
  mutate(Conformity_STD = factor(Conformity_STD), Boundary_STD = factor(Boundary_STD))

.height = 29
.width = 47

# Drawing:
heat_map_facets(.var = "ESBG_sd", .title = "Complex Heat Map 2nd Step") %>%
  ggsave("Pics/s02map22.png", plot = ., units = "cm", height = .height, width = .width)
heat_map_facets(.var = "ESBG_mean", .title = "Complex Heat Map 2nd Step") %>%
  ggsave("Pics/s02map21.png", plot = ., units = "cm", height = .height, width = .width)
heat_map_facets(.var = "extremness_sd", .title = "Complex Heat Map 2nd Step") %>%
  ggsave("Pics/s02map24.png", plot = ., units = "cm", height = .height, width = .width)
heat_map_facets(.var = "extremness_mean", .title = "Complex Heat Map 2nd Step") %>%
  ggsave("Pics/s02map23.png", plot = ., units = "cm", height = .height, width = .width)
heat_map_facets(.var = "diversity_sd", .title = "Complex Heat Map 2nd Step") %>%
  ggsave("Pics/s02map26.png", plot = ., units = "cm", height = .height, width = .width)
heat_map_facets(.var = "diversity_mean", .title = "Complex Heat Map 2nd Step") %>%
  ggsave("Pics/s02map25.png", plot = ., units = "cm", height = .height, width = .width)
# OK, Conformity is not that much important



