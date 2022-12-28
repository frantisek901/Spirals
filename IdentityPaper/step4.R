#### Script for processing data from Identity paper experiments
#### Now we follow with sensitivity analysis of HK model with
#### individual identity and heterogeneous parameters
####


## Encoding: windows-1250
## Created:  2022-11-15 FrK
## Edited:   2022-12-28 FrK

## Notes:
##
##


# Head --------------------------------------------------------------------

# Clearing all
rm(list=ls())

# Packages
library(readr)
library(dplyr)
library(tidyr)
library(tibble)
library(forcats)
library(ggplot2)
library(rstatix)



# Loading and processing data ------------------------------------------------------------

## Here we need to get together data from Step 2 (reduced in Conformity) and Step 3

## Step 4
# Creating object 'raw' (tibble): Loading....
raw4 = read_csv("Step4_indID-hetPar_RS01-05.csv", skip = 6) %>%
  add_row(read_csv("Step4_indID-hetPar_RS06-10.csv", skip = 6))
for (i in seq(11, 56, 5)) {
  raw4 = raw4 %>%
    add_row(read_csv(paste0("Step4_indID-hetPar_RS", i, "-", i + 4, ".csv"), skip = 6))
}

# Transforming 'raw4' to clean 'ts'
ts = raw4 %>%
  # Selecting and renaming...
  select(HK_distribution = 4, 2, Use_identity = 12, 13:16,
         N = 3, Boundary = 7, 8, Conformity = 10, 11,
         39, diversity = 40, extremness = 41, ESBG = 42,
         boundary_mean = 43, boundary_sd = 44,
         conformity_mean = 45, conformity_sd = 46) %>%
  # Dropping duplicate observations due to experiment runs dubling:
  distinct() %>%

  # Processing.
  mutate(
    even_N = ((N %%2) == 0),
    across(.cols = 4:7, factor))


# Checking set factors and their real mean and SD ------------------------------------------------------

hist(ts$boundary_mean)
bm = group_by(ts, Boundary, Use_identity) %>% filter(Boundary_STD > 0) %>% get_summary_stats(boundary_mean)
bm
ts %>% filter(Boundary_STD > 0) %>%
  ggplot() +
  aes(x = boundary_mean) +
  facet_wrap(vars(Boundary), scales = "free_x", nrow = 3) +
  geom_histogram(binwidth = 0.01) +
  scale_x_continuous(breaks = seq(0.09, 0.35, 0.02)) +
  theme_light()


hist(ts$boundary_sd)
bs = group_by(ts, Boundary_STD, Use_identity) %>% filter(Boundary_STD > 0) %>% get_summary_stats(boundary_sd)
bs
ts %>% filter(Boundary_STD > 0) %>%
  ggplot() +
  aes(x = boundary_sd) +
  facet_wrap(vars(Boundary_STD), scales = "free_x", nrow = 2) +
  geom_histogram(binwidth = 0.001) +
  theme_light()


hist(ts$conformity_mean)
cm = group_by(ts, Conformity, Use_identity) %>% filter(Conformity_STD > 0) %>% get_summary_stats(conformity_mean)
cm
ts %>% filter(Conformity_STD > 0) %>%
  ggplot() +
  aes(x = conformity_mean) +
  facet_wrap(vars(Conformity), scales = "free_x", nrow = 2) +
  geom_histogram(binwidth = 0.001) +
  theme_light()


hist(ts$conformity_sd)
cs = group_by(ts, Conformity_STD, Use_identity) %>% filter(Conformity_STD > 0) %>% get_summary_stats(conformity_sd)
cs
ts %>% filter(Conformity_STD > 0) %>%
  ggplot() +
  aes(x = conformity_sd) +
  facet_wrap(vars(Conformity_STD), scales = "free_x", nrow = 2) +
  geom_histogram(binwidth = 0.001) +
  theme_light()

# OK, my conclusion is, that for mapping purposes the real means and SDs are well
# mapped to parameters values, but for regression models it would be safer
# to use real measured values.

# Summary statistics ------------------------------------------------------

# Basic
ts %>% get_summary_stats(ESBG, extremness, diversity)

# According all factors:
for (i in names(ts)[c(1, 3, 6:12, 21)]) {
  print(i)
  ts %>%
    group_by(eval(str2lang(i))) %>%
    get_summary_stats(ESBG) %>%
    print()
  ts %>%
    group_by(eval(str2lang(i))) %>%
    get_summary_stats(extremness) %>%
    print()
  ts %>%
    group_by(eval(str2lang(i))) %>%
    get_summary_stats(diversity) %>%
    print()
}



# Regression --------------------------------------------------------------

# Full model -- just additive factors
m1 = lm(ESBG ~ SPIRO_Mean + SPIRO_STD + Boundary + Boundary_STD + Conformity + Conformity_STD + even_N * HK_distribution, ts)
m2 = lm(extremness ~ SPIRO_Mean + SPIRO_STD + Boundary + Boundary_STD + Conformity + Conformity_STD + even_N * HK_distribution, ts)
m3 = lm(diversity ~ SPIRO_Mean + SPIRO_STD + Boundary + Boundary_STD + Conformity + Conformity_STD + even_N * HK_distribution, ts)
stargazer::stargazer(m1, m2, m3, type = "text")

# Full model with interaction of three main factors:
m71 = lm(ESBG ~ SPIRO_Mean * SPIRO_STD * (Boundary) * Boundary_STD + Conformity + Conformity_STD + even_N * HK_distribution, ts)
m72 = lm(extremness ~ SPIRO_Mean * SPIRO_STD * (Boundary) * Boundary_STD + Conformity + Conformity_STD + even_N * HK_distribution, ts)
m73 = lm(diversity ~ SPIRO_Mean * SPIRO_STD * (Boundary) * Boundary_STD + Conformity + Conformity_STD + even_N * HK_distribution, ts)
stargazer::stargazer(m71, m72, m73, type = "text")

# Just two main factors including their interaction:
m81 = lm(ESBG ~ SPIRO_Mean * SPIRO_STD * (Boundary) * (Boundary_STD), ts)
m82 = lm(extremness ~ SPIRO_Mean * SPIRO_STD * (Boundary) * (Boundary_STD), ts)
m83 = lm(diversity ~ SPIRO_Mean * SPIRO_STD * (Boundary) * (Boundary_STD), ts)
stargazer::stargazer(m81, m82, m83, type = "text")

stargazer::stargazer(m1, m71, m81, type = "text", omit.stat = c("f", "ser"))
stargazer::stargazer(m2, m72, m82, type = "text", omit.stat = c("f", "ser"))
stargazer::stargazer(m3, m73, m83, type = "text", omit.stat = c("f", "ser"))
stargazer::stargazer(m1, m71, m81, m2, m72, m82, m3, m73, m83, type = "text", omit.stat = c("f", "ser"))





# Maps ------------------------------------------------------------------

heat_map_facets = function(.data = tm, .var = "ESBG_mean", .x = "SPIRO_Mean", .y = "Boundary",
                           .x.facet = "SPIRO_STD", .y.facet = "Boundary_STD",
                           .title = "Complex Heat Map 4th Step") {
  .data %>%
    ggplot() +
    aes(x = eval(str2lang(.x)), y = eval(str2lang(.y)), col = eval(str2lang(.var)),
        label = round(eval(str2lang(.var)), 2)) +
    facet_grid(cols = vars(SPIRO_STD),
               rows = vars(Boundary_STD),
               labeller = "label_both") +
    geom_point(shape = 15, size = 15) +
    geom_text(col = "black", size = 4) +
    scale_color_viridis_c() +
    #scale_x_continuous(breaks = seq(0.05, 0.35, 0.05)) +
    # scale_y_continuous(breaks = seq(0.1, 1, 0.1)) +
    labs(x = .x, y = .y, col = .var, title = paste(.title, .var, sep = ": ")) +
    theme_light() +
    theme(legend.position = "top")
}

# Data preparation on boundary and conformity:
tm = ts %>%
  group_by(Boundary, Boundary_STD, SPIRO_Mean, SPIRO_STD) %>%
  summarise(across(.cols = diversity:ESBG,
                   list(mean = mean, sd = sd),
                   .names = "{.col}_{.fn}")) %>% ungroup() %>%
  mutate(Boundary = factor(Boundary), Boundary_STD = factor(Boundary_STD) %>% fct_rev())

.height = 19.5
.width = 49.5

# Drawing:
heat_map_facets(.var = "ESBG_sd") %>%
  ggsave("Pics/s04map02.png", plot = ., units = "cm", height = .height, width = .width)
heat_map_facets(.var = "ESBG_mean") %>%
  ggsave("Pics/s04map01.png", plot = ., units = "cm", height = .height, width = .width)
heat_map_facets(.var = "extremness_sd") %>%
  ggsave("Pics/s04map04.png", plot = ., units = "cm", height = .height, width = .width)
heat_map_facets(.var = "extremness_mean") %>%
  ggsave("Pics/s04map03.png", plot = ., units = "cm", height = .height, width = .width)
heat_map_facets(.var = "diversity_sd") %>%
  ggsave("Pics/s04map06.png", plot = ., units = "cm", height = .height, width = .width)
heat_map_facets(.var = "diversity_mean") %>%
  ggsave("Pics/s04map05.png", plot = ., units = "cm", height = .height, width = .width)





