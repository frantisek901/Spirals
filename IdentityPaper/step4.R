#### Script for processing data from Identity paper experiments
#### Now we follow with sensitivity analysis of HK model with
#### individual identity and heterogeneous parameters
####


## Encoding: windows-1250
## Created:  2022-11-15 FrK
## Edited:   2023-01-02 FrK

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
for (i in c(seq(30, 55, 5), 39, 49, "35B", "45B")) {
  raw4 = raw4 %>%
    add_row(read_csv(paste0("Step4_indID-hetPar_RS", i, ".csv"), skip = 6))
}
for (i in seq(11, 56, 5)) {
  raw4 = raw4 %>%
    add_row(read_csv(paste0("Step4_indID-hetPar_RS", i, "-", i + 4, ".csv"), skip = 6))
}


# Creating object 'raw' (tibble): Loading....
spiro = read_csv("SPIROdistribution_indID-hetPar_RS01-10.csv", skip = 6) %>%
  add_row(read_csv("SPIROdistribution_indID-hetPar_RS20.csv", skip = 6))
for (i in seq(11, 51, 10)) {
  spiro = spiro %>%
    add_row(read_csv(paste0("SPIROdistribution_indID-hetPar_RS", i, "-", i + 9, ".csv"), skip = 6))
}
spiro = select(spiro, -1) %>% distinct()

# Transforming 'raw4' to clean 'ts'
ts = left_join(raw4, spiro, by = names(raw4)[c(2:32, 34:37)]) %>%
  # Selecting and renaming...
  select(HK_distribution = 4, 2, Use_identity = 12, 13:16,
         N = 3, Boundary = 7, 8, Conformity = 10, 11,
         39, diversity = 40, extremness = 41, ESBG = 42,
         boundary_mean = 43, boundary_sd = 44,
         conformity_mean = 45, conformity_sd = 46, 49:56) %>%
  # Dropping duplicate observations due to experiment runs dubling:
  distinct() %>% drop_na() %>%

  # Processing.
  mutate(
    even_N = ((N %%2) == 0),
    across(.cols = 4:7, factor),
    spiro_mean = ((SPIRO_0.15 * .15) + (SPIRO_0.25 * .25) + (SPIRO_0.35 * .35) + (SPIRO_0.45 * .45) + (SPIRO_0.55 * .55) + (SPIRO_0.65 * .65) + (SPIRO_0.75 * .75) + (SPIRO_0.85 * .85)) / N,
    spiro_sd = sqrt((SPIRO_0.15 * ((spiro_mean - 0.15) ^ 2) + SPIRO_0.25 * ((spiro_mean - 0.25) ^ 2) + SPIRO_0.35 * ((spiro_mean - 0.35) ^ 2) + SPIRO_0.45 * ((spiro_mean - 0.45) ^ 2) + SPIRO_0.55 * ((spiro_mean - 0.55) ^ 2) + SPIRO_0.65 * ((spiro_mean - 0.65) ^ 2) + SPIRO_0.75 * ((spiro_mean - 0.75) ^ 2) + SPIRO_0.85 * ((spiro_mean - 0.85) ^ 2)) / (N - 1)) )


## Step 4.2
# Creating object 'raw' (tibble): Loading....
raw42 = read_csv("Step4-2_indID-hetPar_RS01-05.csv", skip = 6) %>%
  add_row(read_csv("Step4-2_indID-hetPar_RS06-10.csv", skip = 6)) %>%
  add_row(read_csv("Step4-2_indID-hetPar_RS01-25.csv", skip = 6))
for (i in seq(11, 56, 5)) {
  raw42 = raw42 %>%
    add_row(read_csv(paste0("Step4-2_indID-hetPar_RS", i, "-", i + 4, ".csv"), skip = 6))
}

# Transforming 'raw4' to clean 'ts'
ts42 = raw42 %>%
  # Selecting and renaming...
  select(HK_distribution = 4, 2, Use_identity = 12, 13:16,
         N = 3, Boundary = 7, 8, Conformity = 10, 11,
         39, diversity = 40, extremness = 41, ESBG = 42,
         boundary_mean = 43, boundary_sd = 44,
         conformity_mean = 45, conformity_sd = 46, 47:54) %>%
  # Dropping duplicate observations due to experiment runs dubling:
  distinct() %>% drop_na() %>%

  # Processing.
  mutate(
    even_N = ((N %%2) == 0),
    across(.cols = 4:7, factor),
    spiro_mean = ((SPIRO_0.15 * .15) + (SPIRO_0.25 * .25) + (SPIRO_0.35 * .35) + (SPIRO_0.45 * .45) + (SPIRO_0.55 * .55) + (SPIRO_0.65 * .65) + (SPIRO_0.75 * .75) + (SPIRO_0.85 * .85)) / N,
    spiro_sd = sqrt((SPIRO_0.15 * ((spiro_mean - 0.15) ^ 2) + SPIRO_0.25 * ((spiro_mean - 0.25) ^ 2) + SPIRO_0.35 * ((spiro_mean - 0.35) ^ 2) + SPIRO_0.45 * ((spiro_mean - 0.45) ^ 2) + SPIRO_0.55 * ((spiro_mean - 0.55) ^ 2) + SPIRO_0.65 * ((spiro_mean - 0.65) ^ 2) + SPIRO_0.75 * ((spiro_mean - 0.75) ^ 2) + SPIRO_0.85 * ((spiro_mean - 0.85) ^ 2)) / (N - 1)))


# Checking RS distribution ------------------------------------------------

ts %>% count(RS) %>% #filter( n == 7680) %>% nrow()
  ggplot(aes(x = RS, y = n)) +
  geom_point(size = 3, alpha = 0.4) +
  scale_x_continuous(breaks = seq(5, 60, 5))+
  scale_y_continuous(limits = c(0, 8000), breaks = seq(0, 8000, 1000)) +
  labs(title = "Step 4: RS distribution") +
  theme_light()

ts42 %>% count(RS) %>%
  ggplot(aes(x = RS, y = n)) +
  geom_point(size = 3, alpha = 0.4) +
  scale_x_continuous(breaks = seq(5, 60, 5))+
  scale_y_continuous(limits = c(0, 5250), breaks = seq(0, 5000, 500)) +
  labs(title = "Step 4.2: RS distribution") +
  theme_light()

# What is missing?
ts42 %>% filter(RS <= 25, SPIRO_STD == 0.15) %>% count(SPIRO_Mean) %>%
  ggplot(aes(y = SPIRO_Mean, x = n, label = n)) + geom_col(alpha = 0.4) + geom_text(nudge_x = +45) + theme_light()


# Checking set factors and their real mean and SD ------------------------------------------------------

## STEP 4.2
# SPIRO SD
ggplot(ts42) +
  aes(x = SPIRO_STD, y = spiro_sd) +
  geom_boxplot() +
  geom_jitter(col = "skyblue", alpha = 0.015) +
  theme_light()

# SPIRO mean
ggplot(ts42) +
  aes(x = SPIRO_Mean, y = spiro_mean) +
  geom_boxplot() +
  geom_jitter(col = "skyblue", alpha = 0.015) +
  scale_y_continuous(limits = c(0, 0.85), breaks = seq(0, 0.85, 0.05)) +
  theme_light()


## STEP 4
# SPIRO SD
ggplot(ts) +
  aes(x = SPIRO_STD, y = spiro_sd) +
  geom_boxplot() +
  geom_jitter(col = "skyblue", alpha = 0.015) +
  theme_light()

# SPIRO mean
ggplot(ts) +
  aes(x = SPIRO_Mean, y = spiro_mean) +
  geom_boxplot() +
  geom_jitter(col = "skyblue", alpha = 0.015) +
  theme_light()

# Boundary
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

# Conformity
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
                           .x.facet = "SPIRO_STD", .y.facet = "Boundary_STD", .title = "") {
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
    labs(x = .x, y = .y, col = .var, title = paste0(.title, ": ", .var)) +
    theme_light() +
    theme(legend.position = "top")
}


## STEP 4
# Data preparation on boundary and conformity:
tm = ts %>%
  group_by(Boundary, Boundary_STD, SPIRO_Mean, SPIRO_STD) %>%
  summarise(across(.cols = diversity:ESBG,
                   list(mean = mean, sd = sd),
                   .names = "{.col}_{.fn}")) %>% ungroup() %>%
  mutate(Boundary = factor(Boundary), Boundary_STD = factor(Boundary_STD) %>% fct_rev())


# Drawing:
.height = 19.5
.width = 49.5
.tit = paste0("Complex Heat Map 4th Step (N = ", nrow(ts), " simulations)")
heat_map_facets(.var = "ESBG_sd", .title = .tit) %>%
  ggsave("Pics/s04map02.png", plot = ., units = "cm", height = .height, width = .width)
heat_map_facets(.var = "ESBG_mean", .title = .tit) %>%
  ggsave("Pics/s04map01.png", plot = ., units = "cm", height = .height, width = .width)
heat_map_facets(.var = "extremness_sd", .title = .tit) %>%
  ggsave("Pics/s04map04.png", plot = ., units = "cm", height = .height, width = .width)
heat_map_facets(.var = "extremness_mean", .title = .tit) %>%
  ggsave("Pics/s04map03.png", plot = ., units = "cm", height = .height, width = .width)
heat_map_facets(.var = "diversity_sd", .title = .tit) %>%
  ggsave("Pics/s04map06.png", plot = ., units = "cm", height = .height, width = .width)
heat_map_facets(.var = "diversity_mean", .title = .tit) %>%
  ggsave("Pics/s04map05.png", plot = ., units = "cm", height = .height, width = .width)


## STEP 4.2
# Data preparation on boundary and conformity:
tm = ts42 %>%
  group_by(Boundary, Boundary_STD, SPIRO_Mean, SPIRO_STD) %>%
  summarise(across(.cols = diversity:ESBG,
                   list(mean = mean, sd = sd),
                   .names = "{.col}_{.fn}")) %>% ungroup() %>%
  mutate(Boundary = factor(Boundary), Boundary_STD = factor(Boundary_STD) %>% fct_rev())


# Drawing:
.height = 38
.width = 110
.tit = paste0("Complex Heat Map 4th Step (N = ", nrow(ts42), " simulations)")
heat_map_facets(.var = "ESBG_sd", .title = .tit) %>%
  ggsave("Pics/s04map52.png", plot = ., units = "cm", height = .height, width = .width)
heat_map_facets(.var = "ESBG_mean", .title = .tit) %>%
  ggsave("Pics/s04map51.png", plot = ., units = "cm", height = .height, width = .width)
heat_map_facets(.var = "extremness_sd", .title = .tit) %>%
  ggsave("Pics/s04map54.png", plot = ., units = "cm", height = .height, width = .width)
heat_map_facets(.var = "extremness_mean", .title = .tit) %>%
  ggsave("Pics/s04map53.png", plot = ., units = "cm", height = .height, width = .width)
heat_map_facets(.var = "diversity_sd", .title = .tit) %>%
  ggsave("Pics/s04map56.png", plot = ., units = "cm", height = .height, width = .width)
heat_map_facets(.var = "diversity_mean", .title = .tit) %>%
  ggsave("Pics/s04map55.png", plot = ., units = "cm", height = .height, width = .width)





