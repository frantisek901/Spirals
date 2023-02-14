#### Script for processing data from Identity paper experiments
#### Now we follow with sensitivity analysis of HK model with
#### individual identity and heterogeneous parameters
####


## Encoding: windows-1250
## Created:  2022-11-15 FrK
## Edited:   2023-02-04 FrK

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
library(stringr)



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


# Preparing additional data:
raw41 = tibble()
for (i in c(13, 16, 19, 23, 26, 29)) {
  raw41 = read_csv(paste0("Step4.1_b0", i, "_indID-hetPar_RS01-10-table.csv"), skip = 6) %>%
    add_row(raw41)
  if (i %in% c(13, 16, 19)) {
    raw41 = read_csv(paste0("Step4.1_b0", i, "_indID-hetPar_RS01-10_B-table.csv"), skip = 6) %>%
      add_row(raw41)
  }
}
for (i in c(13, 16, 19, 23, 26, 29)) {
  for (j in seq(11, 51, 10)) {
  raw41 = read_csv(paste0("Step4.1_b0", i, "_indID-hetPar_RS", j,"-", j + 9, "-table.csv"), skip = 6) %>%
    add_row(raw41)
  }
}
for (i in c(13, 16, 19, 23, 26, 29)) {
  for (j in seq(11, 51, 10)) {
    raw41 = read_csv(paste0("Step4.1_b0", i, "_indID-hetPar_RS", j,"-", j + 9, "-table.csv"), skip = 6) %>%
      add_row(raw41)
  }
  f = paste0("Step4.1_b0", i, "_indID-hetPar_RS", j,"-", j + 9, "_B-table.csv")
  if (file.exists(f)) raw41 = read_csv(f, skip = 6) %>% add_row(raw41)
}



# Transforming 'raw4' to clean 'ts'
ts = left_join(raw4, spiro, by = names(raw4)[c(2:32, 34:37)]) %>%

  # Adding extra data from raw41:
  add_row(select(raw41, -`max-ticks`, -`[step]`)) %>%

  # Selecting and renaming...
  select(HK_distribution = 4, 2, Use_identity = 12, 13:16,
         N = 3, Boundary = 7, 8, Conformity = 10, 11,
         39, diversity = 40, extremness = 41, ESBG = 42,
         boundary_mean = 43, boundary_sd = 44,
         conformity_mean = 45, conformity_sd = 46, 49:56) %>%
  # Dropping duplicate observations due to experiment runs dubling:
  drop_na() %>% distinct() %>%

  # Processing.
  mutate(
    even_N = ((N %%2) == 0),
    across(.cols = 4:7, factor),
    spiro_mean = ((SPIRO_0.15 * .15) + (SPIRO_0.25 * .25) + (SPIRO_0.35 * .35) + (SPIRO_0.45 * .45) + (SPIRO_0.55 * .55) + (SPIRO_0.65 * .65) + (SPIRO_0.75 * .75) + (SPIRO_0.85 * .85)) / N,
    spiro_sd = sqrt((SPIRO_0.15 * ((spiro_mean - 0.15) ^ 2) + SPIRO_0.25 * ((spiro_mean - 0.25) ^ 2) + SPIRO_0.35 * ((spiro_mean - 0.35) ^ 2) + SPIRO_0.45 * ((spiro_mean - 0.45) ^ 2) + SPIRO_0.55 * ((spiro_mean - 0.55) ^ 2) + SPIRO_0.65 * ((spiro_mean - 0.65) ^ 2) + SPIRO_0.75 * ((spiro_mean - 0.75) ^ 2) + SPIRO_0.85 * ((spiro_mean - 0.85) ^ 2)) / (N - 1)) )

# For Joining script we need these data stored under different name:
ts41 = ts


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




# ## Step 4.3 -------------------------------------------------------------

# Creating object 'raw' (tibble): Loading....
raw43 = read_csv("Step4-3_indID-hetPar_RS01a.csv", skip = 6) %>%
  rename(ESBG = 41) %>% mutate(across(57:58, ~as.character(.x)))
for (i in 2:9) {
  print(i)
  raw43 = raw43 %>%
    add_row(read_csv(paste0("Step4-3_indID-hetPar_RS0", i, "a.csv"), skip = 6) %>%
              rename(ESBG = 41) %>% mutate(across(57:58, ~as.character(.x))))
}
for (i in c(10:29, 32:33, 35:40)) {
  print(i)
  raw43 = raw43 %>%
    add_row(read_csv(paste0("Step4-3_indID-hetPar_RS", i, "a.csv"), skip = 6) %>%
              rename(ESBG = 41) %>% mutate(across(57:58, ~as.character(.x))))
}
for (i in 1:9) {
  print(i)
  raw43 = raw43 %>%
    add_row(read_csv(paste0("Step4-3_indID-hetPar_RS0", i, ".csv"), skip = 6) %>%
              rename(ESBG = 41) %>% mutate(across(57:58, ~as.character(.x))))
}
for (i in 10:40) {
  print(i)
  raw43 = raw43 %>%
    add_row(read_csv(paste0("Step4-3_indID-hetPar_RS", i, ".csv"), skip = 6) %>%
              rename(ESBG = 41) %>% mutate(across(57:58, ~as.character(.x))))
}


# Transforming 'raw4' to clean 'ts'
ts43 = raw43 %>%
  # Selecting and renaming...
  select(HK_distribution = 4, 2, Use_identity = 12, 13:16,
         N = 3, Boundary = 7, 8, Conformity = 10, 11,
         steps = 38, diversity = 39, extremness = 40, ESBG,
         boundary_mean = 42, boundary_sd = 43,
         conformity_mean = 44, conformity_sd = 45, 46:55, Entropy = 56,
         Fractal_dimension = 57, Fractal_R2 = 58) %>%
  rename_with(.cols = starts_with("SPIRO_count"), ~paste0(str_sub(.x, end = 6), str_sub(.x, start = -4))) %>%

  # Dropping duplicate observations due to experiment runs dubling:
  distinct() %>% drop_na() %>%

  # Processing.
  mutate(
    even_N = ((N %%2) == 0),
    across(.cols = 4:5, factor),
    Fractal_dimension = parse_number(Fractal_dimension, na = c("<RuntimePrimitiveException>", "NaN")),
    Fractal_R2 = parse_number(Fractal_R2, na = c("<RuntimePrimitiveException>", "NaN")),
    spiro_mean = ((SPIRO_0.05 * .05) + (SPIRO_0.15 * .15) + (SPIRO_0.25 * .25) + (SPIRO_0.35 * .35) + (SPIRO_0.45 * .45) + (SPIRO_0.55 * .55) + (SPIRO_0.65 * .65) + (SPIRO_0.75 * .75) + (SPIRO_0.85 * .85) + (SPIRO_0.95 * .95)) / N,
    spiro_sd = sqrt(((SPIRO_0.05 * ((spiro_mean - 0.05) ^ 2)) + (SPIRO_0.15 * ((spiro_mean - 0.15) ^ 2)) + (SPIRO_0.25 * ((spiro_mean - 0.25) ^ 2)) + (SPIRO_0.35 * ((spiro_mean - 0.35) ^ 2)) + (SPIRO_0.45 * ((spiro_mean - 0.45) ^ 2)) + (SPIRO_0.55 * ((spiro_mean - 0.55) ^ 2)) + (SPIRO_0.65 * ((spiro_mean - 0.65) ^ 2)) + (SPIRO_0.75 * ((spiro_mean - 0.75) ^ 2)) + (SPIRO_0.85 * ((spiro_mean - 0.85) ^ 2)) + (SPIRO_0.95 * ((spiro_mean - 0.95) ^ 2))) / (N - 1))
    ) %>%
  replace_na(list(Fractal_dimension = 0, Fractal_R2 = 1))




# Checking RS distribution ------------------------------------------------

ts43 %>% count(RS) %>%
  ggplot(aes(x = RS, y = n)) +
  geom_point(size = 3, alpha = 0.4) +
  geom_hline(yintercept = 98000, color = "steelblue") +
  scale_x_continuous(breaks = seq(0, 60, 2))+
  scale_y_continuous(limits = c(50000, 105000), breaks = seq(50000, 105000, 5000), minor_breaks = F) +
  labs(title = "Step 4.3: RS distribution") +
  theme_light()

ts43 %>% filter(RS %in% 2:8, as.numeric(as.character(SPIRO_STD)) > 0.05) %>%
  count(RS, SPIRO_STD)

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

heat_map_facets = function(.data = tm, .var = "ESBG_mean", .x = "Boundary_STD", .y = "Boundary",
                           .y.facet = "SPIRO_STD", .x.facet = "SPIRO_Mean", .title = "") {
  .data %>%
    ggplot() +
    aes(x = eval(str2lang(.x)), y = eval(str2lang(.y)), col = eval(str2lang(.var)),
        label = round(eval(str2lang(.var)), 2)) +
    facet_grid(cols = vars(eval(str2lang(.x.facet))),
               rows = vars(eval(str2lang(.y.facet))),
               labeller = "label_value") +  # labeller = "label_both") +
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
  mutate(Boundary = factor(Boundary), Boundary_STD = factor(Boundary_STD),
         SPIRO_STD = factor(SPIRO_STD) %>% fct_rev())


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

.height = 19.5
.width = 49.5
.tit = paste0("Complex Heat Map 4th Step reorganized (N = ", nrow(ts), " simulations)")
heat_map_facets(.var = "ESBG_sd", .title = .tit) %>%
  ggsave("Pics/s04map12.png", plot = ., units = "cm", height = .height, width = .width)
heat_map_facets(.var = "ESBG_mean", .title = .tit) %>%
  ggsave("Pics/s04map11.png", plot = ., units = "cm", height = .height, width = .width)
heat_map_facets(.var = "extremness_sd", .title = .tit) %>%
  ggsave("Pics/s04map14.png", plot = ., units = "cm", height = .height, width = .width)
heat_map_facets(.var = "extremness_mean", .title = .tit) %>%
  ggsave("Pics/s04map13.png", plot = ., units = "cm", height = .height, width = .width)
heat_map_facets(.var = "diversity_sd", .title = .tit) %>%
  ggsave("Pics/s04map16.png", plot = ., units = "cm", height = .height, width = .width)
heat_map_facets(.var = "diversity_mean", .title = .tit) %>%
  ggsave("Pics/s04map15.png", plot = ., units = "cm", height = .height, width = .width)



## STEP 4.2
# Data preparation on boundary and conformity:
tm = ts42 %>%
  group_by(Boundary, Boundary_STD, SPIRO_Mean, SPIRO_STD) %>%
  summarise(across(.cols = diversity:ESBG,
                   list(mean = mean, sd = sd),
                   .names = "{.col}_{.fn}")) %>% ungroup() %>%
  mutate(Boundary = factor(Boundary) %>% fct_rev(), Boundary_STD = factor(Boundary_STD),
         SPIRO_STD = factor(SPIRO_STD))


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


.height = 38
.width = 110
.tit = paste0("Complex Heat Map 4th Step reordered (N = ", nrow(ts42), " simulations)")
heat_map_facets(.var = "ESBG_sd", .title = .tit) %>%
  ggsave("Pics/s04map62.png", plot = ., units = "cm", height = .height, width = .width)
heat_map_facets(.var = "ESBG_mean", .title = .tit) %>%
  ggsave("Pics/s04map61.png", plot = ., units = "cm", height = .height, width = .width)
heat_map_facets(.var = "extremness_sd", .title = .tit) %>%
  ggsave("Pics/s04map64.png", plot = ., units = "cm", height = .height, width = .width)
heat_map_facets(.var = "extremness_mean", .title = .tit) %>%
  ggsave("Pics/s04map63.png", plot = ., units = "cm", height = .height, width = .width)
heat_map_facets(.var = "diversity_sd", .title = .tit) %>%
  ggsave("Pics/s04map66.png", plot = ., units = "cm", height = .height, width = .width)
heat_map_facets(.var = "diversity_mean", .title = .tit) %>%
  ggsave("Pics/s04map65.png", plot = ., units = "cm", height = .height, width = .width)


.height = 38
.width = 110
.tit = paste0("Complex Heat Map 4th Step reordered again (N = ", nrow(ts42), " simulations)")
heat_map_facets(.var = "ESBG_sd", .y = "Boundary_STD", .y.facet = "Boundary",
                .x = "SPIRO_STD", .x.facet = "SPIRO_Mean", .title = .tit) %>%
  ggsave("Pics/s04map72.png", plot = ., units = "cm", height = .height, width = .width)
heat_map_facets(.var = "ESBG_mean", .y = "Boundary_STD", .y.facet = "Boundary",
                .x = "SPIRO_STD", .x.facet = "SPIRO_Mean", .title = .tit) %>%
  ggsave("Pics/s04map71.png", plot = ., units = "cm", height = .height, width = .width)
heat_map_facets(.var = "extremness_sd", .y = "Boundary_STD", .y.facet = "Boundary",
                .x = "SPIRO_STD", .x.facet = "SPIRO_Mean", .title = .tit) %>%
  ggsave("Pics/s04map74.png", plot = ., units = "cm", height = .height, width = .width)
heat_map_facets(.var = "extremness_mean", .y = "Boundary_STD", .y.facet = "Boundary",
                .x = "SPIRO_STD", .x.facet = "SPIRO_Mean", .title = .tit) %>%
  ggsave("Pics/s04map73.png", plot = ., units = "cm", height = .height, width = .width)
heat_map_facets(.var = "diversity_sd", .y = "Boundary_STD", .y.facet = "Boundary",
                .x = "SPIRO_STD", .x.facet = "SPIRO_Mean", .title = .tit) %>%
  ggsave("Pics/s04map76.png", plot = ., units = "cm", height = .height, width = .width)
heat_map_facets(.var = "diversity_mean", .y = "Boundary_STD", .y.facet = "Boundary",
                .x = "SPIRO_STD", .x.facet = "SPIRO_Mean", .title = .tit) %>%
  ggsave("Pics/s04map75.png", plot = ., units = "cm", height = .height, width = .width)



# STEP 4.3 ----------------------------------------------------------------

heat_map_facets = function(.data = tm, .var = "ESBG_mean", .x = "Boundary_STD", .y = "Boundary",
                           .y.facet = "SPIRO_STD", .x.facet = "SPIRO_Mean", .title = "") {
  .data %>%
    ggplot() +
    aes(x = eval(str2lang(.x)), y = eval(str2lang(.y)), col = eval(str2lang(.var)),
        label = round(eval(str2lang(.var)), 2)) +
    facet_grid(cols = vars(eval(str2lang(.x.facet))),
               rows = vars(eval(str2lang(.y.facet))),
               labeller = "label_value") +  # labeller = "label_both") +
    geom_point(shape = 15, size = 3) +
    geom_text(col = "white", size = 1) +
    scale_color_viridis_c() +
    scale_x_continuous(breaks = seq(0.25, 0.75, 0.05)) +
    scale_y_continuous(breaks = seq(0.1, 0.4, 0.05)) +
    labs(x = .x, y = .y, col = .var, title = paste0(.title, ": ", .var)) +
    theme_light() +
    theme(legend.position = "top")
}

# Data preparation on boundary and conformity:
ts43 %>% count(SPIRO_STD)
tm = ts43 %>% filter(as.numeric(as.character(SPIRO_STD)) <= 0.1) %>%
  group_by(Boundary, Boundary_STD, SPIRO_Mean, SPIRO_STD) %>%
  summarise(across(.cols = diversity:ESBG,
                   list(mean = mean, sd = sd),
                   .names = "{.col}_{.fn}")) %>% ungroup() %>%
  mutate(Boundary_STD = factor(Boundary_STD) %>% fct_rev(),
         SPIRO_STD = factor(SPIRO_STD),
         SPIRO_Mean = as.character(SPIRO_Mean) %>% as.numeric())

.height = 50
.width = 57#80
.tit = paste0("Complex Heat Map 4.3rd Step (N = ", nrow(filter(ts43, as.numeric(as.character(SPIRO_STD)) <= 0.1)), " simulations)")
heat_map_facets(.var = "ESBG_sd", .y = "Boundary", .y.facet = "Boundary_STD",
                .x = "SPIRO_Mean", .x.facet = "SPIRO_STD", .title = .tit) %>%
  ggsave("Pics/s04map82.png", plot = ., units = "cm", height = .height, width = .width)
heat_map_facets(.var = "ESBG_mean", .y = "Boundary", .y.facet = "Boundary_STD",
                .x = "SPIRO_Mean", .x.facet = "SPIRO_STD", .title = .tit) %>%
  ggsave("Pics/s04map81.png", plot = ., units = "cm", height = .height, width = .width)
heat_map_facets(.var = "extremness_sd", .y = "Boundary", .y.facet = "Boundary_STD",
                .x = "SPIRO_Mean", .x.facet = "SPIRO_STD", .title = .tit) %>%
  ggsave("Pics/s04map84.png", plot = ., units = "cm", height = .height, width = .width)
heat_map_facets(.var = "extremness_mean", .y = "Boundary", .y.facet = "Boundary_STD",
                .x = "SPIRO_Mean", .x.facet = "SPIRO_STD", .title = .tit) %>%
  ggsave("Pics/s04map83.png", plot = ., units = "cm", height = .height, width = .width)
heat_map_facets(.var = "diversity_sd", .y = "Boundary", .y.facet = "Boundary_STD",
                .x = "SPIRO_Mean", .x.facet = "SPIRO_STD", .title = .tit) %>%
  ggsave("Pics/s04map86.png", plot = ., units = "cm", height = .height, width = .width)
heat_map_facets(.var = "diversity_mean", .y = "Boundary", .y.facet = "Boundary_STD",
                .x = "SPIRO_Mean", .x.facet = "SPIRO_STD", .title = .tit) %>%
  ggsave("Pics/s04map85.png", plot = ., units = "cm", height = .height, width = .width)


ts43 %>% drop_na() %>%
  ggplot() +
  aes(x = (1 + Entropy), y = (1 + Fractal_dimension), col = (SPIRO_Mean)) +
  geom_jitter(alpha = 0.01) +
  scale_x_log10() +
  scale_y_log10() +
  theme_light()

ts43 %>% drop_na() %>%
  ggplot() +
  aes(x = (1 + Entropy), y = (1 + Fractal_dimension), col = (SPIRO_Mean)) +
  geom_density_2d() +
  scale_x_log10() +
  scale_y_log10() +
  theme_light()

ts43 %>% drop_na() %>% filter(RS > 32, Fractal_dimension > 0.05) %>%
  ggplot() +
  aes(x = (1 + Entropy), y = (1 + Fractal_dimension), col = as.numeric(as.character(SPIRO_Mean))) +
  #geom_density_2d(aes(x = (1 + Entropy), y = (1 + Fractal_dimension))) +
  geom_bin2d(aes(x = (1 + Entropy), y = (1 + Fractal_dimension))) +
  geom_jitter(alpha = 0.01) +
  scale_x_log10() +
  scale_y_log10() +
  theme_light()

ts43 %>% drop_na() %>% #filter(RS > 32) %>%
  ggplot() +
  aes(x = (Entropy), y = (Fractal_dimension)) +
  geom_bin2d(bins = 100) +
  #geom_jitter(alpha = 0.01) +
  # scale_x_log10() +
  # scale_y_log10() +
  scale_color_continuous() +
  labs(title = "You were right Ashwin! vast majority of cases ends up with 'Entropy' close to 0.3 and 'Fractal dimension' close to 0, just one tile containing almost everything!") +
  theme_light()

.height = 30
.width = 40
ggsave("Pics/s04map87.png", units = "cm", height = .height, width = .width)





# Step 4.3 -- structured/random start effect ------------------------------

heat_map_facets2 = function(.data = tm, .var = "ESBG_mean", .y = "SPIRO_Mean", .x = "Boundary",
                           .x.facet1 = "Boundary_STD", .x.facet2 = "HK_distribution",
                           .y.facet = "SPIRO_STD", .title = "") {
  .data %>%
    ggplot() +
    aes(x = eval(str2lang(.x)), y = eval(str2lang(.y)), col = eval(str2lang(.var)),
        label = round(eval(str2lang(.var)), 2)) +
    facet_grid(rows = vars(eval(str2lang(.y.facet))),
               cols = vars(eval(str2lang(.x.facet1)), eval(str2lang(.x.facet2))),
               labeller = "label_value") +  # labeller = "label_both") +
    geom_point(shape = 15, size = 2) +
    # geom_text(col = "white", size = 1) +
    scale_color_viridis_c() +
    scale_y_continuous(breaks = seq(0.25, 0.75, 0.05)) +
    scale_x_continuous(breaks = seq(0.1, 0.4, 0.1)) +
    labs(x = .x, y = .y, col = .var, title = paste0(.title, ": ", .var)) +
    theme_light() +
    theme(legend.position = "top")
}

tm = ts43 %>%
  group_by(HK_distribution, Boundary, Boundary_STD, SPIRO_Mean, SPIRO_STD) %>%
  summarise(across(.cols = diversity:ESBG,
                   list(mean = mean, sd = sd),
                   .names = "{.col}_{.fn}")) %>% ungroup() %>%
  mutate(Boundary_STD = factor(Boundary_STD),
         SPIRO_STD = factor(SPIRO_STD) %>% fct_rev(),
         SPIRO_Mean = as.character(SPIRO_Mean) %>% as.numeric()) %>%
  filter(ESBG_mean <= 0.55, ESBG_sd <= 0.35,
         diversity_mean <= 0.59, diversity_sd <= 0.35,
         as.numeric(as.character(SPIRO_STD)) <= 0.15)

.height = 57
.width = 68
.tit = paste0("Complex Heat Map 4.3rd Step (N = ", nrow(ts43), " simulations); Panels are organized by 'HK_distribution' (TRUE = structured/ FALSE = random start), 'Boundary_STD' and 'SPIRO_STD'")
heat_map_facets2(.var = "ESBG_mean", .x = "Boundary",
                .x.facet2 = "Boundary_STD", .x.facet1 = "HK_distribution",
                .y = "SPIRO_Mean", .y.facet = "SPIRO_STD", .title = .tit) %>%
  ggsave("Pics/s04map91.png", plot = ., units = "cm", height = .height, width = .width)
heat_map_facets2(.var = "ESBG_sd", .x = "Boundary",
                 .x.facet2 = "Boundary_STD", .x.facet1 = "HK_distribution",
                 .y = "SPIRO_Mean", .y.facet = "SPIRO_STD", .title = .tit) %>%
  ggsave("Pics/s04map92.png", plot = ., units = "cm", height = .height, width = .width)
heat_map_facets2(.var = "extremness_mean", .x = "Boundary",
                 .x.facet2 = "Boundary_STD", .x.facet1 = "HK_distribution",
                 .y = "SPIRO_Mean", .y.facet = "SPIRO_STD", .title = .tit) %>%
  ggsave("Pics/s04map93.png", plot = ., units = "cm", height = .height, width = .width)
heat_map_facets2(.var = "extremness_sd", .x = "Boundary",
                 .x.facet2 = "Boundary_STD", .x.facet1 = "HK_distribution",
                 .y = "SPIRO_Mean", .y.facet = "SPIRO_STD", .title = .tit) %>%
  ggsave("Pics/s04map94.png", plot = ., units = "cm", height = .height, width = .width)
heat_map_facets2(.var = "diversity_mean", .x = "Boundary",
                 .x.facet2 = "Boundary_STD", .x.facet1 = "HK_distribution",
                 .y = "SPIRO_Mean", .y.facet = "SPIRO_STD", .title = .tit) %>%
  ggsave("Pics/s04map95.png", plot = ., units = "cm", height = .height, width = .width)
heat_map_facets2(.var = "diversity_sd", .x = "Boundary",
                 .x.facet2 = "Boundary_STD", .x.facet1 = "HK_distribution",
                 .y = "SPIRO_Mean", .y.facet = "SPIRO_STD", .title = .tit) %>%
  ggsave("Pics/s04map96.png", plot = ., units = "cm", height = .height, width = .width)





# Heat map function:
heat_map = function(.data = tm, .var = "ESBG_mean", .y = "Boundary", .x = "SPIRO_Mean", .title = "") {
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
tm = ts43 %>% #filter(!HK_distribution) %>%
  group_by(Boundary,SPIRO_Mean) %>%
  summarise(across(.cols = diversity:ESBG,
                   list(mean = mean, sd = sd),
                   .names = "{.col}_{.fn}")) %>%
  ungroup()

# Drawing:
.height = 7.5
.width = 29.5
heat_map(.var = "ESBG_mean", .title = "Heat Map 4.3 Step")
heat_map(.var = "ESBG_sd", .title = "Heat Map 4.3 Step")
heat_map(.var = "extremness_mean", .title = "Heat Map 4.3 Step")
heat_map(.var = "extremness_sd", .title = "Heat Map 4.3 Step")
heat_map(.var = "diversity_mean", .title = "Heat Map 4.3 Step")
heat_map(.var = "diversity_sd", .title = "Heat Map 4.3 Step")





heat_map_facet = function(.data = tm, .var = "ESBG_mean", .y = "SPIRO_Mean", .x = "Boundary",
                          .facet = "HK_distribution", .title = "Complex 3-way Heat Map Step 4.3") {
  # .Boundary_STD = eval(str2lang(.facet))
  .data %>%
    ggplot() +
    aes(x = eval(str2lang(.x)), y = eval(str2lang(.y)), col = eval(str2lang(.var)),
        label = round(eval(str2lang(.var)), 2)) +
    facet_wrap(vars(eval(str2lang(.facet))), nrow = 1, labeller = "label_value") + #labeller = "label_both"
    geom_point(shape = 15, size = 9) +
    # geom_text(color = "black", size = 2.5) +
    scale_color_viridis_c() +
    scale_x_continuous(labels = function(x) round(x, 2)) +
    # scale_y_continuous(breaks = seq(0.1, 1, 0.1)) +
    labs(x = .x, y = .y, col = .var, title = paste(.title, .var, sep = ": ")) +
    theme_light() +
    theme(legend.position = "top")
}

# Data preparation on boundary and conformity:
tm = ts43 %>%
  group_by(Boundary, SPIRO_Mean, HK_distribution) %>%
  summarise(across(.cols = diversity:ESBG,
                   list(mean = mean, sd = sd),
                   .names = "{.col}_{.fn}")) %>% ungroup()


# Drawing:
.height = 19.7
.width = 27.5
heat_map_facet(.var = "ESBG_mean")
heat_map_facet(.var = "ESBG_sd")
heat_map_facet(.var = "extremness_mean")
heat_map_facet(.var = "extremness_sd")
heat_map_facet(.var = "diversity_mean")
heat_map_facet(.var = "diversity_sd")
