#### Script for processing data from Identity paper experiments
#### Now we follow with sensitivity analysis of HK model with
#### global identity and heterogeneous parameters
####


## Encoding: windows-1250
## Created:  2022-11-15 FrK
## Edited:   2022-12-27 FrK

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

## Step 3
# Creating object 'raw' (tibble): Loading....
raw3 = read_csv("GlobalIdentityHK_heterogenousParameters_RS01-05-table.csv", skip = 6) %>%
  add_row(read_csv("GlobalIdentityHK_heterogenousParameters_RS06-10-table.csv", skip = 6)) %>%
  add_row(read_csv("GlobalIdentityHK_heterogenousParameters_RS20.csv", skip = 6)) %>%
  add_row(read_csv("GlobalIdentityHK_heterogenousParameters_RS55.csv", skip = 6))
for (i in seq(11, 56, 5)) {
  raw3 = raw3 %>%
    add_row(read_csv(paste0("GlobalIdentityHK_heterogenousParameters_RS", i, "-", i + 4, "-table.csv"), skip = 6))
}

# Transforming 'raw3' to clean 't3'
t3 = raw3 %>%
  # Selecting and renaming...
  select(HK_distribution = 4, 2, Use_identity = 12, 13:15,
         N = 3, Boundary = 7, 8, Conformity = 10, 11,
         39, diversity = 40, extremness = 41, ESBG = 42,
         boundary_mean = 43, boundary_sd = 44,
         conformity_mean = 45, conformity_sd = 46) %>%
  # Dropping duplicate observations due to experiment runs dubling:
  distinct()


## Step 2
# Creating object 'raw' (tibble): Loading....
raw2 = read_csv("ClassicalHK_heterogenousParameters_RS01-05.csv", skip = 6) %>%
  add_row(read_csv("ClassicalHK_heterogenousParameters_RS06-10.csv", skip = 6))
for (i in seq(11, 56, 5)) {
  raw2 = raw2 %>%
    add_row(read_csv(paste0("ClassicalHK_heterogenousParameters_RS", i, "-", i + 4, ".csv"), skip = 6))
}

# Transforming 'raw2' to clean 't2'
t2 = raw2 %>%
  # Selecting and renaming...
  select(HK_distribution = 4, 2, Use_identity = 12, 13:15,
         N = 3, Boundary = 7, 8, Conformity = 10, 11,
         39, diversity = 40, extremness = 41, ESBG = 42,
         boundary_mean = 43, boundary_sd = 44,
         conformity_mean = 45, conformity_sd = 46) %>%

  # Filtering to only needed cases...
  filter((Conformity == 0.2 | Conformity == 0.8) & (Conformity_STD == 0 | Conformity_STD == 0.1)) %>%

  # Partial processing.
  mutate(
    Identity_Type = "No identity",
    SPIRO_Distribution = "No identity",
    SPIRO_Mean = 0)

## Joining files from both steps together:
ts = t3 %>% add_row(t2) %>%

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
for (i in names(ts)[c(1, 3, 6:11, 20)]) {
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
m1 = lm(ESBG ~ SPIRO_Mean + Boundary + Boundary_STD + Conformity + Conformity_STD + even_N * HK_distribution, ts)
m2 = lm(extremness ~ SPIRO_Mean + Boundary + Boundary_STD + Conformity + Conformity_STD + even_N * HK_distribution, ts)
m3 = lm(diversity ~ SPIRO_Mean + Boundary + Boundary_STD + Conformity + Conformity_STD + even_N * HK_distribution, ts)
stargazer::stargazer(m1, m2, m3, type = "text")

# Full model with interaction of three main factors:
m71 = lm(ESBG ~ SPIRO_Mean * (Boundary) * Boundary_STD + Conformity + Conformity_STD + even_N * HK_distribution, ts)
m72 = lm(extremness ~ SPIRO_Mean * (Boundary) * Boundary_STD + Conformity + Conformity_STD + even_N * HK_distribution, ts)
m73 = lm(diversity ~ SPIRO_Mean * (Boundary) * Boundary_STD + Conformity + Conformity_STD + even_N * HK_distribution, ts)
stargazer::stargazer(m71, m72, m73, type = "text")

# Just two main factors including their interaction:
m81 = lm(ESBG ~ SPIRO_Mean * (Boundary) * (Boundary_STD), ts)
m82 = lm(extremness ~ SPIRO_Mean * (Boundary) * (Boundary_STD), ts)
m83 = lm(diversity ~ SPIRO_Mean * (Boundary) * (Boundary_STD), ts)
stargazer::stargazer(m81, m82, m83, type = "text")

stargazer::stargazer(m1, m71, m81, type = "text", omit.stat = c("f", "ser"))
stargazer::stargazer(m2, m72, m82, type = "text", omit.stat = c("f", "ser"))
stargazer::stargazer(m3, m73, m83, type = "text", omit.stat = c("f", "ser"))
stargazer::stargazer(m1, m71, m81, m2, m72, m82, m3, m73, m83, type = "text", omit.stat = c("f", "ser"))





# Maps ------------------------------------------------------------------

# Heat map function:
heat_map = function(.data = tm, .var = "ESBG_mean", .y = "SPIRO_Mean", .x = "Boundary",
                    .title = "Boundary VS SPIRO") {
  .data %>%
    ggplot() +
    aes(x = eval(str2lang(.x)), y = eval(str2lang(.y)), col = eval(str2lang(.var)),
        label = round(eval(str2lang(.var)), 2)) +
    geom_point(shape = 15, size = 15) +
    geom_text(col = "white") +
    scale_color_viridis_c() +
    #scale_x_continuous(breaks = seq(0.05, 0.35, 0.05)) +
    # scale_y_continuous(breaks = seq(0.1, 1, 0.1)) +
    labs(x = .x, y = .y, col = .var, title = paste(.title, .var, sep = ": ")) +
    theme_light()
}


# Data preparation on Boundary:
tm = ts %>%
  group_by(Boundary, SPIRO_Mean) %>%
  summarise(across(.cols = diversity:ESBG,
                 list(mean = mean, sd = sd),
                 .names = "{.col}_{.fn}")) %>% ungroup()

.height = 10.2
.width = 29.4

# Drawing:
heat_map(.var = "ESBG_sd") %>%
  ggsave("Pics/s03map02.png", plot = ., units = "cm", height = .height, width = .width)
heat_map(.var = "ESBG_mean") %>%
  ggsave("Pics/s03map01.png", plot = ., units = "cm", height = .height, width = .width)
heat_map(.var = "extremness_sd") %>%
  ggsave("Pics/s03map04.png", plot = ., units = "cm", height = .height, width = .width)
heat_map(.var = "extremness_mean") %>%
  ggsave("Pics/s03map03.png", plot = ., units = "cm", height = .height, width = .width)
heat_map(.var = "diversity_sd") %>%
  ggsave("Pics/s03map06.png", plot = ., units = "cm", height = .height, width = .width)
heat_map(.var = "diversity_mean") %>%
  ggsave("Pics/s03map05.png", plot = ., units = "cm", height = .height, width = .width)



heat_map_facets = function(.data = tm, .var = "ESBG_mean", .y = "SPIRO_Mean", .x = "Boundary",
                           .y.facet = "HK_distribution", .x.facet = "N",
                           .title = "Complex Heat Map 3nd Step") {
  .data %>%
    ggplot() +
    aes(x = eval(str2lang(.x)), y = eval(str2lang(.y)), col = eval(str2lang(.var)),
        label = round(eval(str2lang(.var)), 2)) +
    facet_grid(rows = vars(HK_distribution),
               cols = vars(N),
               labeller = "label_both") +
    geom_point(shape = 15, size = 6) +
    geom_text(col = "black", size = 2) +
    scale_color_viridis_c() +
    #scale_x_continuous(breaks = seq(0.05, 0.35, 0.05)) +
    # scale_y_continuous(breaks = seq(0.1, 1, 0.1)) +
    labs(x = .x, y = .y, col = .var, title = paste(.title, .var, sep = ": ")) +
    theme_light()
}

# Data preparation on boundary and conformity:
tm = ts %>%
  group_by(Boundary, SPIRO_Mean, HK_distribution, N) %>%
  summarise(across(.cols = diversity:ESBG,
                   list(mean = mean, sd = sd),
                   .names = "{.col}_{.fn}")) %>% ungroup()
.height = 9.4
.width = 25.6

# Drawing:
heat_map_facets(.var = "ESBG_sd") %>%
  ggsave("Pics/s03map22.png", plot = ., units = "cm", height = .height, width = .width)
heat_map_facets(.var = "ESBG_mean") %>%
  ggsave("Pics/s03map21.png", plot = ., units = "cm", height = .height, width = .width)
heat_map_facets(.var = "extremness_sd") %>%
  ggsave("Pics/s03map24.png", plot = ., units = "cm", height = .height, width = .width)
heat_map_facets(.var = "extremness_mean") %>%
  ggsave("Pics/s03map23.png", plot = ., units = "cm", height = .height, width = .width)
heat_map_facets(.var = "diversity_sd") %>%
  ggsave("Pics/s03map26.png", plot = ., units = "cm", height = .height, width = .width)
heat_map_facets(.var = "diversity_mean") %>%
  ggsave("Pics/s03map25.png", plot = ., units = "cm", height = .height, width = .width)





heat_map_facet = function(.data = tm, .var = "ESBG_mean", .y = "SPIRO_Mean", .x = "Boundary",
                          .facet = "Boundary_STD", .title = "Complex 3-way Heat Map 3nd Step") {
  # .Boundary_STD = eval(str2lang(.facet))
  .data %>%
    ggplot() +
    aes(x = eval(str2lang(.x)), y = eval(str2lang(.y)), col = eval(str2lang(.var)),
        label = round(eval(str2lang(.var)), 2)) +
    facet_wrap(vars(Boundary_STD), ncol = 1, labeller = "label_both") +
    geom_point(shape = 15, size = 6) +
    geom_text(color = "black", size = 1.5) +
    scale_color_viridis_c() +
    #scale_x_continuous(breaks = seq(0.05, 0.35, 0.05)) +
    # scale_y_continuous(breaks = seq(0.1, 1, 0.1)) +
    labs(x = .x, y = .y, col = .var, title = paste(.title, .var, sep = ": ")) +
    theme_light()
}

# Data preparation on boundary and conformity:
tm = ts %>%
  group_by(Boundary, SPIRO_Mean, Boundary_STD) %>%
  summarise(across(.cols = diversity:ESBG,
                   list(mean = mean, sd = sd),
                   .names = "{.col}_{.fn}")) %>% ungroup()
.height = 22.5
.width = 14.4

# Drawing:
heat_map_facet(.var = "ESBG_sd") %>%
  ggsave("Pics/s03map32.png", plot = ., units = "cm", height = .height, width = .width)
heat_map_facet(.var = "ESBG_mean") %>%
  ggsave("Pics/s03map31.png", plot = ., units = "cm", height = .height, width = .width)
heat_map_facet(.var = "extremness_sd") %>%
  ggsave("Pics/s03map34.png", plot = ., units = "cm", height = .height, width = .width)
heat_map_facet(.var = "extremness_mean") %>%
  ggsave("Pics/s03map33.png", plot = ., units = "cm", height = .height, width = .width)
heat_map_facet(.var = "diversity_sd") %>%
  ggsave("Pics/s03map36.png", plot = ., units = "cm", height = .height, width = .width)
heat_map_facet(.var = "diversity_mean") %>%
  ggsave("Pics/s03map35.png", plot = ., units = "cm", height = .height, width = .width)


heat_map_facet = function(.data = tm, .var = "ESBG_mean", .y = "SPIRO_Mean", .x = "Boundary",
                          .facet = "N", .title = "Differences by HK-distribution (T - F)") {
  # .Boundary_STD = eval(str2lang(.facet))
  .data %>%
    ggplot() +
    aes(x = eval(str2lang(.x)), y = eval(str2lang(.y)), col = eval(str2lang(.var)),
        label = round(eval(str2lang(.var)), 2)) +
    facet_wrap(vars(N), ncol = 1, labeller = "label_both") +
    geom_point(shape = 15, size = 6) +
    geom_text(color = "white", size = 1.5) +
    scale_color_viridis_c() +
    #scale_x_continuous(breaks = seq(0.05, 0.35, 0.05)) +
    # scale_y_continuous(breaks = seq(0.1, 1, 0.1)) +
    labs(x = .x, y = .y, col = .var, title = paste(.title, .var, sep = ": ")) +
    theme_light()
}


# Data preparation on boundary and conformity:
tm = ts %>%
  group_by(Boundary, SPIRO_Mean, HK_distribution, N) %>%
  summarise(across(.cols = diversity:ESBG,
                   list(mean = mean, sd = sd),
                   .names = "{.col}_{.fn}")) %>% ungroup() %>%
  # Reshaping to long and then to wide form and then back:
  pivot_longer(diversity_mean:ESBG_sd) %>%
  pivot_wider(id_cols = c(Boundary, SPIRO_Mean, N, name), names_from = HK_distribution, values_from = value) %>%
  mutate(value = `TRUE` - `FALSE`) %>%
  pivot_wider(id_cols = c(Boundary, SPIRO_Mean, N), names_from = name, values_from = value)


.height = 10.1
.width = 14.4

# Drawing:
heat_map_facet(.var = "ESBG_sd") %>%
  ggsave("Pics/s03map12.png", plot = ., units = "cm", height = .height, width = .width)
heat_map_facet(.var = "ESBG_mean") %>%
  ggsave("Pics/s03map11.png", plot = ., units = "cm", height = .height, width = .width)
heat_map_facet(.var = "extremness_sd") %>%
  ggsave("Pics/s03map14.png", plot = ., units = "cm", height = .height, width = .width)
heat_map_facet(.var = "extremness_mean") %>%
  ggsave("Pics/s03map13.png", plot = ., units = "cm", height = .height, width = .width)
heat_map_facet(.var = "diversity_sd") %>%
  ggsave("Pics/s03map16.png", plot = ., units = "cm", height = .height, width = .width)
heat_map_facet(.var = "diversity_mean") %>%
  ggsave("Pics/s03map15.png", plot = ., units = "cm", height = .height, width = .width)



heat_map_facet = function(.data = tm, .var = "ESBG_mean", .y = "SPIRO_Mean", .x = "Boundary",
                          .facet = "HK_distribution", .title = "Differences by N/oddness (`101` - `100`)") {
  # .Boundary_STD = eval(str2lang(.facet))
  .data %>%
    ggplot() +
    aes(x = eval(str2lang(.x)), y = eval(str2lang(.y)), col = eval(str2lang(.var)),
        label = round(eval(str2lang(.var)), 2)) +
    facet_wrap(vars(HK_distribution), ncol = 1, labeller = "label_both") +
    geom_point(shape = 15, size = 6) +
    geom_text(color = "white", size = 1.5) +
    scale_color_viridis_c() +
    #scale_x_continuous(breaks = seq(0.05, 0.35, 0.05)) +
    # scale_y_continuous(breaks = seq(0.1, 1, 0.1)) +
    labs(x = .x, y = .y, col = .var, title = paste(.title, .var, sep = ": ")) +
    theme_light()
}


# Data preparation on boundary and conformity:
tm = ts %>%
  group_by(Boundary, SPIRO_Mean, HK_distribution, N) %>%
  summarise(across(.cols = diversity:ESBG,
                   list(mean = mean, sd = sd),
                   .names = "{.col}_{.fn}")) %>% ungroup() %>%
  # Reshaping to long and then to wide form and then back:
  pivot_longer(diversity_mean:ESBG_sd) %>%
  pivot_wider(id_cols = c(Boundary, SPIRO_Mean, HK_distribution, name), names_from = N, values_from = value) %>%
  mutate(value = `101` - `100`) %>%
  pivot_wider(id_cols = c(Boundary, SPIRO_Mean, HK_distribution), names_from = name, values_from = value)


.height = 10.1
.width = 14.4

# Drawing:
heat_map_facet(.var = "ESBG_sd") %>%
  ggsave("Pics/s03map18.png", plot = ., units = "cm", height = .height, width = .width)
heat_map_facet(.var = "ESBG_mean") %>%
  ggsave("Pics/s03map17.png", plot = ., units = "cm", height = .height, width = .width)
heat_map_facet(.var = "extremness_sd") %>%
  ggsave("Pics/s03map19a.png", plot = ., units = "cm", height = .height, width = .width)
heat_map_facet(.var = "extremness_mean") %>%
  ggsave("Pics/s03map19.png", plot = ., units = "cm", height = .height, width = .width)
heat_map_facet(.var = "diversity_sd") %>%
  ggsave("Pics/s03map19c.png", plot = ., units = "cm", height = .height, width = .width)
heat_map_facet(.var = "diversity_mean") %>%
  ggsave("Pics/s03map19b.png", plot = ., units = "cm", height = .height, width = .width)





