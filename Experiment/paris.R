#### Script for producing results for presentation at Paris ICA Conference

## Encoding: windows-1250
## Created:  2022-05-18
## Edited:   2022-05-23


## Notes:
#
#


# Head --------------------------------------------------------------------

# Cleaning
rm(list = ls())

# Setting working directory
setwd("./Experiment/")
getwd()

# Packages
library(dplyr)
library(readr)
library(tidyr)
library(readxl)
library(writexl)
library(rstatix)
library(stargazer)
library(ggplot2)
library(forcats)



# Loading data ------------------------------------------------------------

sel = function(data, pos = c(2, 5, 8, 16, 19, 22, 28, 31, 33, 44)) select(data, runID = 1, pos,  ESBG = 70)
sel1D = function(data, pos = c(2:4, 7:12, 14)) select(data, runID = 1, pos,  ESBG = 65)
pos1Da = c(2:5, 8:12, 14)
sel4D = function(data, pos = c(2:7, 9:12)) select(data, runID = 1, pos,  ESBG = 80)
pos4Da = c(2:5, 8:12, 14)
#sel1D = function(data) names(data)



# Loading 4D data
res4D = read_csv("experiment01part31.csv", skip = 6) %>% sel4D() %>%
  add_row(read_csv("experiment01part31b.csv", skip = 6) %>% sel4D()) %>%
  add_row(read_csv("experiment01part31c.csv", skip = 6) %>% sel4D()) %>%

  add_row(read_csv("experiment01part33.csv", skip = 6) %>% sel4D()) #%>%
  # add_row(read_csv("experiment01part39a.csv", skip = 6) %>% sel4D()) %>%
  # add_row(read_csv("experiment01part39b.csv", skip = 6) %>% sel4D())



# Loading 1D data
res1D = read_csv("experiment01part21.csv", skip = 6) %>% sel1D() %>%
  add_row(read_csv("experiment01part22.csv", skip = 6) %>% sel1D()) %>%
  add_row(read_csv("experiment01part22b.csv", skip = 6) %>% sel1D()) %>%

  add_row(read_csv("experiment01part22b2.csv", skip = 6) %>% sel1D()) %>%
  add_row(read_csv("experiment01part23.csv", skip = 6) %>% sel1D(pos1Da)) %>%
  add_row(read_csv("experiment01part23b.csv", skip = 6) %>% sel1D(pos1Da)) %>%

  add_row(read_csv("experiment01part27a.csv", skip = 6) %>% sel1D()) %>%
  add_row(read_csv("experiment01part27a2.csv", skip = 6) %>% sel1D()) %>%
  add_row(read_csv("experiment01part27b.csv", skip = 6) %>% sel1D()) %>%

  add_row(read_csv("experiment01part27b2.csv", skip = 6) %>% sel1D()) %>%
  add_row(read_csv("experiment01part28a.csv", skip = 6) %>% sel1D(pos1Da)) %>%
  add_row(read_csv("experiment01part28b.csv", skip = 6) %>% sel1D(pos1Da)) # %>%
  # add_row(read_csv("experiment01part29a.csv", skip = 6) %>% sel1D(pos1Da)) %>%
  # add_row(read_csv("experiment01part29b.csv", skip = 6) %>% sel1D(pos1Da))



# Loading 2D data
res2D = read_csv("experiment01part01.csv", skip = 6) %>% sel() %>%
  add_row(read_csv("experiment01part02.csv", skip = 6) %>% sel()) %>%
  add_row(read_csv("experiment01part03.csv", skip = 6) %>% sel()) %>%

  add_row(read_csv("experiment01part04.csv", skip = 6) %>% sel(pos = c(2, 3, 7, 10, 18, 23, 29, 31, 33, 44))) %>%
  add_row(read_csv("experiment01part05.csv", skip = 6) %>% sel(pos = c(2, 3, 7, 10, 18, 23, 29, 31, 33, 44))) %>%
  add_row(read_csv("experiment01part06.csv", skip = 6) %>% sel(pos = c(2, 3, 7, 10, 18, 23, 29, 31, 33, 44))) %>%

  add_row(read_csv("experiment01part11.csv", skip = 6) %>% sel()) %>%
  add_row(read_csv("experiment01part11b.csv", skip = 6) %>% sel()) %>%
  add_row(read_csv("experiment01part11c.csv", skip = 6) %>% sel()) %>%

  add_row(read_csv("experiment01part13.csv", skip = 6) %>% sel()) %>%
  add_row(read_csv("experiment01part14.csv", skip = 6) %>% sel(pos = c(2:4, 8, 11, 19, 24, 30, 32, 44))) %>%
  add_row(read_csv("experiment01part14b.csv", skip = 6) %>% sel(pos = c(2:4, 8, 11, 19, 24, 30, 32, 44))) %>%

  add_row(read_csv("experiment01part15.csv", skip = 6) %>% sel(pos = c(2:4, 8, 11, 19, 24, 30, 32, 44))) %>%
  add_row(read_csv("experiment01part16.csv", skip = 6) %>% sel(pos = c(2:4, 8, 11, 19, 24, 30, 32, 44))) %>%
  add_row(read_csv("experiment01part17a.csv", skip = 6) %>% sel()) %>%

  add_row(read_csv("experiment01part17b.csv", skip = 6) %>% sel()) %>%
  add_row(read_csv("experiment01part17b2.csv", skip = 6) %>% sel()) %>%
  add_row(read_csv("experiment01part18a.csv", skip = 6) %>% sel()) %>%
  add_row(read_csv("experiment01part18b.csv", skip = 6) %>% sel())



df = res1D %>%

  # Joining data:
  add_row(res2D) %>% add_row(res4D) %>%

  # Filtering:
  filter(mode == "openly-listen") %>%

  # Renaming variables:
  rename(rewiring = `p-random`, neighbors = `n-neis`, speaking = `p-speaking-level`,
         identity = `use_identity?`, conformity = `conformity-level`) %>%

  # Changing values:
  mutate(
    id_threshold = if_else(identity, id_threshold, 0),
    # ESBG = ESBG * sqrt(opinions),
    neighbors = 2 * neighbors)

# Saving and loading data:
save(df, file = "paris.RData")
load("paris.RData")


# Better clearing:
dfx = df %>% select(2:12) %>% unique() #%>%
  # filter(opinions < 4, speaking != 0.95)

dfx2 = dfx %>% mutate(across(c(opinions, id_threshold), ~factor(.x)))

dfx3 = dfx %>% mutate(across(RS:speaking, ~factor(.x)))

# Checking data -----------------------------------------------------------

# t = dfx %>%
#   freq_table(vars = c("speaking", "neighbors", "id_threshold", "boundary", "opinions"))



# Regression --------------------------------------------------------------

# Searching for the best model:
m1 = lm(ESBG~opinions+id_threshold+boundary+conformity+neighbors+rewiring+speaking, dfx)
m2 = lm(ESBG~opinions+id_threshold+boundary+conformity+neighbors+rewiring+speaking, dfx2)
m3 = lm(ESBG~opinions+id_threshold+boundary+conformity+neighbors+rewiring+speaking, dfx3)

stargazer(m1, m2, m3, type = "text" #,
          # add.lines = c(BIC(m1), BIC(m2), BIC(m3))
          )

# Best model:
stargazer(m2, type = "text")



# Graph -------------------------------------------------------------------

# Aggregated data:
dfs = dfx3 %>%
  group_by(opinions, id_threshold, boundary) %>%
  summarise(ESBG = mean(ESBG)) %>%
  ungroup() %>%
  mutate(opinions = recode(opinions, `1` = "1 opinion", `2` = "2 opinions", `4` = "4 opinions"),
         grp = paste0(opinions, "_", id_threshold))


# Graphs:
dfs %>%
  ggplot() +
  aes(x = boundary, y = ESBG, fill = id_threshold) +
  facet_grid(rows = vars(opinions)) +
  geom_col(position = position_dodge()) +
  labs(x = "Openness of communication norms", y = "Average polarization", fill = "Salience\nof identity",
       title = "Fig1: Average polarization caused by\nmain predictors. (N=603,434).") +
  theme_minimal()
ggsave("../Paris/paris1.png")

dfs %>%
  ggplot() +
  aes(x = boundary, y = ESBG, col = id_threshold, group = grp, linetype = opinions) +
  geom_point(size = 5, alpha = 0.35) +
  geom_line() +
  labs(x = "Openness of communication norms", y = "Average polarization", col = "Salience\nof identity:",
       linetype = "Dimensions:", title = "Fig1: Average polarization caused by\nmain predictors. (N=603,434).") +
  theme_minimal()
ggsave("../Paris/paris2.png")




# Phase transitions -------------------------------------------------------

## Loading stored data
load("phase2w.RData")


## Firstly, we have to find, what is the highest complete RS, i.e. set of all parameters' combinations simulated:
RS_complete = (phase2w %>% group_by(RS) %>% summarise(n = n()) %>% filter(n == max(n)))$RS
RS_complete = 1:120

## Preparing individual data 'dfi'
dfi = phase2w %>%
  ## Filtering variables:
  filter(RS %in% RS_complete, identity) %>%  #  | (RS >= 61 & RS <= 96),
  ## Denormalizing ESBG:
  mutate(ESBG = ESBG * sqrt(opinions))

## Summarising 'dfi' into 'dfs':
dfs = dfi %>%
  group_by(opinions, boundary, identity, id_threshold) %>%
  summarise(ESBG = mean(ESBG)) %>% ungroup() %>%

  ## Renaming variables according 2022-05-20 meeting/emails:
  rename(`Opinion dimensions` = 1, `Openness of communication norms` = 2, Identity = 3, `Salience of identity` = 4)



## For presenting variability we try now scatter plot only on individual data (non-aggregated):
dfi %>%
  filter(round(100 * id_threshold, 0) %in% seq(5, 85, 10)) %>%
  # sample_n(2000) %>%
  ## Selecting variables:
  select(opinions, boundary, id_threshold, ESBG) %>%

  ## Changing some variables to factors:
  mutate(id_threshold = factor(id_threshold) %>% fct_rev(),
         opinions = factor(opinions)) %>%

  ## Renaming variables according 2022-05-20 meeting/emails:
  rename(`Opinion dimensions` = 1, `Openness of communication norms` = 2, `Salience of identity` = 3) %>%
  # prejmenuj(1:3, c("Opinion dimensions:", "Acceptability of different opinion:",
  #                  "Narrowness of identity group:")) %>%

  ## Graph itself:
  ggplot(aes(x = `Openness of communication norms`, y = ESBG,
             fill = `Salience of identity`,
             col = `Salience of identity`,
             group = `Openness of communication norms`)) +
  facet_wrap(vars(`Salience of identity`, `Opinion dimensions`), ncol=3) +
  geom_point(alpha = 0.10) +
  scale_x_continuous(breaks = seq(0.05, 0.45, 0.10)) +
  labs(title = "'Polarization' in simulations by 'Opinion dimensions' (1, 2, 4),\n'Salience of identity' (0.05-0.85) and 'Openness of communication norms' (0.05-0.5)",
       x = "Openness of communication norms", y = "Polarization (ESBG)") +
  guides(fill = "none", color = "none") +
  theme_minimal() +
  theme(legend.position = "top")

ggsave("../Paris/paris3.png", width = 7.5, height = 25, limitsize = FALSE, dpi = 300)



## For presenting variability we try now boxplots on individual data (non-aggregated):
dfi %>%
  filter(round(100 * boundary, 0) %in% seq(5, 45, 5)) %>%
  # sample_n(2000) %>%
  ## Selecting variables:
  select(opinions, boundary, id_threshold, ESBG) %>%

  ## Changing some variables to factors:
  mutate(boundary = factor(boundary) %>% fct_rev(),
         opinions = factor(opinions)) %>%

  ## Renaming variables according 2022-05-20 meeting/emails:
  rename(`Opinion dimensions` = 1, `Openness of communication norms` = 2, `Salience of identity` = 3) %>%
  # prejmenuj(1:3, c("Opinion dimensions:", "Acceptability of different opinion:",
  #                  "Narrowness of identity group:")) %>%

  ## Graph itself:
  ggplot(
    # aes(fill = `Acceptability of different opinion:`, y = ESBG,
    #          x = `Narrowness of identity group:`,
    #          col = `Acceptability of different opinion:`)
    ) +
  aes(fill = `Openness of communication norms`, x = ESBG,
      y = `Salience of identity`,
      col = `Openness of communication norms`) +
  facet_wrap(vars(`Openness of communication norms`, `Opinion dimensions`), ncol=3) +
  geom_point(alpha = 0.10) +
  scale_x_continuous(breaks = seq(0.0, 0.8, 0.2)) +
  scale_y_continuous(breaks = seq(0.05, 0.85, 0.2)) +
  labs(title = "'Polarization' in simulations by 'Opinion dimensions' (1, 2, 4),\n'Salience of identity' (0.05-0.85) and 'Openness of communication norms' (0.05-0.45)",
       x = "Polarization (ESBG)", y = "Salience of identity") +
  guides(fill = "none", color = "none") +
  theme_minimal() +
  theme(legend.position = "top")

ggsave("../Paris/paris4.png", width = 7.5, height = 25, limitsize = FALSE, dpi = 300)


dfs %>%
  ggplot() +
  aes(x = `Openness of communication norms`, col = ESBG, fill = ESBG,
      y = `Salience of identity`) +
  facet_wrap(vars(`Opinion dimensions`), ncol=3) +
  geom_point(alpha = 1, size = 2.5, shape = 22) +
  scale_fill_gradient2(low = "green", mid = "red", high = "black", midpoint = 0.3) +
  scale_color_gradient2(low = "green", mid = "red", high = "black", midpoint = 0.3) +
  scale_y_continuous(breaks = seq(0.05, 0.85, 0.05)) +
  scale_x_continuous(breaks = seq(0.05, 0.50, 0.05)) +
  labs(title = "'Polarization' in simulations by 'Opinion dimensions' (1, 2, 4), 'Salience of identity' (0.05-0.85) and 'Openness of communication norms' (0.05-0.5)",
       fill = "Polarization (ESBG)", col = "Polarization (ESBG)") +
  guides(alpha = "none") +
  theme_minimal() +
  theme(legend.position = "top")

ggsave("../Paris/paris5.png", width = 12.5, height = 6.5)

