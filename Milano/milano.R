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
setwd("./Milano/")
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
load("paris.RData")


# Better clearing:
dfx = df %>% select(2:12) %>% unique() #%>%
  # filter(opinions < 4, speaking != 0.95)

dfx2 = dfx %>% mutate(across(c(opinions, id_threshold), ~factor(.x)))

dfx3 = dfx %>% mutate(across(RS:speaking, ~factor(.x)))



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
  aes(x = boundary, y = ESBG, col = id_threshold, group = grp, linetype = opinions) +
  geom_point(size = 5, alpha = 0.35) +
  geom_line() +
  labs(x = "Openness to different opinions", y = "Average polarization (ESBG)", col = "Salience of\nproximity in\nidentity\nrelevant\nopinions\n(SPIRO) :",
       linetype = "Dimensions:", title = "Average polarization caused by\nmain predictors. (N=603,434).") +
  theme_minimal()
ggsave("../Pilsen/milano2.png", height = 4.7, width = 5.7)




# Phase transitions -------------------------------------------------------

## Loading stored data
load("phase2w.RData")


## Firstly, we have to find, what is the highest complete RS, i.e. set of all parameters' combinations simulated:
RS_complete = (phase2w %>% group_by(RS) %>% summarise(n = n()) %>% filter(n == max(n)))$RS
RS_complete = 1:120  # Only 40 seeds are fully completed, let's use also another 80 not fully complete.

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
  rename(`Opinion dimensions` = 1, `Openness to different opinions` = 2, Identity = 3,
         `Salience of distances in identity relevant opinions (SDIRO)` = 4)



## For presenting variability we try now scatter plot only on individual data (non-aggregated):
dfi %>%
  filter(round(100 * id_threshold, 0) %in% seq(5, 85, 20)) %>%
  # sample_n(2000) %>%
  ## Selecting variables:
  select(opinions, boundary, id_threshold, ESBG) %>%

  ## Changing some variables to factors:
  mutate(id_threshold = factor(id_threshold) %>% fct_rev(),
         opinions = factor(opinions)) %>%

  ## Renaming variables according 2022-05-20 meeting/emails:
  rename(`Opinion dimensions` = 1, `Openness to different opinions` = 2,
         `Salience of distances in identity relevant opinions (SDIRO)` = 3) %>%
  # prejmenuj(1:3, c("Opinion dimensions:", "Acceptability of different opinion:",
  #                  "Narrowness of identity group:")) %>%

  ## Graph itself:
  ggplot(aes(x = `Openness to different opinions`, y = ESBG,
             fill = `Salience of distances in identity relevant opinions (SDIRO)`,
             col = `Salience of distances in identity relevant opinions (SDIRO)`,
             group = `Openness to different opinions`)) +
  facet_wrap(vars(`Salience of distances in identity relevant opinions (SDIRO)`, `Opinion dimensions`), ncol=3) +
  geom_point(alpha = 0.10) +
  scale_x_continuous(breaks = seq(0.05, 0.45, 0.10)) +
  labs(title = "'Polarization' in simulations by 'Opinion dimensions' (1, 2, 4),\n'Salience of distances in identity relevant opinions (SDIRO)' (0.05-0.85) and 'Openness to different opinions' (0.05-0.5)",
       x = "Openness of communication norms", y = "Polarization (ESBG)") +
  guides(fill = "none", color = "none") +
  theme_minimal() +
  theme(legend.position = "top")

ggsave("../Milano/milano3.png", width = 9, height = 15, limitsize = FALSE, dpi = 300)



## For presenting variability we try now boxplots on individual data (non-aggregated):
dfi %>%
  filter(round(100 * boundary, 0) %in% seq(5, 45, 5)) %>%
  # sample_n(2000) %>%
  ## Selecting variables:
  select(opinions, boundary, id_threshold, ESBG) %>%

  ## Changing some variables to factors:
  mutate(boundary = factor(boundary), # %>% fct_rev(),
         opinions = factor(opinions)) %>%

  ## Renaming variables according 2022-05-20 meeting/emails:
  rename(`Opinion dimensions` = 1, `Openness to different opinions` = 2, `Salience of identity` = 3) %>%
  # prejmenuj(1:3, c("Opinion dimensions:", "Acceptability of different opinion:",
  #                  "Narrowness of identity group:")) %>%

  ## Graph itself:
  ggplot(
    # aes(fill = `Acceptability of different opinion:`, y = ESBG,
    #          x = `Narrowness of identity group:`,
    #          col = `Acceptability of different opinion:`)
    ) +
  aes(fill = `Openness to different opinions`, x = ESBG,
      y = `Salience of identity`,
      col = `Openness to different opinions`) +
  facet_wrap(vars(`Opinion dimensions`, `Openness to different opinions`), nrow = 3) +
  geom_point(alpha = 0.10) +
  scale_x_continuous(breaks = seq(0.0, 0.8, 0.2)) +
  scale_y_continuous(breaks = seq(0.05, 0.85, 0.2)) +
  labs(title = "'Polarization' in simulations by 'Opinion dimensions' (1, 2, 4), 'Salience of distances in identity relevant opinions (SDIRO)' (0.05-0.85) and\n'Openness to different opinions' (0.05-0.45)",
       x = "Polarization (ESBG)", y = "Salience of distances in identity relevant opinions (SDIRO)") +
  guides(fill = "none", color = "none") +
  theme_minimal() +
  theme(legend.position = "top")

ggsave("../Milano/milano4.png", width = 25, height = 15, limitsize = FALSE, dpi = 300)


dfs %>%
  ggplot() +
  aes(x = `Openness to different opinions`, col = ESBG, fill = ESBG,
      y = `Salience of distances in identity relevant opinions (SDIRO)`) +
  facet_wrap(vars(`Opinion dimensions`), ncol=3) +
  geom_point(alpha = 1, size = 2.5, shape = 22) +
  scale_fill_gradient2(low = "green", mid = "red", high = "black", midpoint = 0.3) +
  scale_color_gradient2(low = "green", mid = "red", high = "black", midpoint = 0.3) +
  scale_y_continuous(breaks = seq(0.05, 0.85, 0.05)) +
  scale_x_continuous(breaks = seq(0.05, 0.50, 0.05)) +
  labs(title = "'Polarization' in simulations by 'Opinion dimensions' (1, 2, 4), 'Salience of distances in identity relevant opinions (SDIRO)' (0.05-0.85) and\n'Openness to different opinions' (0.05-0.5)",
       fill = "Polarization (ESBG)", col = "Polarization (ESBG)") +
  guides(alpha = "none") +
  theme_minimal() +
  theme(legend.position = "top")

ggsave("../Milano/milano5.png", width = 12.5, height = 8)

