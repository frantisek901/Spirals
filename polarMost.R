#### Script for finding the most polarizing parameters of Spirals' project model

## Encoding: windows-1250
## Created:  2022-03-02 Francesco
## Edited:   2022-03-02 Francesco


## NOTES:
#
# We need to get the most polarizing settings to test then,
# how polarizing it is, when we change parameter's value by +/- 20%.
#



# Head --------------------------------------------------------------------

# Clearing all
rm(list = ls())

# Setting working directory
setwd("./BehaviorSearch/")
getwd()

# Packages
library(dplyr)
library(readr)
library(tidyr)
library(readxl)
library(writexl)
library(ggplot2)
library(knitr)


# My own functon for renaming in Tidyverse
prejmenuj = function(data, positions, new.names) {
  names(data)[positions] = new.names
  data
}


# Loading and cleaning search data -----------------------------------------------------

best = read_csv("eighthSearch.bestHistory.csv") %>%
  prejmenuj(c(1, 47), c("ID", "Rechecked")) %>%
  group_by(ID) %>% filter(Rechecked == max(Rechecked)) %>% ungroup()

full = read_csv("eighthSearch.modelRunHistory.csv") %>%
  prejmenuj(c(1, 50), c("ID", "Final")) %>%
  group_by(ID) %>% filter(Final == max(Final)) %>% ungroup()



# Graphs ------------------------------------------------------------------

# Plotting dependencies of maximum polarization on iteration number
best %>%
  ggplot(aes(x = Rechecked, y = evaluation)) +
  geom_point(size = 3, alpha = 0.2) +
  theme_minimal()

full %>%
  ggplot(aes(x = Final, y = evaluation)) +
  geom_point(size = 3, alpha = 0.2) +
  theme_minimal()

# Checking the dependency of maximal polarization on iteration number via linear regression
lm(Rechecked~evaluation, best) %>% summary()
lm(Final~evaluation, full) %>% summary()
# NOTE: OK, for the best the dependency is there, for the full it is not.



# Finding final parameters ------------------------------------------------

# Getting final parameters of the best runs
final.best = best %>% filter(Rechecked == max(Rechecked)) %>%
  mutate(across(.fns = as.character)) %>%
  pivot_longer(cols = everything(), names_to = "Parameters") %>%
  slice(c(3:5, 7, 9:10, 12, 16:20, 45, 47))

final.full = full %>% filter(Final == max(Final)) %>%
  mutate(across(.fns = as.character)) %>%
  pivot_longer(cols = everything(), names_to = "Parameters") %>%
  slice(c(4:6, 8, 10:11, 13, 17:21, 50))


# Printing the main parameters:
kable(final.best)
kable(final.full)
