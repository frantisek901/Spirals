#### Checking effect of population size on effect of boundary







# Head --------------------------------------------------------------------

# clear all
rm(list = ls())

# pacakges
library(tidyverse)


# ReaDING OLD DATA --------------------------------------------------------

load("DATA/Step02data.RData")



# reading new data --------------------------------------------------------


raw = read_csv("DATA/HKidentity_N1000_RS01-10-table.csv", skip = 6)

for (pn in c("11-20", "21-30", "31-40", "41-50", "51-60")) {
  raw = raw %>%
    add_row(read_csv(paste0("DATA/HKidentity_N1000_RS", pn, "-table.csv"), skip = 6))
}

completed = raw %>% count(RS) %>% filter(n == 160) %>% pull(RS)

tb = raw %>%
  filter(RS %in% completed) %>%
  select(HK_distribution = 4, Present_opinion = 5, RS, Use_identity = 12, N = 3, Boundary = 7, Boundary_STD,
         Conformity = 10, Conformity_STD, ticks, diversity = 40, extremness = 41, ESBG= 42) %>%
  mutate(across(c(1, 2, 4), ~factor(.x)))



# Joining data ------------------------------------------------------------

tc = ts20 %>% select(-c(14:18)) %>%
  add_row(tb)



# First fast check --------------------------------------------------------

ts = tc %>%
  group_by(N, Boundary) %>%
  summarise(across(c(ESBG, diversity, extremness), list(mean = mean, sd = sd))) %>%
  ungroup() %>%
  mutate(across(c(N, Boundary), ~factor(.x)))


ts %>%
  ggplot() +
  aes(y = ESBG_mean, x = Boundary, col = N, group = N) +
  geom_point() +
  geom_line() +
  theme_classic()



ts %>%
  ggplot() +
  aes(y = ESBG_sd, x = Boundary, col = N, group = N) +
  geom_point() +
  geom_line() +
  theme_classic()



ts %>%
  ggplot() +
  aes(y = diversity_mean, x = Boundary, col = N, group = N) +
  geom_point() +
  geom_line() +
  theme_classic()



ts %>%
  ggplot() +
  aes(y = diversity_sd, x = Boundary, col = N, group = N) +
  geom_point() +
  geom_line() +
  theme_classic()



ts %>%
  ggplot() +
  aes(y = extremness_mean, x = Boundary, col = N, group = N) +
  geom_point() +
  geom_line() +
  theme_classic()



ts %>%
  ggplot() +
  aes(y = extremness_sd, x = Boundary, col = N, group = N) +
  geom_point() +
  geom_line() +
  theme_classic()



