#### ESBG simulation
#### Here we want controllably create different polarized and
#### consensual distributions and measure their ESBG


## Encoding: windows-1250
## Created:  2023-09-27 FranÈesko
## Edited:   2023-09-28 FranÈesko


## NOTES:
##




# Head --------------------------------------------------------------------

# clearing all
rm(list = ls())

# Packages
library(tidyverse)

## Functions
# Calculation of ESBG
esbg100 = function(vec) {
  x = sort(vec)
  l = length(vec)
  x1 = x[1:round(l / 2)]
  x2 = x[(round(l / 2) + 1):l]
  round(abs(mean(x1) - mean(x2)) / (1 + sd(x1) + sd(x2)), 3)
}

# Generation of opinion distribution for 1 or 2 camps
dstr100 = function(concentrated = 0.8, sizeOfBiggerCamp = 0.6, positionOfBigger = 0.1, lngth = 2520) {
  ## Note: 2520 is 'superdivisible' number, divisible by all numbers from 2 to 10 and many more!
  # rands
  rands = round((1 - concentrated) * lngth) %>% runif() %>% round(3)
  # maj
  maj = round(concentrated * sizeOfBiggerCamp * lngth) %>% rep(positionOfBigger, .)
  # mnr
  mnr = round(concentrated * (1 - sizeOfBiggerCamp) * lngth) %>% rep((1 - positionOfBigger), .)
  sort(c(rands, maj, mnr))
}

# Generation of opinion distr for more than 2 camps
dstrN  = function(concentrated = 0.8, N = 3, lngth = 2520) {
  # rands
  rands = round((1 - concentrated) * lngth) %>% runif() %>% round(digits = 3)
  # consts
  grpSz = floor(concentrated * lngth / N)
  posts = seq(0, 1, 1 / (N - 1)) %>% round(digits =  3)
  # final vector
  c(rands, map(posts, ~rep(.x, times = grpSz)) %>% unlist()) %>% sort()
}

# Graph function
grph = function(vec = dstrN(.9, 3), grps = 3, cntr = .9, rat = 0.33, dstn = 0.5, pstn = 0.1) {
  tibble(Opinions = vec) %>%
    ggplot() +
    aes(x = Opinions) +
    geom_density() +
    geom_jitter(aes(y = 0.25), height = .25, alpha = 0.3, color = "blue") +
    labs(y = "Density", title = "Density and parameters of artificial opinion distributions",
         subtitle = paste0("ESBG: ", esbg100(vec), ", N: ", length(vec) , ", Number of camps: ", grps,
                           ", Concentration: ", cntr, ", Majority size: ", rat,
                           ", Camp distance: ", dstn, ", Position: ", pstn)) +
    theme_classic()
}

# Function for creation of data row
drow = function(vec = dstrN(.9, 3), grps = 3, cntr = .9, rat = 0.33, dstn = 0.5, pstn = 0.1) {
  tibble(
    N = length(vec),
    Camps = grps,
    Concentration = cntr,
    Biggest = rat,
    Distance = dstn,
    Position = pstn,
    ESBG = esbg100(vec),
    Diversity = round(sd(vec), 3),
    Manhattan = round(mean(abs(vec - 0.5)), 3))
}



# Creating idealized graphs ---------------------------------------------------------

# Common constants:
commonLength = 1000

# Creating empty tibble
tb = tibble(N = NA_integer_, Camps = NA_integer_, Concentration = NA_integer_,
            Biggest = NA_integer_, Distance = NA_integer_, Position = NA_integer_,
            ESBG = NA_integer_, Diversity = NA_integer_, Manhattan = NA_integer_) %>%
  drop_na()


# Cycle for 1 camp
for (c in c(1, 0.8, 0.6, 0.4, 0.2, 0)) {
  for (p in c(0, 0.1, 0.2, 0.3, 0.4, 0.5)) {
    vctr = dstr100(concentrated = c, sizeOfBiggerCamp = 1, positionOfBigger = p, lngth = commonLength)
    vctr %>% grph(grps = 1, cntr = c, rat = 1, dstn = 0, pstn = p) %>% print()
    tb = tb %>% add_row(drow(vctr, grps = 1, cntr = c, rat = 1, dstn = 0, pstn = p))
  }
}


# Cycle for 2 camps
for (c in c(1, 0.8, 0.6, 0.4, 0.2)) {
  for (p in c(0, 0.1, 0.2, 0.3, 0.4)) {
    for (r in c(.5, .6, .7, .8, .9)) {
      vctr = dstr100(concentrated = c, sizeOfBiggerCamp = r, positionOfBigger = p, lngth = commonLength)
      vctr %>% grph(grps = 2, cntr = c, rat = r, dstn = (1 - (2 * p)), pstn = p) %>% print()
      tb = tb %>% add_row(drow(vctr, grps = 2, cntr = c, rat = r, dstn = (1 - (2 * p)), pstn = p))
    }
  }
}


# Cycle for 3+ camps
for (g in 3:11) {
  for (c in c(1, 0.8, 0.6, 0.4, 0.2)) {
      vctr = dstrN(concentrated = c, N = g, lngth = commonLength)
      vctr %>% grph(grps = g, cntr = c, rat = round(1 / g, 3),
                    dstn = round(1 / (g - 1), 3), pstn = 0) %>% print()
      tb = tb %>% add_row(drow(vctr, grps = g, cntr = c, rat = round(1 / g, 3),
                               dstn = round(1 / (g - 1), 3), pstn = 0))
  }
}



# Graphs ------------------------------------------------------------------

## Graph of resulting measures
# Reshaping data
tbr = tb %>% pivot_longer(ESBG:Manhattan, names_to = "Measure", values_to = "Value") %>%
  mutate(Measure = factor(Measure, levels = c("ESBG", "Diversity", "Manhattan")))

# Drawing graph
tbr %>%
  ggplot() +
  aes(x = Value, color = Measure, fill = Measure) +
  facet_wrap(vars(Measure), scales = "free_x") +
  geom_density(alpha = 0.2) +
  geom_jitter(aes(y = 0.25), height = .25, alpha = 0.2, size = 2) +
  labs(title = "Comparison of different polarization measures", y = "Density") +
  theme_light() +
    theme(legend.position = "bottom")


## Graphs showing relationship between different polarization measures
# ESBG vs. Diversity
tb %>% #filter(Camps < 4) %>%
  ggplot() +
  aes(x = ESBG, y = Diversity) +
  geom_point(alpha = 0.3, color = "steelblue") +
  theme_classic()

# ESBG vs. Manhattan
tb %>% #filter(Camps < 4) %>%
  ggplot() +
  aes(x = ESBG, y = Manhattan) +
  geom_point(alpha = 0.3, color = "steelblue") +
  theme_classic()

# Diversity vs. Manhattan
tb %>% #filter(Camps < 4) %>%
  ggplot() +
  aes(x = Diversity, y = Manhattan) +
  geom_point(alpha = 0.3, color = "steelblue") +
  theme_classic()






