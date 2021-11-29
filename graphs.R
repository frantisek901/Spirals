#### Script for producing graphs for extended ICA abstract we are writting with Mike and Ashley

## Encoding: windows-1250
## Created:  2021-10-25 Francesco
## Edited:   2021-10-30 Francesco



# Head --------------------------------------------------------------------

# We might run data preparation script for sure:
# source("loading.R")

# Clear the memory
rm(list = ls())

# Packages
library(dplyr)
library(forcats)
library(tidyr)
library(readr)
library(readxl)
library(writexl)
library(sjmisc)
library(ggplot2)

# My own functon for renaming in Tidyverse
prejmenuj = function(data, positions, new.names) {
  names(data)[positions] = new.names
  data
}



# Loading data ------------------------------------------------------------

# We load both data files
com = read_xlsx("components.xlsx")  %>%

  # We need many numeric variables in form of factors, so we do so here in one place:
  mutate(neis = neis * 2, across(N_agents:speaking, ~factor(.x)))



sim = read_xlsx("simulations.xlsx") %>%

  # Also in 'sim' we need many numeric variables in form of factors, so we do so here in one place:
  mutate(
    neis = neis * 2,
    random_links = paste0(random_links * 100, "%"),
    across(N_agents:speaking, ~factor(.x)),
    reason = factor(reason) %>% fct_relevel("Stable", after = 0) %>% fct_rev() %>% recode(Stable = "Cohesive"),
    random_links = fct_rev(random_links)
    )



# Preliminary thoughts -------------------------------------------------------------------

## X-variables:
# N_agents
# random_links
# neis
# opinions
# updating
# boundary
# drawn
# speaking

## Y-variables:
# reason (simulation)
# step (simulation)
# size (component)
# opinion.x (component)


## Preliminary thougths:
# I think that the most important results are on the side of whole simulation,
# so I will focus on whole simulation firstly:
#   1) reason -- whether simulation reaches equilibrium, or is fractured, or doesn't reach it
#   2) step -- length of simulation
#   3) size -- I will take only the size of biggest and second biggest components,
#              use them and compare them as proxy for dominance and polarization,
#              but in case of polarization I need also compute opinion distance of these components
#   4) components -- number of components of size 6+
#   5) opinion -- opinion distance between the first and the second component iff there is 2 or more components
#

fpw = 8  # Final picture width
fph = 5.3  # Final picture height


# Effects on 'reason' -----------------------------------------------------

## Network properties -- size, neighbors, re-wiring
ggplot(sim, aes(y = neis, fill = reason)) +
  facet_grid(rows = vars(N_agents), cols = vars(random_links)) +
  geom_bar(position = "fill") +
  labs(
    title = "Effects of network properties on final state of simulation",
    subtitle = "Y-axis defined by Neigborhood size: {8, 32, 128}, panels are defined by\nRe-wiring (X: {5%, 25%}) and Population size (Y: {129, 257, 513}).",
    y = "Neighborhood size", x = "Relative frequency", fill = "Final state:",
    caption = "Fractured state is more or less constant. Low re-wiring and small neighborhood intensify\nthe probability of non-equilibrium state which is also intensified by larger population."
  ) +
  scale_x_continuous(breaks = seq(.2, 1, .2), labels = function(x) paste0(100*x, "%")) +
  theme_minimal() +
  guides(fill = guide_legend(reverse = T)) +
  theme(legend.position = "bottom", legend.direction = "horizontal")
ggsave("Figs/Fig01-PopulationNeighborhoodRewiringOnFinalState.jpg", width = fpw, height = fph)


# Same relationship, but only on frequency of fractured states
sim %>% filter(reason == "Fractured") %>%
  group_by(neis, random_links, N_agents) %>% summarise(frq = n(), lbl = paste0(round(n()/21.6, 1), "%")) %>%
  ggplot(aes(x = neis, y = frq, fill = random_links, label = lbl)) +
  facet_grid(cols = vars(N_agents)) +
  geom_col(position = "dodge") +
  geom_text(position = position_dodge2(width = .9), hjust = 1.5) +
  labs(
    title = "Effects of network properties (population size, neighborhood size and re-wiring)\non final fractured state of simulation",
    subtitle = "Panels are defined by population size (X: {129, 257, 513})",
    x = "Neighborhood size", y = "Frequency (out of 2,160 simulations)", fill = "Re-wiring:",
    caption = "Probability of Fractured state is higher only in the smallest neighborhood (n=8) and is intensified by smaller population and smaller re-wiring."
  ) +
  theme_minimal() +
  guides(fill = guide_legend(reverse = F)) +
  theme(legend.position = "bottom", legend.direction = "horizontal") +
  coord_flip()
ggsave("Figs/Fig02-PopulationNeighborhoodRewiringOnFracturedState.jpg", width = fpw, height = fph)


# Same relationship, but only on frequency of steady states
sim %>% filter(reason == "Stable") %>%
  group_by(neis, random_links, N_agents) %>% summarise(frq = n(), lbl = paste0(round(n()/21.6, 1), "%")) %>%
  ggplot(aes(x = neis, y = frq, fill = random_links, label = lbl)) +
  facet_grid(cols = vars(N_agents)) +
  geom_col(position = "dodge") +
  geom_text(position = position_dodge2(width = .9), hjust = 1.5) +
  labs(
    title = "Effects of network properties (population size, neighborhood size and re-wiring)\non final fractured state of simulation",
    subtitle = "Panels are defined by population size (X: {129, 257, 513})",
    x = "Neighborhood size", y = "Frequency (out of 2,160 simulations)", fill = "Re-wiring:",
    caption = "Probability of Stable state is the higher the bigger neighborhood and is intensified by smaller re-wiring and smaller population."
  ) +
  theme_minimal() +
  guides(fill = guide_legend(reverse = F)) +
  theme(legend.position = "bottom", legend.direction = "horizontal") +
  coord_flip()
ggsave("Figs/Fig03-PopulationNeighborhoodRewiringOnStableState.jpg", width = fpw, height = fph)



# Same relationship, but only on frequency of steady states
sim %>% filter(reason == "Non-equilibrium") %>%
  group_by(neis, random_links, N_agents) %>% summarise(frq = n(), lbl = paste0(round(n()/21.6, 1), "%")) %>%
  ggplot(aes(x = neis, y = frq, fill = random_links, label = lbl)) +
  facet_grid(cols = vars(N_agents)) +
  geom_col(position = "dodge") +
  geom_text(position = position_dodge2(width = .9), hjust = -.05) +
  labs(
    title = "Effects of network properties (population size, neighborhood size and re-wiring)\non final fractured state of simulation",
    subtitle = "Panels are defined by population size (X: {129, 257, 513})",
    x = "Neighborhood size", y = "Frequency (out of 2,160 simulations)", fill = "Re-wiring:",
    caption = "Probability of Non-equilibrium state is the higher the smaller neighborhood and is intensified by higher re-wiring and bigger population."
  ) +
  scale_y_continuous(breaks = seq(0, 800, 200), limits = c(0, 1150))+
  theme_minimal() +
  guides(fill = guide_legend(reverse = F)) +
  theme(legend.position = "bottom", legend.direction = "horizontal") +
  coord_flip()
ggsave("Figs/Fig04-PopulationNeighborhoodRewiringOnNonequilibriumState.jpg", width = fpw, height = fph)



# Boundary detailed on final state ----------------------------------------

## Boundary properties -- method, average size
ggplot(sim, aes(x = boundary, fill = reason)) +
  facet_grid(cols = vars(drawn)) +
  geom_bar(position = position_dodge2(reverse = T), alpha = .75) +
  labs(
    title = "Effects of boundary properties on final state of simulation",
    subtitle = "X-axis defined by average Boundary size: {.1, .2, .3}, \npanels are defined by Assignment method: {constant, uniform}.",
    x = "Boundary size (fraction of maximum theoretical distance)", y = "Relative frequency",
    fill = "Final state:",
    caption = "Stable state is the more probable the wider average Boundary size (with non-linear interaction effect of Assignment method).\nFractured state is more probable with constant Assignment method, intensified with narrower (average) Boundary size.\nNon-equilibrium state is more or less constant with constant Assignment method and\nthe more probable the narrower average Boundary size with uniform Assignment method."
  ) +
  scale_y_continuous(breaks = seq(648, 6480, 648), labels = function(x) paste0(x/64.80, "%")) +
  theme_minimal() +
  guides(fill = guide_legend(reverse = T)) +
  theme(legend.position = "bottom", legend.direction = "horizontal")
ggsave("Figs/Fig05-BoundaryAssignmentOnFinalState.jpg", width = fpw, height = fph)


# speaking
# opinions
# updating
#

## Opinion properties -- number of opinions, number of updated opinions, probability of speaking
ggplot(sim, aes(x = speaking, fill = reason)) +
  facet_grid(cols = vars(opinions), rows = vars(updating)) +
  geom_bar(position = "dodge") +
  labs(
    title = "Effects of opinion properties on final state of simulation",
    subtitle = "X-axis defined by Probability of speaking: {.1, .5, .9, 1}, \npanels are defined by Opinion dimensions {X: 1, 2, 16} and Dimensions updated at once {Y: 1, 8, 16}.",
    x = "Probability of speaking", y = "Relative frequency", fill = "Final state:",
    caption = "Fractured state is the more probable the more Opinion dimensions, no effect of Opinions updated at once and Probability of speaking.\nStable state is the more probable the less Opinion dimensions and the lower Probability of speaking and the less Dimensions updated at once.\n Non-equilibrium state is the more probable the more Opinion dimensions and \nthe higher Probability of speaking and the more Dimensions updated at once."
  ) +
  scale_y_continuous(breaks = seq(108, 1080, 432), labels = function(x) paste0(x/10.80, "%")) +
  theme_minimal() +
  guides(fill = guide_legend(reverse = F)) +
  theme(legend.position = "bottom", legend.direction = "horizontal")
ggsave("Figs/Fig06-SpeakingDimensionsUpdatingOnFinalState.jpg", width = fpw, height = fph)



ggplot(filter(sim, opinions == 16), aes(x = speaking, fill = reason)) +
  facet_grid(cols = vars(opinions), rows = vars(updating)) +
  geom_bar(position = "dodge") +
  labs(
    title = "Effects of opinion properties on final state of simulation",
    subtitle = "X-axis defined by Probability of speaking: {.1, .5, .9, 1}, \npanels are defined by Opinion dimensions {X: 16} and Dimensions updated at once {Y: 1, 8, 16}.",
    x = "Probability of speaking", y = "Relative frequency", fill = "Final state:",
    caption = "Fractured state probability is 50% with 16 Opinion dimensions.\n Non-linear interaction effect of Opinions updated at once and Probability of speaking on probability of Stable and Non-equilibrium states."
  ) +
  scale_y_continuous(breaks = seq(108, 1080, 216), labels = function(x) paste0(x/10.80, "%")) +
  theme_minimal() +
  guides(fill = guide_legend(reverse = F)) +
  theme(legend.position = "bottom", legend.direction = "horizontal")
ggsave("Figs/Fig07-SpeakingDimensionsUpdatingOnFinalState.jpg", width = fpw, height = fph)




# Opinion detailed on final state -----------------------------------------

## Opinion properties -- number of opinions, number of updated opinions, probability of speaking
sim %>% mutate(
  boundary_type =
    if_else(((boundary == .2 | boundary == .1) & drawn == "constant") | (boundary == 0.1 & drawn == "uniform"),
            "Discordant", "Consensual")) %>%
  ggplot(aes(x = speaking, fill = reason)) +
  facet_grid(cols = vars(opinions), rows = vars(boundary_type)) +
  geom_bar(position = "dodge", alpha = .75) +
  labs(
    title = "Effects of opinion properties and boundary type on final state of simulation",
    subtitle = "X-axis defined by Probability of speaking: {.1, .5, .9, 1}, \npanels are defined by Opinion dimensions {X: 1, 2, 16} and Boundary type {Y: Consensual, Discordant}.",
    x = "Probability of speaking", y = "Relative frequency", fill = "Final state:",
    caption = "Fractured state is the more probable with discordant Boundary type and the more Opinion dimensions, no effect of Probability of speaking.\nStable state is the more probable with consensual Boundary type and the less Opinion dimensions and the lower Probability of speaking.\n Non-equilibrium state is the more probable the more Opinion dimensions and the higher Probability of speaking and \nwith discordant Boundary type (only with low number of Opinion dimensions)."
  ) +
  scale_y_continuous(breaks = seq(162, 1620, 324), labels = function(x) paste0(x/16.20, "%")) +
  theme_minimal() +
  guides(fill = guide_legend(reverse = F)) +
  theme(legend.position = "bottom", legend.direction = "horizontal")
ggsave("Figs/Fig08-SpeakingDimensionsBoundaryOnFinalState.jpg", width = fpw, height = fph)



# Network on final state --------------------------------------------------

## Network properties -- relative size of neighborhoods, re-wiring
sim %>% mutate(size = round(as.numeric(as.character(neis)) /
                              (as.numeric(as.character(N_agents)) - 1) * 100, 1) %>% factor()) %>%
  ggplot(aes(y = fct_rev(size), fill = reason)) +
  facet_grid(cols = vars(random_links)) +
  geom_bar(position = position_fill(reverse = T), alpha = 0.75) +
  labs(
    title = "Effects of network properties on final state of simulation",
    subtitle = "Y-axis defined by Relative neighborhood size: {1.6%...100%}, panels are defined by Re-wiring: {5%, 25%}.",
    y = "Relative neighborhood size (order reversed)", x = "Relative frequency", fill = "Final state:",
    caption = "Fractured state is more or less constant. Low Re-wiring and the smaller Relative neighborhood intensify the probability of Non-equilibrium state."
  ) +
  scale_x_continuous(breaks = seq(.2, 1, .2), labels = function(x) paste0(100*x, "%")) +
  scale_y_discrete(labels = function(x) paste0(x, "%")) +
  theme_minimal() +
  guides(fill = guide_legend(reverse = T)) +
  theme(legend.position = "bottom", legend.direction = "horizontal") +
  coord_flip()
ggsave("Figs/Fig09-RelativeneighborhoodRewiringOnFinalState.jpg", width = fpw, height = fph)


## Network properties -- relative size of neighborhoods, re-wiring
sim %>% mutate(size = round(as.numeric(as.character(neis)) /
                              (as.numeric(as.character(N_agents)) - 1) * 100, 1) %>% factor(),
               speaking = if_else(speaking == .9 | speaking == 1, "High", "Low"),
               boundary_type =
                 if_else(((boundary == .2 | boundary == .1) & drawn == "constant") | (boundary == 0.1 & drawn == "uniform"),
                         "Discordant", "Consensual")) %>%
  ggplot(aes(y = size, fill = reason)) +
  facet_grid(cols = vars(random_links), rows = vars(fct_rev(speaking))) +
  geom_bar(position = "fill", alpha = 0.75) +
  labs(
    title = "Effects of network properties and Probability of speaking on final state of simulation",
    subtitle = "Y-axis defined by Relative neighborhood size: {1.6%...100%}, panels are defined by \nRe-wiring (X: {5%, 25%}) and  Probability of speaking (Y: Low = {.1, .5}, High = {.9, 1}).",
    y = "Relative neighborhood size", x = "Relative frequency", fill = "Final state:",
    caption = "Fractured state is more or less constant. The lower Re-wiring and the smaller Relative neighborhood and\n the higher Probability of speaking then the higher probability of Non-equilibrium state."
  ) +
  scale_x_continuous(breaks = seq(.2, 1, .2), labels = function(x) paste0(100*x, "%")) +
  scale_y_discrete(labels = function(x) paste0(x, "%")) +
  theme_minimal() +
  guides(fill = guide_legend(reverse = T)) +
  theme(legend.position = "bottom", legend.direction = "horizontal")
ggsave("Figs/Fig10-RelativeneighborhoodRewiringSpeakingOnFinalState.jpg", width = fpw, height = fph)


## Network properties -- relative size of neighborhoods, re-wiring
sim %>% mutate(size = round(as.numeric(as.character(neis)) /
                              (as.numeric(as.character(N_agents)) - 1) * 100, 1) %>% factor(),
               opinions = if_else(opinions == 16, "Many", "Few"),
               boundary_type =
                 if_else(((boundary == .2 | boundary == .1) & drawn == "constant") | (boundary == 0.1 & drawn == "uniform"),
                         "Discordant", "Consensual")) %>%
  ggplot(aes(y = size, fill = reason)) +
  facet_grid(cols = vars(random_links), rows = vars((opinions))) +
  geom_bar(position = "fill", alpha = 0.75) +
  labs(
    title = "Effects of network properties and Opinion dimensions on final state of simulation",
    subtitle = "Y-axis defined by Relative neighborhood size: {1.6%...100%}, \npanels are defined by Re-wiring (X: {5%, 25%}) and Opinion dimensions (Y: Few = {1, 2}, Many = 16).",
    y = "Relative neighborhood size", x = "Relative frequency", fill = "Final state:",
    caption = "The more Opinion dimensions the more probable is the  Fractured state. The higher Re-wiring and the wider Relative neighborhood and\n the fewer Opinion dimensions then the higher probability of Stable state."
  ) +
  scale_x_continuous(breaks = seq(.2, 1, .2), labels = function(x) paste0(100*x, "%")) +
  scale_y_discrete(labels = function(x) paste0(x, "%")) +
  theme_minimal() +
  guides(fill = guide_legend(reverse = T)) +
  theme(legend.position = "bottom", legend.direction = "horizontal")
ggsave("Figs/Fig11-RelativeneighborhoodRewiringOpinionsOnFinalState.jpg", width = fpw, height = fph)


## Network properties -- relative size of neighborhoods, re-wiring
sim %>% mutate(size = round(as.numeric(as.character(neis)) /
                              (as.numeric(as.character(N_agents)) - 1) * 100, 1) %>% factor(),
               boundary_type =
                 if_else(((boundary == .2 | boundary == .1) & drawn == "constant") | (boundary == 0.1 & drawn == "uniform"),
                         "Discordant", "Consensual")) %>%
  ggplot(aes(y = size, fill = reason)) +
  facet_grid(cols = vars(random_links), rows = vars(boundary_type)) +
  geom_bar(position = "fill", alpha = 0.75) +
  labs(
    title = "Effects of network properties and Boundary type on final state of simulation",
    subtitle = "Y-axis defined by Relative neighborhood size: {1.6%...100%}, \npanels are defined by Re-wiring (X: {5%, 25%}) and Boundary type (Y: {Consensual, Discordant}).",
    y = "Relative neighborhood size", x = "Relative frequency", fill = "Final state:",
    caption = "The Fractured state is more probable with the discordant Boundary type. The higher Re-wiring and the wider Relative neighborhood and\n with the consensual Boundary type then the higher probability of Stable state."
  ) +
  scale_x_continuous(breaks = seq(.2, 1, .2), labels = function(x) paste0(100*x, "%")) +
  scale_y_discrete(labels = function(x) paste0(x, "%")) +
  theme_minimal() +
  guides(fill = guide_legend(reverse = T)) +
  theme(legend.position = "bottom", legend.direction = "horizontal")
ggsave("Figs/Fig12-RelativeneighborhoodRewiringBoundaryOnFinalState.jpg", width = fpw, height = fph)


# Opinion on final state --------------------------------------------------

sim %>% mutate(
  boundary_type =
    if_else(((boundary == .2 | boundary == .1) & drawn == "constant") | (boundary == 0.1 & drawn == "uniform"),
            "Discordant boundary", "Consensual boundary"),
  opinions = if_else(opinions == 16, "Many dimensions", "Few dimensions"),
  speaking = if_else(speaking == .9 | speaking == 1, "High", "Low")) %>%
  group_by(speaking, opinions, boundary_type, reason) %>%
  summarise(n = n()) %>%
  group_by(speaking, opinions, boundary_type) %>%
  mutate(N = sum(n), p = round(n / N * 100, 1)) %>%
  ggplot(aes(x = fct_rev(speaking), y = p, fill = reason)) +
  facet_grid(cols = vars(opinions, boundary_type)) +
  geom_col(position = position_dodge2(padding = 0.1, preserve = "single", reverse = T), alpha = .75, width = .9) +
  labs(
    title = "Effects of opinion properties and boundary type on final state of simulation",
    subtitle = "X-axis defined by Probability of speaking (Low = {0.1, 0.5}, High = {0.9, 1.0}), panels are defined by \nOpinion dimensions (Few = {1, 2}, Many = 16) and Boundary type (Consensual, Discordant).",
    x = "Probability of speaking", y = "Relative frequency", fill = "Final state:",
    caption = "Fractured state is the more probable with discordant Boundary type and the more Opinion dimensions, no effect of Probability of speaking.\nStable state is the more probable with consensual Boundary type and the less Opinion dimensions and the lower Probability of speaking.\n Non-equilibrium state is the more probable the more Opinion dimensions and the higher Probability of speaking and \nwith discordant Boundary type (only with few Opinion dimensions, in many dimensions Non-equilibrium state turns into Fractured state)."
  ) +
  scale_y_continuous(breaks = seq(10, 90, 20), labels = function(x) paste0(x, "%")) +
  theme_minimal() +
  guides(fill = guide_legend(reverse = T)) +
  theme(legend.position = "bottom", legend.direction = "horizontal")
ggsave("Figs/Fig13-SpeakingDimensionsBoundaryOnFinalState.jpg", width = fpw, height = fph)



# Master graph: final state -----------------------------------------------

sim %>% mutate(
  size = round(as.numeric(as.character(neis)) /
                 (as.numeric(as.character(N_agents)) - 1) * 100, 0) %>% factor(),
  boundary_type =
    if_else(((boundary == .2 | boundary == .1) & drawn == "constant") | (boundary == 0.1 & drawn == "uniform"),
            "Close system", "Open system") %>% factor() %>% fct_rev(),
  opinions = if_else(opinions == 16, "Many dimensions", "Few dimensions"),
  random_links = if_else(random_links == "5%", "More F2F", "Less F2F"), # %>% factor() %>% fct_rev(),
  speaking = if_else(speaking == .9 | speaking == 1, "High spkn'", "Low spkn'")) %>%
  group_by(speaking, opinions, boundary_type, size, random_links, reason) %>%
  summarise(n = n()) %>%
  group_by(speaking, opinions, boundary_type, size, random_links) %>%
  mutate(N = sum(n), p = round(n / N * 100, 1)) %>% ungroup %>% ungroup %>%
  # frq (N) %>%
  ggplot(aes(x = fct_rev(size), y = p, fill = reason)) +
  facet_grid(cols = vars(opinions, boundary_type), rows = vars((speaking), fct_rev(random_links))) +
  # geom_col(position = position_dodge2(padding = 0.1, preserve = "single", reverse = T), alpha = .75, width = .9) +
  geom_col(position = position_fill(reverse = T), alpha = .75, width = .9) +
  labs(
    title = "Effects of opinion and system properties on final state of simulation",
    subtitle = "X-axis defined by Relative size of neighborhood, panels are defined by \nX: Opinion dimensions (Few = {1, 2}, Many = 16) and System type (Open, Close),\nY: Out-spokeness (Low = {0.1, 0.5}, High = {0.9, 1.0}) and Interpersonal communication (More F2F, Less F2F).",
    x = "Relative size of neighborhood (%, order reversed)", y = "Relative frequency", fill = "Final state:",
    caption = "Fractured state is almost exclusive with the discordant Boundary type and many Opinion dimensions.\nStable state is more probable with consensual Boundary, few Opinion dimensions, low Probability of speaking and higher Re-wiring."
  ) +
  scale_y_continuous(breaks = seq(10, 90, 20), labels = function(x) paste0(x, "%")) +
  # scale_x_discrete(labels = function(x) paste0(x, "%")) +
  theme_minimal() +
  guides(fill = guide_legend(reverse = T)) +
  theme(legend.position = "bottom", legend.direction = "horizontal")
ggsave("Figs/Fig14-NetworksAndOpinionsOnFinalState.jpg",
       # width = fpw, height = fph
       width = 9, height = 6
       )

#












