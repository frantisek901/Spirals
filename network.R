#### Reading in network (and other) data recorded in NetLogo and saving them back in cleaner form

## Encoding: windows-1250
## Created:  2021-11-19 Francesco
## Edited:   2022-01-02 Francesco


## NOTES:
#  Data production in NetLogo is changed, but not the script here! Look out!
#


# Head --------------------------------------------------------------------

# Clear the memory
rm(list = ls())

# Packages
library(dplyr)
library(tidyr)
library(readr)
library(stringr)
library(readxl)
library(writexl)
library(sjmisc)
library(forcats)
library(ggplot2)

# My own functon for renaming in Tidyverse
prejmenuj = function(data, positions, new.names) {
  names(data)[positions] = new.names
  data
}



# Reading in and storing data from experiment -----------------------------------------

d = tibble(inFile = dir(path = "Sims")) %>%  # Turning filenames into a tibble
  slice(-(1:3)) %>%  # Slicing out the first three rows containing "strange" files -- check it out in the final code!
  separate(inFile, sep = "_", remove = F, into = c(  # We separate meta info into distinctive variables:
    "Type", "Seed","Population", "Random_links","Close_links", "Dimensions",
    "Updating", "Boundary","Boundary_method", "P_speaking","P_speaking_method", "Mode", "Steps")) %>%
  mutate(across(.cols = c(Seed:Boundary, P_speaking, Steps), ~ parse_number(.x))) %>%  # We parse respective vars as numbers
  group_by(Seed, Population, Random_links, Close_links, Dimensions,
           Boundary, Boundary_method, P_speaking,  P_speaking_method, Mode) %>%
  mutate(Sim_ID = cur_group_id(),
         outFile = paste0("RData/", str_replace(inFile, ".csv", ".RData")),
         inFile = paste0("Sims/", inFile)) %>%
  relocate(c(outFile, Sim_ID), .before = Type) %>%
  arrange(Sim_ID, Steps, Type)

save(d, file = "dataCatalogue.RData")



# The first graphs --------------------------------------------------------

d %>% filter(Steps == max(Steps), Type == "Nodes01") %>% ungroup() %>%
  mutate(
    Finished = if_else(Steps < 1000, "Yes", "No") %>%  factor(levels = c("Yes", "No")),
    across(Population:Mode, ~factor(.x))) %>%
  group_by(P_speaking, Boundary_method, Mode, Population, Boundary, Random_links) %>%
  mutate(N = n()) %>%
  group_by(Finished, P_speaking, Boundary_method, Mode, Population, Boundary, Random_links, N) %>%
  summarise(perc = n()) %>% mutate(perc = round(perc / N * 100, 1)) %>%
  ggplot(aes(y = perc, x = Finished, fill = Finished)) +
  facet_grid(
    # Set according model: glm(Finished ~ P_speaking + Boundary_method + Mode + Population + Random_links + Boundary + P_speaking_method + Close_links + Dimensions, family = binomial, data = d)
    cols = vars(fct_rev(P_speaking), Boundary_method, Boundary, Population#, fct_rev(Dimensions)
    ),
    rows = vars(fct_rev(Mode), fct_rev(Random_links)#, P_speaking_method, Close_links
    )) +
  geom_col() +
  theme_minimal()

ggsave("newExperimentResult.jpg", width = 12, height = 6)



# Moving data from *.csv to *.RData ---------------------------------------

## For sure, loading data catalogue in:
load("dataCatalogue.RData")

## Iterating over the catalogue and transforming all the datafiles:
for (i in 1:nrow(d)) {
  dfi = read_csv(d[[i, "inFile"]])
  save(dfi, file = d[[i, "outFile"]])
}




# Combining data files from one simulation into one file ------------------------

# We define firstly the function for combining:
combine_data = function(ID) {  # We need only to put in ID number of the simulation from data catalogue
  # We load datacatologue and filter out only data on chosen simulation.
  load("dataCatalogue.RData")
  df = filter(d, Sim_ID == ID)

  # We load first file on links and store it as 'links'
  load(df[[1, "outFile"]])
  links = dfi

  # We load second file on nodes and store it as 'nodes'
  load(df[[2, "outFile"]])
  nodes = dfi

  # Combining data on nodes with the data on links/distances
  if (df[[1, "Dimensions"]] == 1) {
    dc = left_join(links, nodes, by = c("ID1" = "ID")) %>%
      rename("Uncertainty_1" = Uncertainty, "pSpeaking_1" = pSpeaking, "Speaks_1" = Speaks, "Opinion1_1" = Opinion1) %>%
      left_join(nodes, by = c("ID2" = "ID")) %>%
      rename("Uncertainty_2" = Uncertainty, "pSpeaking_2" = pSpeaking, "Speaks_2" = Speaks, "Opinion1_2" = Opinion1)
  } else {
    dc = left_join(links, nodes, by = c("ID1" = "ID")) %>%
      rename("Uncertainty_1" = Uncertainty, "pSpeaking_1" = pSpeaking, "Speaks_1" = Speaks, "Opinion1_1" = Opinion1,
             "Opinion2_1" = Opinion2, "Opinion3_1" = Opinion3, "Opinion4_1" = Opinion4, "Opinion5_1" = Opinion5,
             "Opinion6_1" = Opinion6, "Opinion7_1" = Opinion7, "Opinion8_1" = Opinion8) %>%
      left_join(nodes, by = c("ID2" = "ID")) %>%
      rename("Uncertainty_2" = Uncertainty, "pSpeaking_2" = pSpeaking, "Speaks_2" = Speaks, "Opinion1_2" = Opinion1,
             "Opinion2_2" = Opinion2, "Opinion3_2" = Opinion3, "Opinion4_2" = Opinion4, "Opinion5_2" = Opinion5,
             "Opinion6_2" = Opinion6, "Opinion7_2" = Opinion7, "Opinion8_2" = Opinion8)
  }

  # Adding 'Steps' variable:
  dc$Steps = df[[1, "Steps"]]



  ## Resting files we load in and add via cycle:
  for (i in 2:(nrow(df)/2)) {
    # We load first file on links and store it as 'links'
    load(df[[(2 * i) - 1, "outFile"]])
    links = dfi

    # We load second file on nodes and store it as 'nodes'
    load(df[[(2 * i), "outFile"]])
    nodes = dfi

    # Combining data on nodes with the data on links/distances
    if (df[[1, "Dimensions"]] == 1) {
      dcx = left_join(links, nodes, by = c("ID1" = "ID")) %>%
        rename("Uncertainty_1" = Uncertainty, "pSpeaking_1" = pSpeaking, "Speaks_1" = Speaks, "Opinion1_1" = Opinion1) %>%
        left_join(nodes, by = c("ID2" = "ID")) %>%
        rename("Uncertainty_2" = Uncertainty, "pSpeaking_2" = pSpeaking, "Speaks_2" = Speaks, "Opinion1_2" = Opinion1)
    } else {
      dcx = left_join(links, nodes, by = c("ID1" = "ID")) %>%
        rename("Uncertainty_1" = Uncertainty, "pSpeaking_1" = pSpeaking, "Speaks_1" = Speaks, "Opinion1_1" = Opinion1,
               "Opinion2_1" = Opinion2, "Opinion3_1" = Opinion3, "Opinion4_1" = Opinion4, "Opinion5_1" = Opinion5,
               "Opinion6_1" = Opinion6, "Opinion7_1" = Opinion7, "Opinion8_1" = Opinion8) %>%
        left_join(nodes, by = c("ID2" = "ID")) %>%
        rename("Uncertainty_2" = Uncertainty, "pSpeaking_2" = pSpeaking, "Speaks_2" = Speaks, "Opinion1_2" = Opinion1,
               "Opinion2_2" = Opinion2, "Opinion3_2" = Opinion3, "Opinion4_2" = Opinion4, "Opinion5_2" = Opinion5,
               "Opinion6_2" = Opinion6, "Opinion7_2" = Opinion7, "Opinion8_2" = Opinion8)
    }

    # Adding 'Steps' variable:
    dcx$Steps = df[[(2 * i), "Steps"]]

    # Combining 'dc' and 'dcx' into new 'dc', i.e. adding right processed data:
    dc = dc %>% add_row(dcx)
  }

  # Returning combined file back:
  dc
}


# Testing of function 'combine_data()'
dc = combine_data(19)
save(dc, file = "RData/dc19.RData")

# Using the function on the largest example (1001 nodes, 8 opinion dimensions, 6 time points)
dy = d %>% filter(Dimensions == 8, Population == 1001, Close_links == 50, Steps > 800)
dc = combine_data(205)
save(dc, file = "RData/SimID205_101_1001_0_50_8_1_0.1_uniform_0.5_uniform_vaguely-speak.RData")

dy = d %>% filter(Dimensions == 8, Population == 101, Close_links == 50, Steps > 800)
dc = combine_data(348)
save(dc, file = "RData/dc348.RData")



# ====== CODE is updated up to here, code bellow is still old =====

# Measuring entropy of distances ------------------------------------------

# 1D ----------------------------------------------------------------------

entropy1D = function(v1, bin.width = 0.1) {
  # preparation of empty vector of distances
  l = length(v1)  # Length of file/number of agents in simulation
  ld = (l * (l - 1)) / 2  # Length of 'd' is number of pair of agents, according equation: N * (N - 1) / 2
  d = rep(NA, ld)  # We construct empty vector of length 'ld'
  c = 0  # Counter which will be updated inside of the following cycles so to be able points to the first empty cell in 'd'.
  # Note: we start 'c' with 0 because we will start with updating 'c'.

  # We need two for cycles for going over all pairs of agents and measures their distance in 2D space
  for (i in 1:(l - 1)) {  # We start first cycle with the first agent and end with the next to the last.
    for (j in (i + 1):l) {
      c = c + 1  # Updating counter
      d[c] = abs(v1[i] - v1[j])  # Euclidean distance
    }
  }

  # Now we use histogram function to separate vector into bins and
  # transform it into counts of these bins.
  cnts = hist(d, breaks = seq(0, 2.11, bin.width))$counts  # Transformation
  cnts = cnts[cnts>0]  # Filtering out 0s.

  # Now we transform counts into relative freqencies:
  frqs = cnts / ld

  # Now we transform them into product of relative frequency and its logarithm:
  ents = frqs * log(frqs, 2)

  # Now we need the value of maximum entropy with given granularity:
  bin.num = ceiling(2 / bin.width)
  max.ent = log(1 / bin.num, 2)

  # Here in clean environment we prepare all coefficients:
  p_ent =  round(100 - (sum(ents) / (max.ent) * 100), 1)
  p_oneGroup = round(sqrt(cnts[1] * 2) / l * 100, 1)
  n_groups = round(l * l / (l + (2 * cnts[1])), 2)
  p_lengths =  round(frqs[1] * 100, 1)
  SD = round(sd(d), 3)
  DI = round(mean(d) / 2 * 100, 3)

  # Returning back four parameters as one string:
  #print(paste(p_ent, p_oneGroup, n_groups, p_lengths, SD, sep = "_"))
  paste(p_ent, p_oneGroup, n_groups, p_lengths, SD, DI, sep = "_")
}

# dfp = dfb %>% filter(Dimensions == 1, Sim_ID <= 4) %>% group_by(Sim_ID)
# de = entropy1D(dfp$Opinion1_Final, 0.1)
# de



a = Sys.time()
dfy = dfb %>%
  # filter(Sim_ID <= 10) %>%
  filter(Dimensions == 1) %>%
  group_by(Sim_ID, Seed, Population, Random_links, Close_links, Dimensions, Boundary, Boundary_method, P_speaking,  Mode) %>%
  summarise(ent = entropy1D(Opinion1_Final, 0.1)) %>%
  separate(ent,
           into = c("Far_from_entropy", "One_group_size",
                    "Number_of_equal_groups", "Zero_lenghts", "SD", "DI"),
           sep = "_", convert = T)
b = Sys.time()
b - a

# Saving processed meta indicators:
write_csv(dfy, "Sims02_processed_1D.csv")



# 2D ----------------------------------------------------------------------

entropy2D = function(v1, v2, bin.width = 0.1) {
  # preparation of empty vector of distances
  l = length(v1)  # Length of file/number of agents in simulation
  ld = (l * (l - 1)) / 2  # Length of 'd' is number of pair of agents, according equation: N * (N - 1) / 2
  d = rep(NA, ld)  # We construct empty vector of length 'ld'
  c = 0  # Counter which will be updated inside of the following cycles so to be able points to the first empty cell in 'd'.
  # Note: we start 'c' with 0 because we will start with updating 'c'.

  # We need two for cycles for going over all pairs of agents and measures their distance in 2D space
  for (i in 1:(l - 1)) {  # We start first cycle with the first agent and end with the next to the last.
    for (j in (i + 1):l) {
      c = c + 1  # Updating counter
      d1 = v1[i] - v1[j]  # distance in D1
      d2 = v2[i] - v2[j]  # distance in D2
      d[c] = sqrt((d1 ^ 2) + (d2 ^ 2))  # Euclidean distance
    }
  }

  # Now we use histogram function to separate vector into bins and
  # transform it into counts of these bins.
  cnts = hist(d, breaks = seq(0, 3, bin.width))$counts  # Transformation
  cnts = cnts[cnts>0]  # Filtering out 0s.

  # Now we transform counts into relative freqencies:
  frqs = cnts / ld

  # Now we transform them into product of relative frequency and its logarithm:
  ents = frqs * log(frqs, 2)

  # Now we need the value of maximum entropy with given granularity:
  bin.num = ceiling(2.82 / bin.width)
  max.ent = log(1 / bin.num, 2)

  # Here in clean environment we prepare all coefficients:
  p_ent =  round(100 - (sum(ents) / (max.ent) * 100), 1)
  p_oneGroup = round(sqrt(cnts[1] * 2) / l * 100, 1)
  n_groups = round(l * l / (l + (2 * cnts[1])), 2)
  p_lengths =  round(frqs[1] * 100, 1)
  SD = round(sd(d), 3)
  DI = round(mean(d) / sqrt(8) * 100, 3)

  # Returning back four parameters as one string:
  #print(paste(p_ent, p_oneGroup, n_groups, p_lengths, SD, sep = "_"))
  paste(p_ent, p_oneGroup, n_groups, p_lengths, SD, DI, sep = "_")
}

a = Sys.time()
dfy = dfb %>%
  # filter(Sim_ID <= 72) %>%
  filter(Dimensions == 2) %>%
  group_by(Sim_ID, Seed, Population, Random_links, Close_links, Dimensions, Boundary, Boundary_method, P_speaking,  Mode) %>%
  summarise(ent = entropy2D(Opinion1_Final, Opinion2_Final, 0.1)) %>%
  separate(ent,
           into = c("Far_from_entropy", "One_group_size",
                    "Number_of_equal_groups", "Zero_lenghts", "SD", "DI"),
           sep = "_", convert = T)
b = Sys.time()
b - a

# Saving processed meta indicators:
write_csv(dfy, "Sims02_processed_2D.csv")


