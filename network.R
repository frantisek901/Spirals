#### Reading in network (and other) data recorded in NetLogo and saving them back in cleaner form

## Encoding: windows-1250
## Created:  2021-11-19 Francesco
## Edited:   2022-01-14 Francesco


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




# Measuring entropy of distances ------------------------------------------

# General entropy measure, i.e., independent of the data dimensions -- easy, we produce in NetLogo dimensionless data:
entropy = function(distances, bin.width = 0.02) {
  # Function computes entropy of distances out of vector of opinion distances

  # Test whwther we put in right data
  if (class(distances) != "numeric") {
    print("Function expects vector of values of class 'numeric'.")
    break
  }

  # Now we use histogram function to separate vector into bins and
  # transform it into counts of these bins.
  cnts = hist(distances, breaks = seq(0, 1, bin.width))$counts  # Transformation

  # Now we transform counts into relative freqencies:
  frqs = cnts[cnts>0] / length(distances) # Filtering out 0s.

  # Now we transform them into product of relative frequency and its logarithm:
  ents = frqs * log(frqs, 2)

  # Now we need the value of maximum entropy with given granularity:
  bin.num = ceiling(1 / bin.width)
  max.ent = log(bin.num, 2)

  # Here in clean environment we prepare all coefficients:
  entropy = abs(round(sum(ents), 3))
  normalized_entropy =  abs(round(sum(ents) / (max.ent) * 100, 1))
  SD = round(sd(distances), 3)
  AD = round(mean(1 - distances) * 100, 3)

  # Returning back four parameters as one string:
  #print(paste(p_ent, p_oneGroup, n_groups, p_lengths, SD, sep = "_"))
  paste(entropy, normalized_entropy, SD, AD, sep = "_")
}


# Function for finding dataset and processing it into entropy measures:
find_and_process = function(catalog, file.order, granularity, path = "D:/!Spirals/") {

  load(paste0(path, catalog$outFile[file.order]))
  entropy(dfi[["Distance"]], granularity)
}

# Entropy of distances makes sense only for "Links01" data:
dl = filter(d, Type == "Links01") %>% mutate(ent.measures = NA_character_)


for (i in 765:nrow(dl)) {
 dl$ent.measures[i] = (find_and_process(dl, i, 0.02))
}


dlp = separate(dl, ent.measures,
               into = c("Entropy", "Normalized.entropy", "SD", "Average.distance"),
               sep = "_", convert = T)

dlp %>% mutate(Average.distance = 2 * Average.distance) %>%
  ggplot(aes(x = Normalized.entropy, y = Average.distance)) +
  geom_point() +
  lab(x = "Normalized entropy", y = "Average distance") +
  theme_minimal()

ggsave("EntropyVsDistance.jpg", height = 6, width = 6)
