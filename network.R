#### Reading in network (and other) data recorded in NetLogo and saving them back in cleaner form

## Encoding: windows-1250
## Created:  2021-11-19 Francesco
## Edited:   2021-12-06 Francesco


## NOTES:
#  Nothing now...
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
library(ggplot2)

# My own functon for renaming in Tidyverse
prejmenuj = function(data, positions, new.names) {
  names(data)[positions] = new.names
  data
}



# Reading in and storing data from experiment -----------------------------------------

# Since data-files are still huge, we have to combine separately data with POPULATIONS 513, 257, 129,
# so we firstly need cycle iterating over ´Population´:
for (k in c(129, 257, 513)) {

  # Since data are long we have to do them seed by seed.
  # So we need secondly the cycle iterating over seeds:
  for (j in 1:10) {

    # We take a vector with filenames of all simulation results from directory 'Sims':
    d = dir(path = "Sims", pattern = paste0("Sims02_", j,"_", k, "_*"))

    # We read the first file and transform it accordingly for the storing in master datafile:
    df = read_csv(paste0("Sims/", d[37])) %>%  # NOTE: 37 is important -- we have to find the first simulation with 2 dimensions of opinion space to more smoothly add cases of data-files from other simulations
      select(-starts_with("Nei")) %>%  # NOTE: without comp. cluster or better code or data storing politics we can't read in also the network data, files are then too huge and personal computer is not able to red, add and store updated dataframes, so for now, we have to information on networks cut out. Now the option is to agree on what we need to prapare from the network data, prepare it and store by each agent or by each simulation.
      mutate(Meta = d[37]) %>%  # We store file name as variable, since the name contains all meta data
      separate(
        Meta, sep = "_",
        into = c(  # We separate meta info into distinctive variables:
          "Sim", "Seed","Population", "Random_links","Close_links", "Dimensions",
          "Updating", "Boundary","Boundary_method", "P_speaking", "Mode")) %>%
      mutate(Mode = str_sub(Mode, end = -5),  # We cut out last 4 characters from new variable mode, since they are ".csv"
             Phase = if_else(Step==0, "Start", "Final"),  # We prepare variable coining final and starting state of agent/simulation
             Step = max(Step)) %>%  # For the next step of data-management we need var ´Step´ constant
      pivot_wider(id_cols = c(ID:Uncertainty, Step:Phase),  # Since a lot of data on final state are redundant -- they are constant for whole simulation -- we create opinion variables for starting and final state and the rest of data is same, so we save a lot of space!
                  names_from = Phase, names_glue = "{.value}_{Phase}",  # We identify cases by all data, except Opinion vars.
                  values_from = starts_with("Opinion")) %>%
      relocate(Sim:Mode, .after = ID) %>% relocate(Step, .after = Mode) %>%
      relocate(starts_with("Opinion"), .after = Uncertainty)  # We logically relocate vars.

    # We save processed data-file.
    # NOTE: since this is the first record, by 'append = F' we erase potential older versin of .csv file.
    write_csv(df, paste0("Sims02_", j,"_", k, ".csv"), append = F)

    # We iterate over whole vector -- we prepare next file for inclusion and
    # then append it into .csv file.
    # NOTE: Code is the very same as the code for the first record,
    #       annotated are only differing parts of code.
    for (i in c(1:36, 38:length(d))) {
    # for (i in c(34:36, 38:39)) {   #
      df2 = read_csv(paste0("Sims/", d[i])) %>%
        select(-starts_with("Nei")) %>%  # NOTE: After commenting this line out the data on network structure stay in the file.
        mutate(Meta = d[i]) %>%
        separate(
          Meta, sep = "_",
          into = c(
            "Sim", "Seed","Population", "Random_links","Close_links", "Dimensions",
            "Updating", "Boundary","Boundary_method", "P_speaking", "Mode")) %>%
        mutate(Mode = str_sub(Mode, end = -5),
               Phase = if_else(Step==0, "Start", "Final"),
               Step = max(Step)) %>%
        pivot_wider(id_cols = c(ID:Uncertainty, Step:Phase),
                    names_from = Phase, names_glue = "{.value}_{Phase}",
                    values_from = starts_with("Opinion")) %>%
        relocate(Sim:Mode, .after = ID) %>% relocate(Step, .after = Mode) %>%
        relocate(starts_with("Opinion"), .after = Uncertainty)

      # Storing record into main .csv file:
      # NOTE: Note 'append = T', it means that we append next record to the existing .csv file.
      write_csv(df2, paste0("Sims02_", j,"_", k, ".csv"), append = T)
    }
  }
}



# Reading files back and joining them: ------------------------------------

# We do it for sure step by step.
# Firstly we make a file for each population size and them join them into one master file.

# Joining files for population size 129:
df129 = read_csv("Sims/Sims02_1_129.csv") %>%
  add_row(read_csv("Sims/Sims02_2_129.csv")) %>%
  add_row(read_csv("Sims/Sims02_3_129.csv")) %>%
  add_row(read_csv("Sims/Sims02_4_129.csv")) %>%
  add_row(read_csv("Sims/Sims02_5_129.csv")) %>%
  add_row(read_csv("Sims/Sims02_6_129.csv")) %>%
  add_row(read_csv("Sims/Sims02_7_129.csv")) %>%
  add_row(read_csv("Sims/Sims02_8_129.csv")) %>%
  add_row(read_csv("Sims/Sims02_9_129.csv")) %>%
  add_row(read_csv("Sims/Sims02_10_129.csv"))
write_csv(df129, "Sims02_pop129.csv")


# Joining files for population size 257:
df257 = read_csv("Sims/Sims02_1_257.csv") %>%
  add_row(read_csv("Sims/Sims02_2_257.csv")) %>%
  add_row(read_csv("Sims/Sims02_3_257.csv")) %>%
  add_row(read_csv("Sims/Sims02_4_257.csv")) %>%
  add_row(read_csv("Sims/Sims02_5_257.csv")) %>%
  add_row(read_csv("Sims/Sims02_6_257.csv")) %>%
  add_row(read_csv("Sims/Sims02_7_257.csv")) %>%
  add_row(read_csv("Sims/Sims02_8_257.csv")) %>%
  add_row(read_csv("Sims/Sims02_9_257.csv")) %>%
  add_row(read_csv("Sims/Sims02_10_257.csv"))
write_csv(df257, "Sims02_pop257.csv")


# Joining files for population size 513:
df513 = read_csv("Sims/Sims02_1_513.csv") %>%
  add_row(read_csv("Sims/Sims02_2_513.csv")) %>%
  add_row(read_csv("Sims/Sims02_3_513.csv")) %>%
  add_row(read_csv("Sims/Sims02_4_513.csv")) %>%
  add_row(read_csv("Sims/Sims02_5_513.csv")) %>%
  add_row(read_csv("Sims/Sims02_6_513.csv")) %>%
  add_row(read_csv("Sims/Sims02_7_513.csv")) %>%
  add_row(read_csv("Sims/Sims02_8_513.csv")) %>%
  add_row(read_csv("Sims/Sims02_9_513.csv")) %>%
  add_row(read_csv("Sims/Sims02_10_513.csv") %>% filter(Dimensions < 3) %>% mutate(Opinion2_Final = as.numeric(Opinion2_Final)))
write_csv(df513, "Sims02_pop513.csv")


# Joining population files into the master file:
dfa = read_csv("Sims02_pop129.csv") %>% select(-starts_with("Nei")) %>%
  add_row(read_csv("Sims02_pop257.csv")) %>%
  add_row(read_csv("Sims02_pop513.csv"))
write_csv(dfa, "Sims02_all.csv")


# Reding the master file back:
dfa = read_csv("Sims02_all.csv")


# Preparing data for the analysis:
dfx = dfa %>% select(-Sim) %>%
  group_by(Seed, Population, Random_links, Close_links, Dimensions, Boundary, Boundary_method, P_speaking,  Mode) %>%
  mutate(Sim_ID = cur_group_id()) %>% relocate(Sim_ID, .before = ID)
write_csv(dfx, "Sims02_all_new.csv")

# Reading back:
dfb = read_csv("Sims02_all_new.csv")
dfx = dfb %>% filter(Dimensions == 2) %>% group_by(Sim_ID) %>%
  summarise(op1 = sd(Opinion1_Final), op2 = sd(Opinion2_Final))
ggplot(dfx, aes(x = op1, y = op2)) + geom_point(alpha = 0.1) + theme_minimal()
ggplot(dfx, aes(x = op1, y = op2)) + geom_density_2d() + theme_minimal()
ggplot(dfx, aes(x = op1)) + geom_density() + theme_minimal()
ggplot(dfx, aes(x = op2)) + geom_density() + theme_minimal()

dfy = dfb %>% filter(Sim_ID == 39)

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

  # Returning back four parameters as one string:
  #print(paste(p_ent, p_oneGroup, n_groups, p_lengths, SD, sep = "_"))
  paste(p_ent, p_oneGroup, n_groups, p_lengths, SD, sep = "_")
}

dfp = dfb %>% filter(Dimensions == 1, Sim_ID <= 4) %>% group_by(Sim_ID)
de = entropy1D(dfp$Opinion1_Final, 0.1)
de



a = Sys.time()
dfy = dfb %>%
  # filter(Sim_ID <= 10) %>%
  filter(Dimensions == 1) %>%
  group_by(Sim_ID, Seed, Population, Random_links, Close_links, Dimensions, Boundary, Boundary_method, P_speaking,  Mode) %>%
  summarise(ent = entropy1D(Opinion1_Final, 0.1)) %>%
  separate(ent,
           into = c("Far_from_entropy", "One_group_size",
                    "Number_of_equal_groups", "Zero_lenghts", "SD"),
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

  # Returning back four parameters as one string:
  #print(paste(p_ent, p_oneGroup, n_groups, p_lengths, SD, sep = "_"))
  paste(p_ent, p_oneGroup, n_groups, p_lengths, SD, sep = "_")
}

dfp = dfb %>% filter(Dimensions == 2, Sim_ID < 40) %>% group_by(Sim_ID)
de = entropy2D(dfp$Opinion1_Final, dfp$Opinion2_Final, 0.1)
de

a = Sys.time()
dfy = dfb %>%
  # filter(Sim_ID <= 72) %>%
  filter(Dimensions == 2) %>%
  group_by(Sim_ID, Seed, Population, Random_links, Close_links, Dimensions, Boundary, Boundary_method, P_speaking,  Mode) %>%
  summarise(ent = entropy2D(Opinion1_Final, Opinion2_Final, 0.1)) %>%
  separate(ent,
           into = c("Far_from_entropy", "One_group_size",
                    "Number_of_equal_groups", "Zero_lenghts", "SD"),
           sep = "_", convert = T)
b = Sys.time()
b - a

# Saving processed meta indicators:
write_csv(dfy, "Sims02_processed_2D.csv")


# Joining files -----------------------------------------------------------

# Firstly, we need to join files/add rows from different entropy dimensions,
# later we add/join also files fro HK-benchmarking experiment
df = dfy %>% ungroup() %>% add_row(read_csv("Sims02_processed_2D.csv")) %>%
  filter(P_speaking > 0.4)



# Entropy graphs ----------------------------------------------------------

ggplot(df, aes(x = Far_from_entropy, y = Zero_lenghts)) +
  facet_grid(cols = vars(Dimensions), rows = vars(Boundary_method)) +
  geom_point(alpha = 0.1) +
  theme_minimal()

ggplot(df, aes(x = Far_from_entropy, y = One_group_size)) +
  facet_grid(cols = vars(Dimensions), rows = vars(Mode)) +
  geom_point(alpha = 0.1) +
  scale_x_log10() +
  theme_minimal()

ggplot(df, aes(x = Far_from_entropy, y = Number_of_equal_groups)) +
  facet_grid(cols = vars(Dimensions), rows = vars(Boundary)) +
  geom_point(alpha = 0.1) +
  scale_x_log10() +
  scale_y_log10(breaks = c(seq(1, 7, 2), 10, 20, 30)) +
  theme_minimal()

ggplot(df, aes(x = Far_from_entropy, y = SD)) +
  facet_grid(cols = vars(Dimensions), rows = vars(Population)) +
  geom_point(alpha = 0.1) +
  # scale_x_log10() +
  # scale_y_log10() +
  theme_minimal()


# Graph sketches ----------------------------------------------------------

ggplot(filter(dfa, Boundary_method != "constant" ),
       aes(x = Uncertainty, y = abs(Opinion1_Final - Opinion1_Start))) +
  facet_wrap(vars(Mode)) +
  geom_point(alpha = 0.01) +
  theme_minimal()
ggsave("try0.jpg")


dfa %>% filter(Boundary_method != "constant" & Dimensions == 2 & Population == 129) %>%
  ggplot(aes(x = abs(Opinion1_Final - Opinion1_Start),
             y = abs(Opinion2_Final - Opinion2_Start))) +
  facet_grid(cols = vars(Close_links), rows = vars(Boundary)) +
  geom_point(alpha = 0.01) +
  scale_x_log10() +
  scale_y_log10() +
  theme_minimal()
ggsave("try1.jpg", width = 8, height = 8)


dfa %>% filter(Boundary_method != "constant" & Dimensions == 2) %>%
  mutate(Close_links = round(Close_links * 2) / (Population - 1) * 100, 1) %>%
  ggplot(aes(x = abs(Opinion1_Final - Opinion1_Start),
             y = abs(Opinion2_Final - Opinion2_Start))) +
  facet_grid(cols = vars(Close_links), rows = vars(Boundary)) +
  geom_point(alpha = 0.01) +
  theme_minimal()
ggsave("try2.jpg", width = 14, height = 6)


dfa %>% filter(Boundary_method != "constant" & Dimensions == 2) %>%
  ggplot(aes(x = abs(Opinion1_Final - Opinion1_Start),
             y = abs(Opinion2_Final - Opinion2_Start),
             col = Mode)) +
  facet_grid(cols = vars(Close_links, Population), rows = vars(Boundary, Mode)) +
  geom_point(alpha = 0.01) +
  guides(col = "none") +
  theme_minimal()
ggsave("try3.jpg", width = 18, height = 12)



# ### This is the end my only friend, the end! ### ------------------------



# Old code ----------------------------------------------------------------

entropy2D_show = function(v1, v2, bin.width = 0.1) {
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
  #max.ent = log(1 / length(ents))


  #print(cnts)
  #print(frqs)
  print(ggplot(tibble(v1, v2), aes(x = v1, y = v2)) +
          geom_jitter(alpha = 0.2, height = bin.width / 10, width = bin.width / 10) +
          scale_x_continuous(breaks = seq(-1, 1, 0.2)) +
          scale_y_continuous(breaks = seq(-1, 1, 0.2)) +
          theme_minimal())
  #print(bin.num)
  #print(max.ent)
  # print(paste("Effective bins:", length(cnts), "out of", bin.num))
  print(paste("Average percentage of lengths close to 0 from theoretical maximum:", round(frqs[1] * 100, 1)))
  print(paste("Percentage of agents enough to cover in one group all lenghts close to 0:", round(sqrt(cnts[1] * 2) / l * 100, 1)))
  print(paste("Minimum number of equal groups capable to cover all lenghts close to 0:", round(l * l / (l + (2 * cnts[1])), 2)))
  print(paste("Number of lengths:", ld, "; Number of agents:", l, "; Number of lenghts close to 0:", cnts[1]))
  print(paste("Percentage of uniform public sphere:", round(100 - (sum(ents) / (max.ent) * 100), 1)))
  # print(max.ent2)
  # print(sum(ents) / (max.ent2) * 100)
  # print(tibble(cnts, frqs, ents) %>% arrange(desc(cnts)))
  d
}

# Sim_ID == 12140 well illustrates the measures
dfy = dfb %>% filter(Sim_ID == 12140) #%>% group_by(Sim_ID)
de = entropy2D(dfy$Opinion1_Final, dfy$Opinion2_Final, 0.1)
hist(de, breaks = seq(0, 3, 0.1))



