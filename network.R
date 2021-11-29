#### Reading in network (and other) data recorded in NetLogo and saving them back in cleaner form

## Encoding: windows-1250
## Created:  2021-11-19 Francesco
## Edited:   2021-11-29 Francesco


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



