#### Reading first data in! Cleaning it and trying first graphs.

## Encoding: windows-1250
## Created:  2021-10-25 Francesco
## Edited:   2021-11-06 Francesco


## NOTES:
#  Results from NetLogo are lists stored in *.csv file as one string.
#  I have to find a way, how one list transform to several variables.
#  Note: results are stored in file without extender 'resultsV03',
#  it is *.csv file, but I forgot in NetLogo to add extender (it's not automatic).
#
#



# Head --------------------------------------------------------------------

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



# Loading and preparing main data ------------------------------------------------------------

# We load raw data -- we must skip the first 6 lines since there NetLogo writes some meta info, not data itself
data = read_csv("resultsV03.csv", skip = 6) %>%

  # We rename variables -- in NetLogo I use heavily minus ('-') in names and R doesn't like it...
  prejmenuj(c(1, 3:5, 9:10, 20),
            c("rid", "N_agents", "random_links", "neis", "drawn", "speaking", "step")) %>%

  # I find that some cases are impossible for analyses -- they have no components,
  # one instance is natural: after 5000 steps simulation doesn't converge so the components counting algorithm doesn't start
  # but one instance is strange: simulation ends after 99 steps, but it could count only after 100 steps since the record of changes has length 100... it's strange...
  # Ough! Even stranger! Now I see that some simulations end after 199, 200, 201 steps without any component...
  # Hmm... Since I do not store components smaller than 6, it might be that the largest component in simulation
  # is at maximum 5 and that is why its not here. So, for the next experiment I probably should put
  # the threshold lower, omit just components of size 1 or 2, but all others record.
  # It also depends on our questions -- it might be sufficient now replace empty component with
  # component of size 5, or may by simply with NA, so we store info that simulation doesn't produce
  # component larger than 5.
  #
  # And the problem with step 99 -- it probably might be because I put the step-counter after
  # the stopping condition, so, when simulation ends during the step 100 in records there is 99
  # because the simulation stopped just before the counter counted the 100th step.
  # So then it would mean that there was no change during the first 100 steps.
  #
  # TO-DO: Look at simulation program whether it is like that, if so,
  #        then put the step counter at the start or before stopping conditions.
  filter(step < 5000) %>%

  # At least now for data cleaning I need to omit simulations with empty component and positions lists,
  # when these are equal I will omit them.
  filter(components != positions) %>%

  # for purpose of data cleaning I filter out simple cases, in working code this line will be commented
  # filter(opinions > 1) %>%
  # slice_sample(n = 100) %>%

  # We separate list of component sizes into 50 separate variables
  # TO-DO: Check after end of experiment if 50 is enough, now its just 'the first 50 values' and it is sufficient.
  separate(components, sep = " ", into = paste0("size.", 1:50)) %>%

  # We also separate list of components' positions into 50 separate variables
  separate(positions, sep = "\\] \\[", into = paste0("p.", 1:50)) %>%

  # For easy manipulation we reshape data to long form
  # TO-DO: Check the number of variables where we separated positions and components!
  #        Now the code 'thinks' we are separating into 50 variables --> 'p.50'.
  pivot_longer(size.1:p.50) %>%

  # We filter out empty values
  filter(!is.na(value)) %>%

  # Variable 'name' contains now the original var names c.1:p.50,
  # we separate the root of the name from number and the number we store as ID of component.
  separate(name, sep = "\\.", into = c("name", "component_ID")) %>%

  # Now we have values of component size and the opinion position in one var 'value',
  # we reshape date so to have on one line component ID, size and opinion position of such an component.
  # Note: For identification is sufficient 'rid' and 'component_ID', but then R omits all other vars,
  #       so to avoid that and store all needed vars, we use for identification of row also
  #       obsolete vars and specify rows by 'id_cols = c(1:10, 20, 22)'.
  pivot_wider(id_cols = c(1:10, 20, 22), names_from = "name", values_from = "value") %>%

  # Now we have to separate list with opinion position of cluster into separate opinions/coordinates:
  separate(p, sep = " ", into = paste0("opinion.", 1:16)) %>%

  # Finally, we parse numbers as numbers
  mutate(
    across(.cols = opinion.1:opinion.16, ~parse_number(.x)),
    component_ID = parse_number(component_ID),
    size = round(parse_number(size) / N_agents * 100, digits = 2))



# Loading and praparing data on omitted cases -----------------------------

# We load raw data -- we must skip the first 6 lines since there NetLogo writes some meta info, not data itself
reason = read_csv("resultsV03.csv", skip = 6) %>%

  # We rename variables -- in NetLogo I use heavily minus ('-') in names and R doesn't like it...
  prejmenuj(c(1, 3:5, 9:10, 20),
            c("rid", "N_agents", "random_links", "neis", "drawn", "speaking", "step")) %>%

  # Since we know that component size and opinion positions are different types of lists,
  # we might use it: they might be equal only in one case -- when they are empty,
  # so we might identify cases with both empty lists by asking whether they are equal.
  #
  # Now we could construct easy indicator of reason of omitting the case from main file:
  mutate(reason =
           if_else(
             components != positions,
             "Stable",
             if_else(step == 5000, "Non-equilibrium","Fractured")
            ) %>% factor()
         ) %>%

  # We need only some vars so we select them:
  select(1:10, 20, 23)



# Storing data ------------------------------------------------------------

# We store data as two Excel files
write_xlsx(data, "components.xlsx")
write_xlsx(reason, "simulations.xlsx")

