#### Script for processing data from Identity paper experiments
#### Now we follow with sensitivity analysis of HK model with heterogenous parameters


## Encoding: windows-1250
## Created:  2022-11-15 FrK
## Edited:   2022-12-06 FrK

## Notes:
##
##


# Head --------------------------------------------------------------------

# Clearing all
rm(list=ls())

# Packages
library(readr)
library(dplyr)
library(tibble)
library(forcats)
library(ggplot2)



# Loading and processing data ------------------------------------------------------------

# Creating object 'raw' (tibble): Loading....
raw = read_csv("ClassicalHK_heterogenousParameters_RS01-05.csv", skip = 6) %>%
  add_row(read_csv("ClassicalHK_heterogenousParameters_RS06-10.csv", skip = 6))
for (i in seq(11, 56, 5)) {
  raw = raw %>%
    add_row(read_csv(paste0("ClassicalHK_heterogenousParameters_RS", i, "-", i + 4, ".csv"), skip = 6))
}

# Transforming 'raw' to clean 'ts'
ts = raw %>%
  # Selecting and renaming...
  select(HK_distribution = 4, Present_opinion = 5,
         2, Use_identity = 12,
         N = 3, Boundary = 7, 8, Conformity = 10, 11,
         39, diversity = 40, extremness = 41, ESBG = 42) %>%

  # Processing.
  mutate(
    across(.cols = c(1:2, 4), factor),
    even_N = ((N %%2) == 0))


