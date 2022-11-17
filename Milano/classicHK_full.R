#### Script for analyzing full-factor search of classical HK model

## Encoding: windows-1250
## Created:  2022-08-14 FrK
## Edited:   2022-08-14 FrK

## Notes:
#



# Head --------------------------------------------------------------------

# Clear all
rm(list=ls())


# Packages
library(ggplot2)
library(tidyr)
library(readr)
library(dplyr)
library(stargazer)



# Loading and preparing data ----------------------------------------------

# Loading...
df = read_csv("bSearch/RS001-200_HKClassicFullFactorial.csv", skip = 6) %>% 
  
  # Selecting...
  select(openness = boundary, conformity = `conformity-level`, ESBG = ESBG_polarisation)



# Aggregatted graph -------------------------------------------------------

# Preparing data
dfa = df %>% 
  
  # Grouping by main drivers of polarization
  group_by(openness, conformity) %>% 
  
  # Summarizing -- mean and median ESBG
  summarise(mean = mean(ESBG), median = median(ESBG))
  

# Drawing graph:
dfa %>% ggplot() +
  aes(y = conformity, x = openness, color = mean, fill = mean, label = round(mean, 2)) +
  geom_point(shape = 22, size = 11) +
  geom_text(color = "white", size = 3) +
  theme_minimal()



# Regression --------------------------------------------------------------

# Formulating regression model:
m = lm(mean ~ conformity + openness, data = dfa)

# Printing nice output:
stargazer(m, type = "text")
# OK, it's evident, that openness is primary factor, conformity is the secondary.
