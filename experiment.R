#### Script analyzing results from the first real experiment with Spirals' project model

## Encoding: windows-1250
## Created:  2022-03-03 Francesco
## Edited:   2022-03-03 Francesco


## NOTES:
#
#



# Head --------------------------------------------------------------------

# Clearing all
rm(list = ls())

# Setting working directory
setwd("./Experiment/")
getwd()

# Packages
library(dplyr)
library(readr)
library(tidyr)
library(readxl)
library(writexl)
library(ggplot2)
library(knitr)


# My own functon for renaming in Tidyverse
prejmenuj = function(data, positions, new.names) {
  names(data)[positions] = new.names
  data
}



# Loading and cleaning data -----------------------------------------------

test = read_csv("speedTestingData.csv", skip = 6)

wrng = read_csv("wronglySpecifiedExperiment.csv", skip = 6)








