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

res = read_csv("experiment01part01.csv", skip = 6) %>%
  add_row(read_csv("experiment01part02.csv", skip = 6)) %>%
  add_row(read_csv("experiment01part03.csv", skip = 6))



# Some graphs -------------------------------------------------------------

ggplot(res, aes(x = `n-neis`, y = betweenness_final)) + geom_point()
ggplot(res, aes(x = `n-neis`, y = eigenvector_final)) + geom_point()
ggplot(res, aes(x = `n-neis`, y = clustering_final)) + geom_point()
ggplot(res, aes(x = `n-neis`, y = mean_path_final)) + geom_point()

ggplot(res, aes(y = normalized_polarization_final , x = ESBSG_polarization_final)) + geom_point()
summary(lm(normalized_polarization_final~ESBSG_polarization_final,res))

ggplot(res, aes(y = normalized_polarization_final , x = boundary)) + geom_jitter()
ggplot(res, aes(y = normalized_polarization_final , x = mode)) + geom_boxplot(alpha = 0.2) + geom_jitter()
ggplot(res, aes(y = normalized_polarization_final , x = mode, col = as.factor(boundary))) + geom_boxplot(alpha = 0.2) + geom_jitter()
ggplot(res, aes(y = normalized_polarization_final , x = mode, col = as.factor(id_threshold))) + geom_jitter(alpha = 0.5, size = 2) + theme_minimal()

ggplot(res, aes(x = RS)) + geom_bar()



# Some regressions --------------------------------------------------------

summary(lm(normalized_polarization_final~ESBSG_polarization_final,res))




