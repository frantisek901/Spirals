#### Script analyzing results from the first real experiment with Spirals' project model

## Encoding: windows-1250
## Created:  2022-03-03 Francesco
## Edited:   2022-03-07 Francesco


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

res = read_csv("experiment01part01.csv", skip = 6) %>%
  add_row(read_csv("experiment01part02.csv", skip = 6)) %>%
  add_row(read_csv("experiment01part03.csv", skip = 6)) %>%
  add_row(read_csv("experiment01part04.csv", skip = 6) %>% mutate(mean_path_final = as.numeric(mean_path_final))) %>%
  add_row(read_csv("experiment01part05.csv", skip = 6) %>% mutate(mean_path_final = as.numeric(mean_path_final))) %>%
  add_row(read_csv("experiment01part06.csv", skip = 6) %>% mutate(mean_path_final = as.numeric(mean_path_final))) %>%
  add_row(read_csv("experiment01part11.csv", skip = 6) %>% mutate(mean_path_final = as.numeric(mean_path_final))) %>%
  add_row(read_csv("experiment01part13.csv", skip = 6) %>% mutate(mean_path_final = as.numeric(mean_path_final))) %>%
  add_row(read_csv("experiment01part14.csv", skip = 6) %>% mutate(mean_path_final = as.numeric(mean_path_final))) %>%
  select(-c(3:4, 7, 9:12, 14:15, 17:18, 20:21, 23:27, 29, 32, 34:35, 38, 40:43, 45:47)) %>%
  mutate(
    iqr_op1_start = upper_op1_start - lower_op1_start,
    iqr_op2_start = upper_op2_start - lower_op2_start,
    iqr_op1_final = upper_op1_final - lower_op1_final,
    iqr_op2_final = upper_op2_final - lower_op2_final
  ) %>%
  # select(-c(31:34, 47:50)) %>%
  relocate(any_of(c("iqr_op1_start", "iqr_op2_start")), .after = median_op2_start)



long = read_csv("experiment01part03LONG.csv", skip = 6) %>%
  mutate(mean_path_final = as.numeric(mean_path_final)) %>%
  add_row(read_csv("experiment01part01LONG01.csv", skip = 6)) %>%
  add_row(read_csv("experiment01part01LONG02.csv", skip = 6)) %>%
  add_row(read_csv("experiment01part01LONG02b.csv", skip = 6)) %>%
  add_row(read_csv("experiment01part01LONG04.csv", skip = 6)) %>%
  add_row(read_csv("experiment01part03LONG.csv", skip = 6)) %>%
  select(-c(3:4, 7, 9:12, 14:15, 17:18, 20:21, 23:27, 29, 32, 34:35, 38, 40:43, 45:47)) %>%
  mutate(
    iqr_op1_start = upper_op1_start - lower_op1_start,
    iqr_op2_start = upper_op2_start - lower_op2_start,
    iqr_op1_final = upper_op1_final - lower_op1_final,
    iqr_op2_final = upper_op2_final - lower_op2_final
  ) %>%
  # select(-c(31:34, 47:50)) %>%
  relocate(any_of(c("iqr_op1_start", "iqr_op2_start")), .after = median_op2_start)



test = read_csv("speedTestingData.csv", skip = 6)

wrng = read_csv("wronglySpecifiedExperiment.csv", skip = 6)



# Saving data -------------------------------------------------------------

save(res, file = "shortData.RData")

save(long, file = "longData.RData")



# Some graphs -------------------------------------------------------------

ggplot(res, aes(x = iqr_op1_final, y = iqr_op2_final, col = `use_identity?`)) +
  geom_point(alpha = 0.15) +
  theme_minimal()

ggplot(long, aes(x = iqr_op1_final, y = iqr_op2_final, col = `use_identity?`)) +
  geom_point(alpha = 0.15) +
  theme_minimal()

ggplot(res, aes(x = `n-neis`, y = betweenness_final)) + geom_point()
ggplot(res, aes(x = `n-neis`, y = eigenvector_final)) + geom_point()
ggplot(res, aes(x = `n-neis`, y = clustering_final)) + geom_point()
ggplot(res, aes(x = `n-neis`, y = mean_path_final)) + geom_point()

ggplot(res, aes(y = normalized_polarization_final , x = ESBSG_polarization_final)) + geom_point()
summary(lm(normalized_polarization_final~ESBSG_polarization_final,res))

ggplot(res, aes(y = normalized_polarization_final , x = boundary)) + geom_jitter()
ggplot(res, aes(y = normalized_polarization_final , x = mode)) + geom_boxplot(alpha = 0.2) + geom_jitter()
ggplot(res, aes(y = normalized_polarization_final , x = mode, col = as.factor(boundary))) + geom_boxplot(alpha = 0.2) + geom_jitter()
ggplot(res, aes(y = normalized_polarization_final , x = mode, col = as.factor(id_threshold))) +
  facet_wrap(vars(`use_identity?`)) +
  geom_jitter(alpha = 0.05, size = 2) +
  theme_minimal()

ggplot(res, aes(x = RS)) + geom_bar()



ggplot(res, aes(y = normalized_polarization_final , x = mode, col = as.factor(id_threshold))) +
  facet_grid(rows = vars(`use_identity?`), cols = vars(boundary, `tolerance-level`)) +
  geom_jitter(alpha = 0.05, size = 2) +
  theme_minimal()


res %>% group_by(mode, id_threshold, `use_identity?`, boundary, `tolerance-level`) %>%
  summarise(npf = mean(normalized_polarization_final)) %>%
  ggplot(aes(y = npf , x = mode, fill = as.factor(id_threshold))) +
  facet_grid(rows = vars(`use_identity?`), cols = vars(boundary, `tolerance-level`)) +
  geom_col(alpha = 0.5) +
  theme_minimal()



# Some regressions --------------------------------------------------------

# 365 ticks
summary(lm(normalized_polarization_final~boundary+mode+id_threshold+`use_identity?`+`tolerance-level`+`p-speaking-level`+`conformity-level`+`p-random`+`n-neis`,res))
summary(lm(ESBSG_polarization_final~boundary+mode+id_threshold+`use_identity?`+`tolerance-level`+`p-speaking-level`+`conformity-level`+`p-random`+`n-neis`,res))

# 3650 ticks
summary(lm(normalized_polarization_final~boundary+mode+id_threshold+`use_identity?`+`tolerance-level`+`p-speaking-level`+`conformity-level`+`p-random`+`n-neis`, long))
summary(lm(ESBSG_polarization_final~boundary+mode+id_threshold+`use_identity?`+`tolerance-level`+`p-speaking-level`+`conformity-level`+`p-random`+`n-neis`, long))



