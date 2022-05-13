#### Loading and analyzing data from robustness check experiments

## Encoding: windows-1250
## Created:  2022-05-13 Francesco
## Edited:   2022-05-13 Francesco


## NOTES:
#
#



# Head --------------------------------------------------------------------

# Clearing all
rm(list = ls())

# Setting working directory
setwd("./RobustnessCheck/")
getwd()

# Packages
library(dplyr)
library(readr)
library(tidyr)
library(readxl)
library(writexl)
library(ggplot2)
library(knitr)
library(sjmisc)
library(stargazer)


# Loading -----------------------------------------------------------------

df = read_csv("rob010.csv", skip = 6) %>%
  add_row(read_csv("rob020.csv", skip = 6)) %>%
  add_row(read_csv("rob030.csv", skip = 6)) %>%
  add_row(read_csv("rob040.csv", skip = 6)) %>%
  add_row(read_csv("rob050.csv", skip = 6)) %>%
  add_row(read_csv("rob060.csv", skip = 6)) %>%
  add_row(read_csv("rob070.csv", skip = 6)) %>%
  add_row(read_csv("rob080.csv", skip = 6)) %>%
  add_row(read_csv("rob090.csv", skip = 6)) %>%
  add_row(read_csv("rob100.csv", skip = 6)) %>%
  add_row(read_csv("rob110.csv", skip = 6)) %>%
  add_row(read_csv("rob120.csv", skip = 6)) %>%
  add_row(read_csv("rob130.csv", skip = 6)) %>%
  add_row(read_csv("rob140.csv", skip = 6)) %>%
  add_row(read_csv("rob150.csv", skip = 6)) %>%
  add_row(read_csv("robNoID050.csv", skip = 6)) %>%
  add_row(read_csv("robNoID100.csv", skip = 6)) %>%
  add_row(read_csv("robNoID150.csv", skip = 6)) %>%
  select(3:6, 22, 57) %>%
  rename(ESBG = 6) %>%
  mutate(opinions = factor(opinions))


# Regression --------------------------------------------------------------

noInt = lm(ESBG~`use_identity?`+id_threshold+boundary+opinions+`conformity-level`, df)
full = lm(ESBG~`use_identity?`+id_threshold*boundary*opinions+`conformity-level`, df)
summary(full)
stargazer(noInt, full, type = "text")


# Graph -------------------------------------------------------------------
df = df  %>%
  mutate(id_threshold = if_else(`use_identity?`, id_threshold, 0) %>% factor(),
         boundary = factor(boundary))

dfs = df %>% group_by(`use_identity?`, id_threshold, boundary, opinions) %>%
  summarise(ESBG = mean(ESBG)) %>%
  ungroup()

dfs %>%
  ggplot() +
  aes(fill = id_threshold, y = ESBG, x = boundary, color = id_threshold, group = id_threshold) +
  facet_grid(rows = vars(opinions)) +
  # geom_col(position = position_dodge2()) +
  geom_line() +
  geom_point(size = 5, alpha = 0.3) +
  geom_jitter(data = df, size = 1, alpha = 0.05 ) +
  theme_light()
ggsave("graph1.png", width = 8, height = 6)

dfs %>%
  ggplot() +
  aes(fill = id_threshold, y = ESBG, x = boundary, color = id_threshold, group = id_threshold) +
  facet_grid(rows = vars(opinions), cols = vars(id_threshold)) +
  # geom_col(position = position_dodge2()) +
  geom_line() +
  geom_point(size = 5, alpha = 0.3) +
  geom_jitter(data = df, size = 1, alpha = 0.05 ) +
  theme_light()
ggsave("graph2.png", width = 8, height = 6)

