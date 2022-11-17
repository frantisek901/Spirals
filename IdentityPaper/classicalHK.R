#### Script for processing data from Identity paper experiments
#### Now we start with sensitivity analysis of classical HK model


## Encoding: windows-1250
## Created:  2022-11-15 FrK
## Edited:   2022-11-15 FrK

## Notes:
## 1) We have to do very detailed experiments on classical HK and
##    then build on these results the strategy how to proceed more
##    complex experiments with random distribution of boundary & conformity,
##    and also the experiments with identity.
##    We can't do as much experiments (~19k)
##    per one new parameter value and random seed,
##    so now we will carefully explore the simple space and
##    will find interesting areas there. Then we will continue with these areas.
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

# Creating object 'tb' (tibble): Loading....
tb = read_csv("ClassicalHK_moreN.csv", skip = 6) %>%

  # Selecting and renaming...
  select(HK_distribution = 4, Present_opinion = 5,
         2, Use_identity = 12,
         N = 3, Boundary = 7, Conformity = 10,
         39, diversity = 40, extremness = 41, ESBG = 42) %>%

  # Processing.
  mutate(
    across(.cols = 1:4, factor),
    even_N = (N %%2) == 0)


# Creation of aggregated object 'at' (aggregated tibble):
at = tb %>%

  # We group the file according the anly three real variables:
  # group_by(Boundary, Conformity) %>%
  group_by(Boundary, Conformity, even_N) %>%

  # Computing mean and SD for each resulting variable:
  summarise(across(.cols = 6:9, list(mean = mean, sd = sd))) %>%
  ungroup()


# Aggregated tibble for estimation of effect of conformity (tc):
tc = tb %>%

  # We group tibble by N and Boundary, because we are interested in
  # SD in their combinations -- for their combinations there are
  # 10 different values each for one Conformity value,
  # so in fact, we measure SD caused by Conformity, i.e.
  # whether Conformity makes difference or not,
  # we will be able to map out for which combinations of N and Boundary
  # is bigger SD driven by Conformity and for which it is close to 0.
  # We also have to use 'even_N', since we want to store it for later use.
  group_by(Boundary, N, even_N) %>%

  # Now we finally compute SDs:
  summarise(across(.cols = ticks:ESBG, list(sd = sd))) %>% ungroup()


# Heat maps ---------------------------------------------------------------

# Average var value across all N:
heat_map = function(.data = at, .var = "ticks_sd", .y = "Conformity", .x = "Boundary", .title = "") {
  .data %>%
    ggplot() +
    aes(x = eval(str2lang(.x)), y = eval(str2lang(.y)), col = eval(str2lang(.var))) +
    geom_point(shape = 15, size = 13) +
    scale_color_viridis_c() +
    #scale_x_continuous(breaks = seq(0.05, 0.35, 0.05)) +
    scale_y_continuous(breaks = seq(0.1, 1, 0.1)) +
    labs(x = .x, y = .y, col = .var, title = .title) +
    theme_light()
}

# WHole sample
# heat_map() %>% ggsave("map01.png", plot = ., units = "cm", height = 11.5, width = 36)
# heat_map(.var = "ticks_mean") %>% ggsave("map02.png", plot = ., units = "cm", height = 11.5, width = 36)
#
# heat_map(.var = "diversity_sd") %>% ggsave("map03.png", plot = ., units = "cm", height = 11.5, width = 36)
# heat_map(.var = "diversity_mean") %>% ggsave("map04.png", plot = ., units = "cm", height = 11.5, width = 36)
#
# heat_map(.var = "extremness_sd") %>% ggsave("map05.png", plot = ., units = "cm", height = 11.5, width = 36)
# heat_map(.var = "extremness_mean") %>% ggsave("map06.png", plot = ., units = "cm", height = 11.5, width = 36)
#
# heat_map(.var = "ESBG_sd") %>% ggsave("map07.png", plot = ., units = "cm", height = 11.5, width = 36)
# heat_map(.var = "ESBG_mean") %>% ggsave("map08.png", plot = ., units = "cm", height = 11.5, width = 36)

# Odd and even N
# SD
filter(at, !even_N) %>% heat_map(.var = "ticks_sd", .title = "Odd N") %>%
  ggsave("Pics/map21.png", plot = ., units = "cm", height = 11.5, width = 36)

filter(at, even_N) %>% heat_map(.var = "ticks_sd", .title = "Even N") %>%
  ggsave("Pics/map25.png", plot = ., units = "cm", height = 11.5, width = 36)

filter(at, !even_N) %>% heat_map(.var = "diversity_sd", .title = "Odd N") %>%
  ggsave("Pics/map22.png", plot = ., units = "cm", height = 11.5, width = 36)

filter(at, even_N) %>% heat_map(.var = "diversity_sd", .title = "Even N") %>%
  ggsave("Pics/map26.png", plot = ., units = "cm", height = 11.5, width = 36)

filter(at, !even_N) %>% heat_map(.var = "extremness_sd", .title = "Odd N") %>%
  ggsave("Pics/map23.png", plot = ., units = "cm", height = 11.5, width = 36)

filter(at, even_N) %>% heat_map(.var = "extremness_sd", .title = "Even N") %>%
  ggsave("Pics/map27.png", plot = ., units = "cm", height = 11.5, width = 36)

filter(at, !even_N) %>% heat_map(.var = "ESBG_sd", .title = "Odd N") %>%
  ggsave("Pics/map24.png", plot = ., units = "cm", height = 11.5, width = 36)

filter(at, even_N) %>% heat_map(.var = "ESBG_sd", .title = "Even N") %>%
  ggsave("Pics/map28.png", plot = ., units = "cm", height = 11.5, width = 36)

# Mean
filter(at, !even_N) %>% heat_map(.var = "ticks_mean", .title = "Odd N") %>%
  ggsave("Pics/map31.png", plot = ., units = "cm", height = 11.5, width = 36)

filter(at, even_N) %>% heat_map(.var = "ticks_mean", .title = "Even N") %>%
  ggsave("Pics/map35.png", plot = ., units = "cm", height = 11.5, width = 36)

filter(at, !even_N) %>% heat_map(.var = "diversity_mean", .title = "Odd N") %>%
  ggsave("Pics/map32.png", plot = ., units = "cm", height = 11.5, width = 36)

filter(at, even_N) %>% heat_map(.var = "diversity_mean", .title = "Even N") %>%
  ggsave("Pics/map36.png", plot = ., units = "cm", height = 11.5, width = 36)

filter(at, !even_N) %>% heat_map(.var = "extremness_mean", .title = "Odd N") %>%
  ggsave("Pics/map33.png", plot = ., units = "cm", height = 11.5, width = 36)

filter(at, even_N) %>% heat_map(.var = "extremness_mean", .title = "Even N") %>%
  ggsave("Pics/map37.png", plot = ., units = "cm", height = 11.5, width = 36)

filter(at, !even_N) %>% heat_map(.var = "ESBG_mean", .title = "Odd N") %>%
  ggsave("Pics/map34.png", plot = ., units = "cm", height = 11.5, width = 36)

filter(at, even_N) %>% heat_map(.var = "ESBG_mean", .title = "Even N") %>%
  ggsave("Pics/map38.png", plot = ., units = "cm", height = 11.5, width = 36)




# Meaning of conformity ---------------------------------------------------

# Odd and even N: ==> SD
filter(tc, !even_N) %>%
  heat_map(.x = "N", .y = "Boundary", .var = "ticks_sd", .title = "Odd N") %>%
  ggsave("Pics/map51.png", plot = ., units = "cm", height = 11.5, width = 36)

filter(tc, even_N) %>%
  heat_map(.x = "N", .y = "Boundary", .var = "ticks_sd", .title = "Even N") %>%
  ggsave("Pics/map55.png", plot = ., units = "cm", height = 11.5, width = 36)

filter(tc, !even_N) %>%
  heat_map(.x = "N", .y = "Boundary", .var = "diversity_sd", .title = "Odd N") %>%
  ggsave("Pics/map52.png", plot = ., units = "cm", height = 11.5, width = 36)

filter(tc, even_N) %>%
  heat_map(.x = "N", .y = "Boundary", .var = "diversity_sd", .title = "Even N") %>%
  ggsave("Pics/map56.png", plot = ., units = "cm", height = 11.5, width = 36)

filter(tc, !even_N) %>%
  heat_map(.x = "N", .y = "Boundary", .var = "extremness_sd", .title = "Odd N") %>%
  ggsave("Pics/map53.png", plot = ., units = "cm", height = 11.5, width = 36)

filter(tc, even_N) %>%
  heat_map(.x = "N", .y = "Boundary", .var = "extremness_sd", .title = "Even N") %>%
  ggsave("Pics/map57.png", plot = ., units = "cm", height = 11.5, width = 36)

filter(tc, !even_N) %>%
  heat_map(.x = "N", .y = "Boundary", .var = "ESBG_sd", .title = "Odd N") %>%
  ggsave("Pics/map54.png", plot = ., units = "cm", height = 11.5, width = 36)

filter(tc, even_N) %>%
  heat_map(.x = "N", .y = "Boundary", .var = "ESBG_sd", .title = "Even N") %>%
  ggsave("Pics/map58.png", plot = ., units = "cm", height = 11.5, width = 36)



# Slices through data -----------------------------------------------------

## Function for slicing data, N after N...
heat_map_N = function(.data = tb, .N = 101, .var = "ticks", .y = "Conformity", .x = "Boundary", .title = "") {
  .data %>%
    filter(N == .N) %>%
    ggplot() +
    aes(x = eval(str2lang(.x)), y = eval(str2lang(.y)), col = eval(str2lang(.var))) +
    geom_point(shape = 15, size = 13) +
    scale_color_viridis_c() +
    scale_x_continuous(breaks = seq(0.05, 0.35, 0.05)) +
    scale_y_continuous(breaks = seq(0.1, 1, 0.1)) +
    labs(x = .x, y = .y, col = .var, title = .title) +
    theme_light()
}


# # Function, but in form of FOR cycle ----------------
# # for printing slicing graphs for wanted Ns and variables

# Printing
heat_map_slices = function(.populations = unique(tb$N), .titles = c("", "", "", ""),
                           .vars = c("ticks", "diversity", "extremness", "ESBG"),
                           .files = c("Pics/map11_", "Pics/map12_", "Pics/map13_", "Pics/map14_")) {
  if (length(.vars) != length(.files)) stop("Lengths of sliced variables and respective filenames do not match!")
  for (j in 1:length(.vars)) {
    for (i in .populations){
      .i = if_else(i < 100, paste0("0", i), as.character(i))
      heat_map_N(.N = i, .var = .vars[j], .title = .titles[j]) %>%
        ggsave(paste0(.files[j], .i, ".png"), plot = ., units = "cm", height = 11.5, width = 36)
    }
  }
}

# Max N in the 'tb' data:
max(unique(tb$N))

# Systematic use for producing all heatmaps
heat_map_slices(.populations = seq(179, 203, 2), .vars = c("ticks", "diversity", "extremness", "ESBG"),
                .files = paste0("Pics/map", seq(41, 47, 2), "_odd_"),
                .titles = paste("Odd Ns only for", c("'ticks'", "'diversity", "'extremness", "'ESBG'")))
heat_map_slices(.populations = seq(178, 202, 2), .vars = c("ticks", "diversity", "extremness", "ESBG"),
                .files = paste0("Pics/map", seq(42, 48, 2), "_even_"),
                .titles = paste("Even Ns only for", c("'ticks'", "'diversity", "'extremness", "'ESBG'")))



# Regression --------------------------------------------------------------

m1 = lm(ESBG ~ Boundary + Conformity + even_N + log(N), tb)
# summary(m1)

m2 = lm(extremness ~ Boundary + Conformity + even_N + log(N), tb)
# summary(m2)

m3 = lm(diversity ~ Boundary + Conformity + even_N + log(N), tb)
# summary(m3)

m4 = lm(ticks ~ Boundary + Conformity + even_N + log(N), tb)
# summary(m4)

stargazer::stargazer(m1, m2, m3, m4, type = "text")





# Rubish code -------------------------------------------------------------

# The first uses:
heat_map_slices(.populations = seq(17, 155, 2), .vars = c("ESBG"),
                .files = "Pics/map47_odd_", .title = "Odd Ns only for ESBG")
heat_map_slices(.populations = seq(16, 168, 2), .vars = c("ESBG"),
                .files = "Pics/map48_even_", .title = "Even Ns only for ESBG")

heat_map_slices(.populations = seq(17, 167, 2), .vars = c("ticks", "diversity", "extremness"),
                .files = c("Pics/map41_odd_", "Pics/map43_odd_", "Pics/map45_odd_"),
                .titles = paste("Odd Ns only for", c("'ticks'", "'diversity", "'extremness")))
heat_map_slices(.populations = seq(16, 168, 2), .vars = c("ticks", "diversity", "extremness"),
                .files = c("Pics/map42_even_", "Pics/map44_even_", "Pics/map46_even_"),
                .titles = paste("Even Ns only for", c("'ticks'", "'diversity", "'extremness")))

