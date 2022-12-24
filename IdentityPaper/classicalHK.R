#### Script for processing data from Identity paper experiments
#### Now we start with sensitivity analysis of classical HK model


## Encoding: windows-1250
## Created:  2022-11-15 FrK
## Edited:   2022-12-06 FrK

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
    across(.cols = c(1:2, 4), factor),
    even_N = ((N %%2) == 0))


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
  summarise(across(.cols = ticks:ESBG, list(sd = sd))) %>% ungroup() %>%
  mutate(
    r = (diversity_sd - min(diversity_sd)) / (max(diversity_sd) - min(diversity_sd)),
    g = (extremness_sd - min(extremness_sd)) / (max(extremness_sd) - min(extremness_sd)),
    b = (ESBG_sd - min(ESBG_sd)) / (max(ESBG_sd) - min(ESBG_sd)),
    color = (r + g + b) / 3)


# Creating 'ted' (tibble with experiments deviating from classical HK):
td1 = read_csv("ClassicalHK-PresentOpinionRS01-10.csv", skip = 6) %>%
  add_row(read_csv("ClassicalHK-PresentOpinionRS11-20.csv", skip = 6)) %>%
  add_row(read_csv("ClassicalHK-PresentOpinionRS21-30.csv", skip = 6)) %>%
  add_row(read_csv("ClassicalHK-PresentOpinionRS31-40.csv", skip = 6)) %>%
  add_row(read_csv("ClassicalHK-PresentOpinionRS41-50.csv", skip = 6)) %>%
  add_row(read_csv("ClassicalHK-PresentOpinionRS51-60.csv", skip = 6)) %>%
  mutate(file = "present opinion")

td2 = read_csv("ClassicalHK-RandomPositionAtStartRS01-10.csv", skip = 6) %>%
  add_row(read_csv("ClassicalHK-RandomPositionAtStartRS11-20.csv", skip = 6)) %>%
  add_row(read_csv("ClassicalHK-RandomPositionAtStartRS21-30.csv", skip = 6)) %>%
  add_row(read_csv("ClassicalHK-RandomPositionAtStartRS31-40.csv", skip = 6)) %>%
  add_row(read_csv("ClassicalHK-RandomPositionAtStartRS41-50.csv", skip = 6)) %>%
  add_row(read_csv("ClassicalHK-RandomPositionAtStartRS51-60.csv", skip = 6)) %>%
  mutate(file = "random position")

td3 = read_csv("ClassicalHK-RandomPositionAtStartPresentOpinionRS01-10.csv", skip = 6) %>%
  add_row(read_csv("ClassicalHK-RandomPositionAtStartPresentOpinionRS11-20.csv", skip = 6)) %>%
  add_row(read_csv("ClassicalHK-RandomPositionAtStartPresentOpinionRS21-30.csv", skip = 6)) %>%
  add_row(read_csv("ClassicalHK-RandomPositionAtStartPresentOpinionRS31-40.csv", skip = 6)) %>%
  add_row(read_csv("ClassicalHK-RandomPositionAtStartPresentOpinionRS41-50.csv", skip = 6)) %>%
  add_row(read_csv("ClassicalHK-RandomPositionAtStartPresentOpinionRS51-60.csv", skip = 6)) %>%
  mutate(file = "random + present")


# Main individual file cleaned:
ted = td1 %>% add_row(td2) %>% add_row(td3) %>%

  # Selecting and renaming...
  select(43, HK_distribution = 4, Present_opinion = 5,
         2, Use_identity = 12,
         N = 3, Boundary = 7, Conformity = 10,
         39, diversity = 40, extremness = 41, ESBG = 42) %>%

  # Processing.
  mutate(
    across(.cols = c(1:3, 5), factor),
    even_N = ((N %%2) == 0))

# Test Whether the factors and files match:
count(ted, HK_distribution, Present_opinion, file)


# Processing individual partial files back from TED:
td1 = ted %>%
  filter(file == "present opinion") %>%
  group_by(N, Boundary, Conformity, even_N) %>%

  # Computing mean and SD for each resulting variable:
  summarise(across(.cols = 6:9, list(mean = mean, sd = sd))) %>%
  ungroup()

td2 = ted %>%
  filter(file == "random position") %>%
  group_by(N, Boundary, Conformity, even_N) %>%

  # Computing mean and SD for each resulting variable:
  summarise(across(.cols = 6:9, list(mean = mean, sd = sd))) %>%
  ungroup()

td3 = ted %>%
  filter(file == "random + present") %>%
  group_by(N, Boundary, Conformity, even_N) %>%

  # Computing mean and SD for each resulting variable:
  summarise(across(.cols = 6:9, list(mean = mean, sd = sd))) %>%
  ungroup()

# And also complete file processed with all conditions (tda):
tda = ted %>%
  add_row(filter(tb, N %in% unique(ted$N)) %>% mutate(file = "classic")) %>%

  group_by(file, N, Boundary, Conformity, even_N) %>%

  # Computing mean and SD for each resulting variable:
  summarise(across(.cols = ticks:ESBG, list(mean = mean, sd = sd))) %>%
  ungroup() %>%

  # Ordering 'file'
  mutate(file = factor(file, levels = c("classic", "present opinion", "random position", "random + present")))


# And hopefully lastly, agregation over different files (only 4 values -- ada):
ada = tda %>%
  # We group against N and file:
  group_by(Boundary, Conformity, even_N) %>%

  # We compute SD across the file and N:
  summarise(across(.cols = c(ticks_mean, diversity_mean, extremness_mean, ESBG_mean), list(sd = sd)))

# Same file, but containing also info on N:
adaN = tda %>%
  # We again group, but now only against file:
  group_by(N, Boundary, Conformity, even_N) %>%

  # Summarise according N, Boundary, Conformity, and Even_N:
  summarise(across(.cols = c(ticks_mean, diversity_mean, extremness_mean, ESBG_mean), list(sd = sd))) %>%
  ungroup()





# Heat maps ---------------------------------------------------------------

# Average var value across all N:
heat_map = function(.data = at, .var = "ticks_sd", .y = "Conformity", .x = "Boundary", .title = "") {
  .data %>%
    ggplot() +
    aes(x = eval(str2lang(.x)), y = eval(str2lang(.y)), col = eval(str2lang(.var))) +
    geom_point(shape = 15, size = 13) +
    scale_color_viridis_c() +
    scale_x_continuous(breaks = seq(0.05, 0.35, 0.05)) +
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


# SD of means
filter(ada, !even_N) %>% heat_map(.var = "ticks_mean_sd", .title = "Odd N") %>%
  ggsave("Pics/map61.png", plot = ., units = "cm", height = 11.5, width = 36)

filter(ada, even_N) %>% heat_map(.var = "ticks_mean_sd", .title = "Even N") %>%
  ggsave("Pics/map65.png", plot = ., units = "cm", height = 11.5, width = 36)

filter(ada, !even_N) %>% heat_map(.var = "diversity_mean_sd", .title = "Odd N") %>%
  ggsave("Pics/map62.png", plot = ., units = "cm", height = 11.5, width = 36)

filter(ada, even_N) %>% heat_map(.var = "diversity_mean_sd", .title = "Even N") %>%
  ggsave("Pics/map66.png", plot = ., units = "cm", height = 11.5, width = 36)

filter(ada, !even_N) %>% heat_map(.var = "extremness_mean_sd", .title = "Odd N") %>%
  ggsave("Pics/map63.png", plot = ., units = "cm", height = 11.5, width = 36)

filter(ada, even_N) %>% heat_map(.var = "extremness_mean_sd", .title = "Even N") %>%
  ggsave("Pics/map67.png", plot = ., units = "cm", height = 11.5, width = 36)

filter(ada, !even_N) %>% heat_map(.var = "ESBG_mean_sd", .title = "Odd N") %>%
  ggsave("Pics/map64.png", plot = ., units = "cm", height = 11.5, width = 36)

filter(ada, even_N) %>% heat_map(.var = "ESBG_mean_sd", .title = "Even N") %>%
  ggsave("Pics/map68.png", plot = ., units = "cm", height = 11.5, width = 36)





# Meaning of conformity ---------------------------------------------------

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



# Three measures together as color ----------------------------------------

filter(tc, !even_N) %>%
  heat_map(.x = "N", .y = "Boundary", .var = "color", .title = "Odd N") %>%
  ggsave("Pics/map59o.png", plot = ., units = "cm", height = 11.5, width = 36)

filter(tc, even_N) %>%
  heat_map(.x = "N", .y = "Boundary", .var = "color", .title = "Even N") %>%
  ggsave("Pics/map59e.png", plot = ., units = "cm", height = 11.5, width = 36)



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
heat_map_slices = function(.dat = tb, .populations = unique(tb$N), .titles = c("", "", "", ""),
                           .vars = c("ticks", "diversity", "extremness", "ESBG"),
                           .files = c("Pics/map11_", "Pics/map12_", "Pics/map13_", "Pics/map14_")) {
  if (length(.vars) != length(.files)) stop("Lengths of sliced variables and respective filenames do not match!")
  for (j in 1:length(.vars)) {
    for (i in .populations){
      .i = if_else(i < 100, paste0("0", i), as.character(i))
      heat_map_N(.data = .dat, .N = i, .var = .vars[j], .title = .titles[j]) %>%
        ggsave(paste0(.files[j], .i, ".png"), plot = ., units = "cm", height = 11.5, width = 36)
    }
  }
}

# Max N in the 'tb' data:
max(unique(tb$N))

# Systematic use for producing all heatmaps
heat_map_slices(.populations = seq(17, 257, 2), .vars = c("ticks", "diversity", "extremness", "ESBG"),
                .files = paste0("Pics/map", seq(41, 47, 2), "_odd_"),
                .titles = paste("Odd Ns only for", c("'ticks'", "'diversity", "'extremness", "'ESBG'")))

heat_map_slices(.populations = seq(16, 256, 2), .vars = c("ticks", "diversity", "extremness", "ESBG"),
                .files = paste0("Pics/map", seq(42, 48, 2), "_even_"),
                .titles = paste("Even Ns only for", c("'ticks'", "'diversity", "'extremness", "'ESBG'")))


# Systematic use for producing all heatmaps
heat_map_slices(.dat = adaN, .populations = c(21, 27, 33,  51, 65, 101, 129, 201, 257),
                .vars = c("ticks_mean_sd", "diversity_mean_sd", "extremness_mean_sd", "ESBG_mean_sd"),
                .files = paste0("Pics/map", seq(71, 77, 2), "_odd_"),
                .titles = paste("Odd Ns only for", c("'ticks'", "'diversity", "'extremness", "'ESBG'")))

heat_map_slices(.dat = adaN, .populations = c(20,  26, 32, 50, 64, 100, 128, 200, 256),
                .vars = c("ticks_mean_sd", "diversity_mean_sd", "extremness_mean_sd", "ESBG_mean_sd"),
                .files = paste0("Pics/map", seq(72, 78, 2), "_even_"),
                .titles = paste("Even Ns only for", c("'ticks'", "'diversity", "'extremness", "'ESBG'")))


# Systematic use for producing all heatmaps, also by file
heat_map_slices(.dat = filter(tda, file == "present opinion"),
                .populations = c(21, 27, 33,  51, 65, 101, 129, 201, 257),
                .vars = c("ticks_mean", "diversity_mean", "extremness_mean", "ESBG_mean"),
                .files = paste0("Pics/map", seq(81, 87, 2), "_odd_"),
                .titles = paste("Present opinion sample: Odd Ns only for", c("'ticks'", "'diversity", "'extremness", "'ESBG'")))

heat_map_slices(.dat = filter(tda, file == "present opinion"),
                .populations = c(20,  26, 32, 50, 64, 100, 128, 200, 256),
                .vars = c("ticks_mean", "diversity_mean", "extremness_mean", "ESBG_mean"),
                .files = paste0("Pics/map", seq(82, 88, 2), "_even_"),
                .titles = paste("Present opinion sample: Even Ns only for", c("'ticks'", "'diversity", "'extremness", "'ESBG'")))

heat_map_slices(.dat = filter(tda, file == "present opinion"),
                .populations = c(21, 27, 33,  51, 65, 101, 129, 201, 257),
                .vars = c("ticks_sd", "diversity_sd", "extremness_sd", "ESBG_sd"),
                .files = paste0("Pics/map", seq(91, 97, 2), "_odd_"),
                .titles = paste("Present opinion sample: Odd Ns only for", c("'ticks'", "'diversity", "'extremness", "'ESBG'")))

heat_map_slices(.dat = filter(tda, file == "present opinion"),
                .populations = c(20,  26, 32, 50, 64, 100, 128, 200, 256),
                .vars = c("ticks_sd", "diversity_sd", "extremness_sd", "ESBG_sd"),
                .files = paste0("Pics/map", seq(92, 98, 2), "_even_"),
                .titles = paste("Present opinion sample: Even Ns only for", c("'ticks'", "'diversity", "'extremness", "'ESBG'")))



# Files/conditions comparisons --------------------------------------------

heat_map_facet = function(.data = tda, .var = "ESBG_mean", .y = "Conformity", .x = "Boundary",
                          .facet = "file", .title = "") {
  .data %>%
    ggplot() +
    aes(x = eval(str2lang(.x)), y = eval(str2lang(.y)), col = eval(str2lang(.var))) +
    facet_wrap(vars(eval(str2lang(.facet)) %>% fct_rev()), ncol = 1) +
    geom_point(shape = 15, size = 13) +
    scale_color_viridis_c() +
    #scale_x_continuous(breaks = seq(0.05, 0.35, 0.05)) +
    scale_y_continuous(breaks = seq(0.1, 1, 0.1)) +
    scale_x_continuous(breaks = seq(0.05, 0.35, 0.02)) +
    labs(x = .x, y = .y, col = .var, title = .title) +
    theme_light() +
    theme(legend.position = "bottom")
}

.height = 46
.width = 34.5

filter(tda, !even_N) %>% heat_map_facet(.var = "ticks_mean", .title = "Odd N")  %>%
  ggsave("Pics/mapa1o.png", plot = ., units = "cm", height = .height, width = .width)

filter(tda, even_N) %>% heat_map_facet(.var = "ticks_mean", .title = "Even N")  %>%
  ggsave("Pics/mapa2e.png", plot = ., units = "cm", height = .height, width = .width)

filter(tda, !even_N) %>% heat_map_facet(.var = "diversity_mean", .title = "Odd N")  %>%
  ggsave("Pics/mapa3o.png", plot = ., units = "cm", height = .height, width = .width)

filter(tda, even_N) %>% heat_map_facet(.var = "diversity_mean", .title = "Even N")  %>%
  ggsave("Pics/mapa4e.png", plot = ., units = "cm", height = .height, width = .width)

filter(tda, !even_N) %>% heat_map_facet(.var = "extremness_mean", .title = "Odd N")  %>%
  ggsave("Pics/mapa5o.png", plot = ., units = "cm", height = .height, width = .width)

filter(tda, even_N) %>% heat_map_facet(.var = "extremness_mean", .title = "Even N")  %>%
  ggsave("Pics/mapa6e.png", plot = ., units = "cm", height = .height, width = .width)

filter(tda, !even_N) %>% heat_map_facet(.title = "Odd N")  %>%
  ggsave("Pics/mapa7o.png", plot = ., units = "cm", height = .height, width = .width)

filter(tda, even_N) %>% heat_map_facet(.title = "Even N")  %>%
  ggsave("Pics/mapa8e.png", plot = ., units = "cm", height = .height, width = .width)


# Regression --------------------------------------------------------------

# Regression just on classical conditions of HK
m1 = lm(ESBG ~ Boundary + Conformity + even_N + log(N), tb)
# summary(m1)

m2 = lm(extremness ~ Boundary + Conformity + even_N + log(N), tb)
# summary(m2)

m3 = lm(diversity ~ Boundary + Conformity + even_N + log(N), tb)
# summary(m3)

m4 = lm(ticks ~ Boundary + Conformity + even_N + log(N), tb)
# summary(m4)

stargazer::stargazer(m1, m2, m3, m4, type = "text")



# Regression compariong effect of different non-classy conditions on MEAN:
m11 = lm(ESBG_mean ~  file + Boundary + Conformity + even_N + log(N), tda)
# summary(m11)

m12 = lm(extremness_mean ~  file + Boundary + Conformity + even_N + log(N), tda)
# summary(m12)

m13 = lm(diversity_mean ~  file + Boundary + Conformity + even_N + log(N), tda)
# summary(m13)

m14 = lm(ticks_mean ~  file + Boundary + Conformity + even_N + log(N), tda)
# summary(m14)

stargazer::stargazer(m11, m12, m13, m14, type = "text")


# Regression compariong effect of different non-classy conditions on SD:
m21 = lm(ESBG_sd ~  file + Boundary + Conformity + even_N + log(N), tda)
# summary(m11)

m22 = lm(extremness_sd ~  file + Boundary + Conformity + even_N + log(N), tda)
# summary(m12)

m23 = lm(diversity_sd ~  file + Boundary + Conformity + even_N + log(N), tda)
# summary(m13)

m24 = lm(ticks_sd ~  file + Boundary + Conformity + even_N + log(N), tda)
# summary(m14)

stargazer::stargazer(m21, m22, m23, m24, type = "text")


# Regression: another model building:
m31 = lm(ESBG_mean ~  log(N) + file, filter(tda, Boundary > 0.15 & Boundary <= 0.25))
m32 = lm(ESBG_mean ~  log(N) + file + even_N, filter(tda, Boundary > 0.15 & Boundary <= 0.25))
m33 = lm(ESBG_mean ~  log(N) + file * even_N, filter(tda, Boundary > 0.15 & Boundary <= 0.25))
m34 = lm(ESBG_mean ~  log(N) + file * even_N + even_N * log(N), filter(tda, Boundary > 0.15 & Boundary <= 0.25))
stargazer::stargazer(m31, m32, m33, type = "text")

m41 = lm(extremness_mean ~  log(N) + file, filter(tda, Boundary > 0.15 & Boundary <= 0.25))
m42 = lm(extremness_mean ~  log(N) + file + even_N, filter(tda, Boundary > 0.15 & Boundary <= 0.25))
m43 = lm(extremness_mean ~  log(N) + file * even_N, filter(tda, Boundary > 0.15 & Boundary <= 0.25))
m44 = lm(extremness_mean ~  log(N) + file * even_N + even_N * log(N), filter(tda, Boundary > 0.15 & Boundary <= 0.25))
stargazer::stargazer(m41, m42, m43, type = "text")

m51 = lm(diversity_mean ~  log(N) + file, filter(tda, Boundary > 0.15 & Boundary <= 0.25))
m52 = lm(diversity_mean ~  log(N) + file + even_N, filter(tda, Boundary > 0.15 & Boundary <= 0.25))
m53 = lm(diversity_mean ~  log(N) + file * even_N, filter(tda, Boundary > 0.15 & Boundary <= 0.25))
m54 = lm(diversity_mean ~  log(N) + file * even_N + even_N * log(N), filter(tda, Boundary > 0.15 & Boundary <= 0.25))
stargazer::stargazer(m51, m52, m53, type = "text")




# How much differ odd and even Ns in different scenarios? -----------------

## Variability analysis -- something Veeeeeeeeeeeeeeeeeeery simple:
tda %>% glimpse()

# Firstly, without distinguishing for even/odd:
nodist = tda %>%
  filter(Boundary > 0.15 & Boundary <= 0.25) %>%
  group_by(file) %>%
  summarise(diversity_sd = sd(diversity_mean),
            extremness_sd = sd(extremness_mean), ESBG_sd = sd(ESBG_mean),
            diversity_mean = mean(diversity_mean),
            extremness_mean = mean(extremness_mean), ESBG_mean = mean(ESBG_mean))

# Secondly, without distinguishing for even/odd:
dist = tda %>%
  filter(Boundary > 0.15 & Boundary <= 0.25) %>%
  group_by(file, even_N) %>%
  summarise(diversity_sd = sd(diversity_mean),
            extremness_sd = sd(extremness_mean), ESBG_sd = sd(ESBG_mean),
            diversity_mean = mean(diversity_mean),
            extremness_mean = mean(extremness_mean), ESBG_mean = mean(ESBG_mean))
nodist
dist

# Thirdly, we try to compute differences for censecutive odd/even roughly same N:
diff = tda %>%
  group_by(file, Boundary, Conformity, even_N) %>%
  arrange(N) %>% mutate(order = row_number()) %>%
  group_by(file, Boundary, Conformity, order) %>%
  summarise(diversity = max(diversity_mean) - min(diversity_mean),
            extremness = max(extremness_mean) - min(extremness_mean),
            ESBG = max(ESBG_mean) - min(ESBG_mean)) %>%
  mutate(N_hi = order >= 6)

# Let's compute some numbers!
diff %>% group_by(file) %>%
  summarise(max_d = max(diversity), max_e = max(extremness), max_E = max(ESBG),
            diversity = mean(diversity), extremness = mean(extremness), ESBG = mean(ESBG))
diff %>% filter(Boundary >= 0.15 & Boundary <= 0.25) %>% group_by(file) %>%
  summarise(max_d = max(diversity), max_e = max(extremness), max_E = max(ESBG),
            diversity = mean(diversity), extremness = mean(extremness), ESBG = mean(ESBG))

# And lastly some maps again!
df = diff %>%
  # filter(Boundary > 0.15 & Boundary <= 0.25) %>%
  group_by(file, Boundary, Conformity) %>%
  summarise(diversity = mean(diversity), extremness = mean(extremness), ESBG = mean(ESBG))

.height = 46
.width = 34.5

df %>% heat_map_facet(.var = "diversity", .title = "Map of differences in DIVERSITY between consecutive even and odd numbers") %>%
  ggsave("Pics/mapb1.png", plot = ., units = "cm", height = .height, width = .width)
df %>% heat_map_facet(.var = "extremness", .title = "Map of differences in EXTREMNESS between consecutive even and odd numbers")  %>%
  ggsave("Pics/mapb2.png", plot = ., units = "cm", height = .height, width = .width)
df %>% heat_map_facet(.var = "ESBG", .title = "Map of differences in ESBG between consecutive even and odd numbers")  %>%
  ggsave("Pics/mapb3.png", plot = ., units = "cm", height = .height, width = .width)



# T-test ------------------------------------------------------------------

# Package for tidy testing
# https://cran.r-project.org/web/packages/infer/
install.packages("infer", dependencies = T)
library(infer)

# Summary statistics:
# Mean ESBG for simulation parameters
tda %>%
  get_summary_stats(ESBG_mean, type = "mean_sd")
# Differences for adjacent odd and even N of populations
diff %>% ungroup() %>%
  get_summary_stats(ESBG, type = "mean_sd")

# Mere T-tests:
# One-sample T-test:
tda %>%
  infer::t_test(response = ESBG_mean, mu = 0)

# Paired T-test:
diff %>%
  infer::t_test(response = ESBG, mu = 0)

# Two-sample T-test:
infer::t_test(x = tda, formula = ESBG_mean ~ even_N, order = c(T, F), alternative = "two-sided")



# Means for more than two groups: ANOVA -----------------------------------

#### (with relaxed assumptions on normality)
install.packages("rstatix", dependencies = T)
library(rstatix)
library(ggpubr)

# One-way ANOVA test
#:::::::::::::::::::::::::::::::::::::::::
# Summary statistics:
diff %>%
  group_by(file) %>%
  get_summary_stats(ESBG, type = "mean_sd")

# ANOVA
res.aov2 = diff %>% ungroup() %>%
  welch_anova_test(ESBG ~ file)
res.aov2

# Post-hoc test
pwc2 = diff %>% ungroup() %>% games_howell_test(ESBG ~ file)
pwc2

# Visualization: box plots with p-values
pwc2 = pwc2 %>% add_xy_position(x = "file", step.increase = .1)
ggboxplot((diff %>% ungroup() #%>% mutate(ESBG = ESBG + 0.00001)
           ),
          x = "file", y = "ESBG", color = "file", palette = "jco") +
  # scale_y_log10() +
  stat_pvalue_manual(pwc2, hide.ns = TRUE) +
  labs(
    subtitle = get_test_label(res.aov2, detailed = TRUE),
    caption = get_pwc_label(pwc2)
  )


## Alternative when ANOVA criteria are not met:
# Kruskal-Wallis test
res.kruskal = diff %>% ungroup() %>% kruskal_test(ESBG ~ file)
res.kruskal
diff %>% ungroup() %>% kruskal_effsize(ESBG ~ file)

# Pairwise comparisons:
# Dunn
pwc = diff %>% ungroup() %>%
  dunn_test(ESBG ~ file, p.adjust.method = "bonferroni")
pwc
# Další možnost je Wilcoxon, ale Dunn má pár vychytávek navíc,
# tak vám nebudu motat hlavu.


# Visualization: box plots with p-values
pwc = pwc %>% add_xy_position(x = "file", step = 1)
ggboxplot(data = (diff %>% ungroup() %>% mutate(ESBG = ESBG + 0.00001)),
          x = "file", y = "ESBG", color = "file", palette = "jco") +
  stat_pvalue_manual(pwc, hide.ns = TRUE) +
  scale_y_log10() +
  labs(
    subtitle = get_test_label(res.kruskal, detailed = TRUE),
    caption = get_pwc_label(pwc)
  )



# Grouped One-way ANOVA test
#:::::::::::::::::::::::::::::::::::::::::
# ANOVA
diff %>% ungroup() %>%
  group_by(N_hi) %>%
  welch_anova_test(ESBG ~ file)

# Grouped post-hoc test
post.hoc = diff %>%
  group_by(N_hi) %>%
  games_howell_test(ESBG ~ file)
post.hoc



# Chi-square --------------------------------------------------------------

# Categorization of data:
df = diff %>% mutate(ESBG_hi = if_else(ESBG > 0.05, "High", "Low")) %>% ungroup()

# Simple counts:
count(df, ESBG_hi, file)

# Table of counts and percents:
count(df, ESBG_hi, file) %>%
  pivot_wider(id_cols = file, names_from = ESBG_hi, values_from = n) %>%
  rowwise() %>%
  mutate(suma = sum(c_across(cols = 2:3)) ,
         across(.cols = 2:4, ~.x / suma * 100, .names = "{.col}_perc"))

# Basic wrapper:
infer::chisq_test(df, ESBG_hi ~ file)




# Two-way ANOVA test ------------------------------------------------------

#:::::::::::::::::::::::::::::::::::::::::
# Summary statistics:
diff %>%
  group_by(N_hi, file) %>%
  get_summary_stats(ESBG, type = "mean_sd")

# Visualisation
bxp = ggboxplot(
  diff, x = "N_hi", y = "ESBG",
  color = "file", palette = "jco"
)
bxp


#### Check assumptions
## outliers/extremes:
diff %>%
  group_by(file, N_hi) %>%
  identify_outliers(ESBG)

## Normality:
# Build the linear model
model = lm(ESBG ~ file * N_hi,
           data = (diff %>% ungroup() %>%  sample_n_by(data = ., file, N_hi, size = 600)))
# Create a QQ plot of residuals
ggqqplot(residuals(model))

# Compute Shapiro-Wilk test of normality for whole sample:
shapiro_test(residuals(model))

# Compute it for groups:
diff %>%
  group_by(file, N_hi) %>%
  shapiro_test(ESBG) %>% pivot_wider(id_cols = c(N_hi), names_from = file, values_from = p)

# Vizualise:
ggqqplot(diff, "ESBG", ggtheme = theme_light()) +
  facet_grid(N_hi ~ file)


## Homogeneity of Variance assumption:
diff %>% ungroup() %>% levene_test(ESBG ~ file * N_hi)


## Two-way ANOVA itself:
res.aov = diff %>% ungroup() %>% anova_test(ESBG ~ file * N_hi)
res.aov

## Post-hoc tests:
# Simple main effect:
model = lm(ESBG ~ file * N_hi, data = diff)
diff %>%
  group_by(file) %>%
  anova_test(ESBG ~ N_hi, error = model)
diff %>%
  group_by(N_hi) %>%
  anova_test(ESBG ~ file, error = model)

# pairwise comparisons:
library(emmeans)
pwc = diff %>%
  group_by(N_hi) %>%
  emmeans_test(ESBG ~ file, p.adjust.method = "bonferroni", model = model)
pwc


# Visualization: box plots with p-values
pwc = pwc %>% add_xy_position(x = "order")
bxp +
  stat_pvalue_manual(pwc) +
  labs(
    subtitle = get_test_label(res.aov, detailed = TRUE),
    caption = get_pwc_label(pwc)
  )



# Correlation -------------------------------------------------------------

diff %>% cor_mat(diversity, extremness, ESBG) %>% cor_gather()
