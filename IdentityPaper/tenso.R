#### Script for transforming opinion positions to vector/matrix/tensor and
#### then for measuring fractal dimension and entropy of agents' positions

## Encoding: windows-1250
## Created:  2022-12-25 FranÈesko
## Edited:   2022-01-03 FranÈesko

## NOTES:
#  1) For more dimensions it's better to pivot tibbles on agents, not tiles/patches --
#  we have 101 agents, so only 10x10 tiles in 2D, 3x3x3x3 in 4D are more efficient
#  than tibble with 101 rows for 101 agents and encoding their positions in the space.
#  And we want to use 360 tiles per dimension.
#


# Head --------------------------------------------------------------------

# Clearing:
rm(list = ls())

# Packages:
library(stringr)
library(tidyverse)
library(rstatix)



# Functions ---------------------------------------------------------------

compute_basic_tibble = function(.tibble, .ratio = 2.002, .max_tile = 360, .delta = 1.0005) {
  # We suppose that:
  # 1) raw data are of class tibble
  # 2) in variables are coded only agents' positions in the opinion space
  # 3) each variable codes position of agents in one dimension
  # 4) each row is a complete position of an respective agent
  # 5) nothing else than opinion in the opinion space is in the tibble
  # 6) opinions are on the continuous scale -1; +1, step is 0.001 (precision 3/ round(., 3))
  #
  .tibble %>%
    mutate(across(.cols = everything(), ~ceiling((.x + .delta) / .ratio * .max_tile)))
}

compute_row = function(.tibble, .reduced_tiles = 360, .original_tiles = 360) {
  # We suppose that:
  # 1) basic data are of class tibble
  # 2) in variables are coded only agents' positions in the opinion space
  # 3) each variable codes position of agents in one dimension
  # 4) each row is a complete position of an respective agent
  # 5) nothing else than opinion in the opinion space is in the tibble
  # 6) basic opinions are on the discrete scale 1; 2; ... 360, step is 1
  #

  # Firstly, we must compute useful constants:
  # a) factor by which the basic tiles will be reduced:
  .reducing_factor = .original_tiles / .reduced_tiles
  # b) number of dimensions:
  .dims = ncol(.tibble)

  # Secondly, we use this factor for reducing tile numbers
  .tibble %>%
    mutate(across(.cols = everything(), ~ceiling(.x / .reducing_factor))) %>%

    # Thirdly, let's compute entropy!
    count(across()) %>%  # Finding all non-zero combinations
    mutate(  # Computing itself
      p = n / sum(n),
      h = -1 * p * log(p, base = if_else(nrow(.tibble) < (.reduced_tiles ^ .dims),
                                         nrow(.tibble),
                                         as.integer(.reduced_tiles ^ .dims)))) %>%

    # Fourthly, let' summarize everything to row and
    # add variables needed for computing fractal dimension:
    summarise(h = sum(h), tiles_covered = n()) %>%
    mutate(dimensions = .dims, tiles_size = as.integer(.reducing_factor),
           tiles_max = as.integer(.reduced_tiles)) %>%
    relocate(dimensions, .before = h)
}

compute_final_tibble = function(.tibble, .sizes = c(360, 180, 120, 90, 72, 60, 45, 40, 36, 30, 24, 20,
                                                    18, 15, 12, 10, 9, 8, 6, 5, 4, 3, 2)) {
  # We prepare the first row of final tibble:
  .df = compute_row(.tibble, .sizes[1])

  # We use FOR cycle for computing 'h' and fractal variables for all tiles' sizes:
  for (.s in .sizes[2:length(.sizes)]) {
    .df = .df %>%
      add_row(compute_row(.tibble, .s))
  }

  # Publishing resulting file:
  .df
}

compute_regressions = function(.tibble, .min_tile_size = 1, .max_tile_size = 180) {
  # We suppose that '.tibble' is the product of the function compute_final_tibble()
  # We also suppose that main produced results will be updated, i.e. what this function produces.

  # Here we compute regressions -- always once for 'tiles_max' divisible by 5 and not
  ho = lm(h ~ log(tiles_max), data = filter(.tibble, tiles_max%%5 == 0 & tiles_size >= .min_tile_size & tiles_size <= .max_tile_size))
  he = lm(h ~ log(tiles_max), data = filter(.tibble, tiles_max%%5 != 0 & tiles_size >= .min_tile_size & tiles_size <= .max_tile_size))
  hdo = lm(log(tiles_covered) ~ log(1 / tiles_size), data = filter(.tibble, tiles_max%%5 == 0 & tiles_size >= .min_tile_size & tiles_size <= .max_tile_size))
  hde = lm(log(tiles_covered) ~ log(1 / tiles_size), data = filter(.tibble, tiles_max%%5 != 0 & tiles_size >= .min_tile_size & tiles_size <= .max_tile_size))

  # Finally, we publish results as a tibble:
  tibble(h_odd = ho$coefficients[2], h_even = he$coefficients[2],
         r_h_odd = summary(ho)$adj.r.squared, r_h_even = summary(he)$adj.r.squared,
         fd_odd = hdo$coefficients[2], fd_even = hde$coefficients[2],
         r_fd_odd = summary(hdo)$adj.r.squared, r_fd_even = summary(hde)$adj.r.squared)
}

compute_entropy = function(.tibble, .min_tile_size = 1, .max_tile_size = 180) {
  # We suppose that '.tibble' is the product of the function compute_final_tibble()

  # Computing regression:
  mc1 = lm(h ~ log(tiles_max),
           data = filter(.tibble, tiles_size >= .min_tile_size & tiles_size <= .max_tile_size))
  ms1 = mc1 %>% summary()
  mc2 = lm(h ~ log(tiles_max) * (tiles_max%%5 == 0) ,
           data = filter(.tibble, tiles_size >= .min_tile_size & tiles_size <= .max_tile_size))
  ms2 = mc2 %>% summary()

  # Printing results:
  print("Simple model:")
  ms1$coefficients %>% print()
  print(paste0("R2: ", round(ms1$adj.r.squared, 3)))
  print(""); print("")
  print("Model with interaction:")
  ms2$coefficients %>% print()
  print(paste0("R2: ", round(ms2$adj.r.squared, 3)))
  print(""); print("")
  print("Difference in BIC: model with interaction - simple model")
  round(BIC(mc2) - BIC(mc1), 1) %>% print()
  # # Here we compute regression
  # lm(h ~ log(tiles_max) * (tiles_max%%5 == 0) ,
  #    data = filter(.tibble, tiles_size >= .min_tile_size)) %>%
  #   summary()
}
# compute_entropy(tn4)

compute_fractal_dimension = function(.tibble, .min_tile_size = 1, .max_tile_size = 360) {
  # We suppose that '.tibble' is the product of the function compute_final_tibble()

  # Here we compute regressions -- always once for 'tiles_max' divisible by 5 and not
  mc1 = lm(log(tiles_covered) ~ log(1 / tiles_size),
           data = filter(.tibble, tiles_size >= .min_tile_size & tiles_size <= .max_tile_size))
  ms1 = mc1 %>% summary()
  mc2 = lm(log(tiles_covered) ~ log(1 / tiles_size) * (tiles_max%%5 == 0) ,
           data = filter(.tibble, tiles_size >= .min_tile_size & tiles_size <= .max_tile_size))
  ms2 = mc2 %>% summary()

  # And here we present them:
  print("Simple model:")
  ms1$coefficients %>% print()
  print(paste0("R2: ", round(ms1$adj.r.squared, 3)))
  print(""); print("")
  print("Model with interaction:")
  ms2$coefficients %>% print()
  print(paste0("R2: ", round(ms2$adj.r.squared, 3)))
  print(""); print("")
  print("Difference in BIC: model with interaction - simple model")
  round(BIC(mc2) - BIC(mc1), 1)
}
# compute_fractal_dimension(tn4, 1, 90)


# Generating test files -------------------------------

tu4 = tibble(
  o1 = runif(101, min = -.99, max = .99) %>% round(3),
  o2 = runif(101, min = -.99, max = .99) %>% round(3),
  o3 = runif(101, min = -.99, max = .99) %>% round(3),
  o4 = runif(101, min = -.99, max = .99) %>% round(3)
) %>% compute_basic_tibble() %>% compute_final_tibble()
tu4

tn4 = tibble(
  o1 = c(rnorm(50, mean = 0.41, sd = 0.2), rnorm(51, mean = -0.42, sd = 0.02)) %>% round(3),
  o2 = c(rnorm(50, mean = 0.41, sd = 0.2), rnorm(51, mean = -0.42, sd = 0.02)) %>% round(3),
  o3 = c(rnorm(50, mean = 0.41, sd = 0.2), rnorm(51, mean = -0.42, sd = 0.02)) %>% round(3),
  o4 = c(rnorm(50, mean = 0.41, sd = 0.2), rnorm(51, mean = -0.42, sd = 0.02)) %>% round(3)
) %>% compute_basic_tibble() %>% compute_final_tibble()
tn4

tu3 = tibble(
  o1 = runif(101, min = -.99, max = .99) %>% round(3),
  o2 = runif(101, min = -.99, max = .99) %>% round(3),
  o3 = runif(101, min = -.99, max = .99) %>% round(3)
) %>% compute_basic_tibble() %>% compute_final_tibble()
tu3

tn3 = tibble(
  o1 = c(rnorm(50, mean = 0.41, sd = 0.2), rnorm(51, mean = -0.42, sd = 0.02)) %>% round(3),
  o2 = c(rnorm(50, mean = 0.41, sd = 0.2), rnorm(51, mean = -0.42, sd = 0.02)) %>% round(3),
  o3 = c(rnorm(50, mean = 0.41, sd = 0.2), rnorm(51, mean = -0.42, sd = 0.02)) %>% round(3)
) %>% compute_basic_tibble() %>% compute_final_tibble()
tn3

tu2 = tibble(
  o1 = runif(101, min = -.99, max = .99) %>% round(3),
  o2 = runif(101, min = -.99, max = .99) %>% round(3)
) %>% compute_basic_tibble() %>% compute_final_tibble()
tu2

tn2 = tibble(
  o1 = c(rnorm(50, mean = 0.41, sd = 0.2), rnorm(51, mean = -0.42, sd = 0.02)) %>% round(3),
  o2 = c(rnorm(50, mean = 0.41, sd = 0.2), rnorm(51, mean = -0.42, sd = 0.02)) %>% round(3)
) %>% compute_basic_tibble() %>% compute_final_tibble()
tn2

tu1 = tibble(
  o1 = runif(101, min = -.99, max = .99) %>% round(3)
) %>% compute_basic_tibble() %>% compute_final_tibble()
tu1

tn1 = tibble(
  o1 = c(rnorm(50, mean = 0.41, sd = 0.2), rnorm(51, mean = -0.42, sd = 0.02)) %>% round(3)
) %>% compute_basic_tibble() %>% compute_final_tibble()
tn1

res = tn4



# Regressions -------------------------------------------------------------

mt = lm(h ~ log(tiles_max) * (tiles_max%%5 == 0) , data = filter(res, tiles_size >= 1))
mt %>% summary()
mc2 = lm(log(tiles_covered) ~ log(1 / tiles_size) * (tiles_max%%5 == 0) , data = filter(res, tiles_size >= 1))
mc2 %>% summary()
mc1 = lm(log(tiles_covered) ~ log(1 / tiles_size), data = filter(res, tiles_size >= 1))
mc1 %>% summary()
BIC(mc2) - BIC(mc1)

compute_regressions(tn2)
compute_regressions(res, 2, 90)
compute_entropy(tn2)
compute_fractal_dimension(tu2, 6)



# Graphs with entropy and fractal dimension -------------------------------

# Shannon's entropy
res = tn1
res %>% filter(tiles_size >= 2, tiles_size <= 90) %>%
  ggplot() +
  aes(x = tiles_max, y = h, group = tiles_max%%5 == 0, color = tiles_max%%5 == 0) +
  geom_line() +
  geom_point(size = 3, alpha = 1) +
  scale_x_log10(breaks = c(2:6, 8:10, 12, 15, 18, 20, 24, 30, 36, 40, 45, 60, 72, 90, 120, 180, 360),
                minor_breaks = NULL) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.1)) +
  labs(x = "Number of tiles", y = "h (normalized entropy)") +
  theme_light()
compute_entropy(res, 2, 90)


# Hausdorff's fractal dimension
res = tn1
res %>% filter(tiles_size >= 2, tiles_size <= 90) %>%
  ggplot() +
  aes(x = (tiles_size), y = tiles_covered, group = tiles_max%%5 == 0, color = tiles_max%%5 == 0) +
  geom_line() +
  geom_point(size = 3) +
  scale_x_log10(breaks = c(1:6, 8:10, 12, 15, 18, 20, 24, 30, 36, 40, 45, 60, 72, 90, 120, 180),
                minor_breaks = NULL) +
  scale_y_log10() +
  theme_light()
compute_fractal_dimension(res, 2, 90)
