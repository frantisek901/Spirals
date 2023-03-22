#### Test of measures of polarization
####

library(tidyverse)

tb1 = tibble(x = seq(0.01, .99, 0.01)) %>%
  mutate(
    dissonance = if_else(x <= .5, (x ^ 3) / ((1 - x) ^ 2), ((1 - x) ^ 3) / (x ^ 2)),
    ratio = x * (1 - x) * 2,
    entropy = 1 - ((log2(x) + log2(1 - x)) / (log2(0.01) + log2(0.99))),
    mean = (dissonance + entropy + ratio) / 3
  ) %>%
  pivot_longer(cols = 2:5, names_to = "Measures")

tb1 %>%
  ggplot() +
  aes(x = x, y = value, group = Measures, col = Measures) +
  geom_line(linewidth = 2) +
  theme_classic() +
  theme(legend.position = "bottom")



tb2 = tibble(x = rep(seq(0.01, .99, 0.01), times = 99), y = rep(seq(0.01, .99, 0.01), each = 99)) %>%
  filter((x + y) <= 1) %>%
  mutate(
    dissonance = if_else(x <= y, (x ^ 3) / (y ^ 2), (y ^ 3) / (x ^ 2)),
    ratio = x * y / (x + y) * 2,
    entropy = 1 - ((log2(x) + log2(y)) / (log2(0.01) + log2(x + y - 0.01))),
    mean = (dissonance + entropy + ratio) / 3
  ) %>%
  pivot_longer(cols = 3:6, names_to = "Measures")

tb2 %>%
  ggplot() +
  aes(x = x, y = y, z = value) +
  facet_wrap(vars(Measures)) +
  geom_contour_filled() +
  scale_fill_viridis_d() +
  theme_classic()

