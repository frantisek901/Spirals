#### Script for transforming opinion positions to vector/matrix/tensor and
#### then for measuring fractal dimension and entropy of agents' positions

## Encoding: windows-1250
## Created:  2022-12-25 FranÈesko
## Edited:   2022-12-28 FranÈesko


# Head --------------------------------------------------------------------

# Clearing:
rm(list = ls())

# Packages:
library(stringr)


# 1D ----------------------------------------------------------------------

# Computing basic vector 't' of 360 tiles
v = seq(-1, 1, 0.001)
l = length(v) / 1000
l = 2.002  # 2.002 gives tile 181 and tile 180 of equal size (n=6).
#l = 2.001 # 2.001 gives tiles 1 and 360 of equal size (n=6).
v_trans = v + 1.0005
t_max = 360
#t = ceiling((v + 1.0005) / 2.002 * 360)
t = ceiling(v_trans / l * t_max)


# Checking 't'
max(t)
t[t == 180]
t[t == 181]
t[t == max(t)]
t[t == min(t)]
hist(t, breaks = 0:t_max)
tibble(t) %>% count(t) %>%
  ggplot(aes(x = n)) + geom_bar() + labs(title = "Frequencies of tiles") + theme_light()



# Testing 1D equation -----------------------------------------------------

te = c(rnorm(50, mean = 0.41, sd = 0.2),
       # rnorm(1, mean = -0.9, sd = 0.002),
       rnorm(51, mean = -0.42, sd = 0.02)) %>% round(3)
te = runif(101, min = -.99, max = .99) %>% round(3)
te
hist(x = te, breaks = seq(-1, 1, (0.0005 / 2.002 * 360)))
min(te)
max(te)
# Converting according the equation...
te360 = tibble(t = ceiling((te + 1.0005) / 2.002 * 360))
hist(x = te360$t, breaks = 0:360)



# Computing re-scaled 1D vectors ------------------------------------------

t360 = te360 %>%
  count(t) %>%
  # Computing Shannon's entropy
  mutate(p = n / sum(n),
         h = -1 * p * log(p, base = if_else(length(te360$t) < 360, length(te360$t), as.integer(360))))
sum(t360$h)


# Converting according the equation...
t180 = tibble(t = ceiling(te360$t / 2)) %>%
  count(t) %>%
  # Computing Shannon's entropy
  mutate(p = n / sum(n),
         h = -1 * p * log(p, base = if_else(length(te360$t) < 180, length(te360$t), as.integer(180))))
sum(t180$h)


# Converting according the equation...
t120 = tibble(t = ceiling(te360$t / 3)) %>%
  count(t) %>%
  # Computing Shannon's entropy
  mutate(p = n / sum(n),
         h = -1 * p * log(p, base = if_else(length(te360$t) < 120, length(te360$t), as.integer(120))))
sum(t120$h)


# Converting according the equation...
t090 = tibble(t = ceiling(te360$t / 4)) %>%
  count(t) %>%
  # Computing Shannon's entropy
  mutate(p = n / sum(n),
         h = -1 * p * log(p, base = if_else(length(te360$t) < 90, length(te360$t), as.integer(90))))
sum(t090$h)


# Converting according the equation...
t072 = tibble(t = ceiling(te360$t / 5)) %>%
  count(t) %>%
  # Computing Shannon's entropy
  mutate(p = n / sum(n),
         h = -1 * p * log(p, base = if_else(length(te360$t) < 72, length(te360$t), as.integer(72))))
sum(t072$h)


# Converting according the equation...
t060 = tibble(t = ceiling(te360$t / 6)) %>%
  count(t) %>%
  # Computing Shannon's entropy
  mutate(p = n / sum(n),
         h = -1 * p * log(p, base = if_else(length(te360$t) < 60, length(te360$t), as.integer(60))))
sum(t060$h)


# Converting according the equation...
t045 = tibble(t = ceiling(te360$t / 8)) %>%
  count(t) %>%
  # Computing Shannon's entropy
  mutate(p = n / sum(n),
         h = -1 * p * log(p, base = if_else(length(te360$t) < 45, length(te360$t), as.integer(45))))
sum(t045$h)


# Converting according the equation...
t040 = tibble(t = ceiling(te360$t / 9)) %>%
  count(t) %>%
  # Computing Shannon's entropy
  mutate(p = n / sum(n),
         h = -1 * p * log(p, base = if_else(length(te360$t) < 40, length(te360$t), as.integer(40))))
sum(t040$h)


# Converting according the equation...
t036 = tibble(t = ceiling(te360$t / 10)) %>%
  count(t) %>%
  # Computing Shannon's entropy
  mutate(p = n / sum(n),
         h = -1 * p * log(p, base = if_else(length(te360$t) < 36, length(te360$t), as.integer(36))))
sum(t036$h)


# Converting according the equation...
t030 = tibble(t = ceiling(te360$t / 12)) %>%
  count(t) %>%
  # Computing Shannon's entropy
  mutate(p = n / sum(n),
         h = -1 * p * log(p, base = if_else(length(te360$t) < 30, length(te360$t), as.integer(30))))
sum(t030$h)


# Converting according the equation...
t024 = tibble(t = ceiling(te360$t / 15)) %>%
  count(t) %>%
  # Computing Shannon's entropy
  mutate(p = n / sum(n),
         h = -1 * p * log(p, base = if_else(length(te360$t) < 24, length(te360$t), as.integer(24))))
sum(t024$h)


# Converting according the equation...
t020 = tibble(t = ceiling(te360$t / 18)) %>%
  count(t) %>%
  # Computing Shannon's entropy
  mutate(p = n / sum(n),
         h = -1 * p * log(p, base = if_else(length(te360$t) < 20, length(te360$t), as.integer(20))))
sum(t020$h)


# Converting according the equation...
t018 = tibble(t = ceiling(te360$t / 20)) %>%
  count(t) %>%
  # Computing Shannon's entropy
  mutate(p = n / sum(n),
         h = -1 * p * log(p, base = if_else(length(te360$t) < 18, length(te360$t), as.integer(18))))
sum(t018$h)


# Converting according the equation...
t015 = tibble(t = ceiling(te360$t / 24)) %>%
  count(t) %>%
  # Computing Shannon's entropy
  mutate(p = n / sum(n),
         h = -1 * p * log(p, base = if_else(length(te360$t) < 15, length(te360$t), as.integer(15))))
sum(t015$h)


# Converting according the equation...
t012 = tibble(t = ceiling(te360$t / 30)) %>%
  count(t) %>%
  # Computing Shannon's entropy
  mutate(p = n / sum(n),
         h = -1 * p * log(p, base = if_else(length(te360$t) < 12, length(te360$t), as.integer(12))))
sum(t012$h)


# Converting according the equation...
t010 = tibble(t = ceiling(te360$t / 36)) %>%
  count(t) %>%
  # Computing Shannon's entropy
  mutate(p = n / sum(n),
         h = -1 * p * log(p, base = if_else(length(te360$t) < 10, length(te360$t), as.integer(10))))
sum(t010$h)


# Converting according the equation...
t009 = tibble(t = ceiling(te360$t / 40)) %>%
  count(t) %>%
  # Computing Shannon's entropy
  mutate(p = n / sum(n),
         h = -1 * p * log(p, base = if_else(length(te360$t) < 9, length(te360$t), as.integer(9))))
sum(t009$h)


# Converting according the equation...
t008 = tibble(t = ceiling(te360$t / 45)) %>%
  count(t) %>%
  # Computing Shannon's entropy
  mutate(p = n / sum(n),
         h = -1 * p * log(p, base = if_else(length(te360$t) < 8, length(te360$t), as.integer(8))))
sum(t008$h)


# Converting according the equation...
t006 = tibble(t = ceiling(te360$t / 60)) %>%
  count(t) %>%
  # Computing Shannon's entropy
  mutate(p = n / sum(n),
         h = -1 * p * log(p, base = if_else(length(te360$t) < 6, length(te360$t), as.integer(6))))
sum(t006$h)


# Converting according the equation...
t005 = tibble(t = ceiling(te360$t / 72)) %>%
  count(t) %>%
  # Computing Shannon's entropy
  mutate(p = n / sum(n),
         h = -1 * p * log(p, base = if_else(length(te360$t) < 5, length(te360$t), as.integer(5))))
sum(t005$h)


# Converting according the equation...
t004 = tibble(t = ceiling(te360$t / 90)) %>%
  count(t) %>%
  # Computing Shannon's entropy
  mutate(p = n / sum(n),
         h = -1 * p * log(p, base = if_else(length(te360$t) < 4, length(te360$t), as.integer(4))))
sum(t004$h)


# Converting according the equation...
t003 = tibble(t = ceiling(te360$t / 120)) %>%
  count(t) %>%
  # Computing Shannon's entropy
  mutate(p = n / sum(n),
         h = -1 * p * log(p, base = if_else(length(te360$t) < 3, length(te360$t), as.integer(3))))
sum(t003$h)


# Converting according the equation...
t002 = tibble(t = ceiling(te360$t / 180)) %>%
  count(t) %>%
  # Computing Shannon's entropy
  mutate(p = n / sum(n),
         h = -1 * p * log(p, base = if_else(length(te360$t) < 2, length(te360$t), as.integer(2))))
sum(t002$h)




# Graphs with entropy and fractal dimension -------------------------------

# Data preparation -- computing results from rescaled vectors
res = tibble(t = NA_real_, h = NA_real_, covered = NA_real_, size = NA_real_)
for (f in c("t002", "t003", "t004", "t005", "t006", "t008", "t009", "t010", "t012", "t015", "t018",
            "t020", "t024", "t030", "t036", "t040", "t045", "t060", "t072", "t090", "t120", "t180", "t360")) {
  x = tibble(
    t = str_sub(f, start = 2) %>% as.numeric(),
    h = eval(as.name(f))$h %>% sum(),
    covered = nrow(eval(as.name(f))),
    size = 360 / t)
  res = res %>% add_row(x) %>% filter(!is.na(t))
}
res

# Regressions
mt = lm(h ~ log(t) * (t%%5 == 0) , data = filter(res, t >= 6))
mt %>% summary()
mc2 = lm(log(covered) ~ log(1 / size) * (t%%5 == 0) , data = filter(res, t >= 6))
mc2 %>% summary()
mc1 = lm(log(covered) ~ log(1 / size), data = filter(res, t >= 6))
mc1 %>% summary()
BIC(mc2) - BIC(mc1)

## Graphs
# Shannon's entropy
res %>%
  ggplot() +
  aes(x = t, y = h, group = t%%5 == 0, color = t%%5 == 0) +
  geom_line() +
  geom_point(size = 3, alpha = 1) +
  scale_x_log10(breaks = c(2:6, 8:10, 12, 15, 18, 20, 24, 30, 36, 40, 45, 60, 72, 90, 120, 180, 360),
                minor_breaks = NULL) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.1)) +
  labs(x = "Number of tiles", y = "h (normalized entropy)") +
  theme_light()

# Hausdorff's fractal dimension
res %>% filter(t >= 6) %>%
  ggplot() +
  aes(x = (size), y = covered, group = t%%5 == 0, color = t%%5 == 0) +
  geom_line() +
  geom_point(size = 3) +
  scale_x_log10(breaks = c(1:6, 8:10, 12, 15, 18, 20, 24, 30, 36, 40, 45, 60),
                minor_breaks = NULL) +
  scale_y_log10() +
  theme_light()

