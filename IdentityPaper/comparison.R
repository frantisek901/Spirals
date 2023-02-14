#### Script for processing data from Identity paper experiments
#### Now we combine data from different steps.
####


## Encoding: windows-1250
## Created:  2023-02-03 FrK
## Edited:   2023-02-03 FrK

## Notes:
## 1) We suppose that scripts from each step were processed
##
## 2) We have to simulate 'Boundary' for Step 4.1 for values 0.11:0.19 and 0.21:0.29,
##    because it is needed for comparisons with previous steps
##



# Processing scripts ------------------------------------------------------

# source('classicalHK.R')
# source('step2.R')
# source('step3.R')
# source('step4.R')



# Head --------------------------------------------------------------------

# Packages:
library(tidyverse)
library(rstatix)



# Combining data-files -----------------------------------------------------

# We start by joining files -- we do just changes needed for joining.
tc = ts41 %>%
  mutate(
    across(c(HK_distribution, Use_identity), ~factor(.x)),
    Step = "Normal SPIRO, Normal Boundary") %>%
  add_row(ts30 %>%
            mutate(
              N = N %>% as.character() %>% as.numeric(),
              SPIRO_STD = factor(0),
              across(c(HK_distribution, Use_identity, SPIRO_Mean), ~factor(.x)),
              Step = "Constant SPIRO, Normal Boundary")) %>%
  add_row(select(ts20, -Present_opinion) %>%
            mutate(Step = "No Identity, Normal Boundary",
                   Identity_Type = "none", SPIRO_Distribution = "none",
                   SPIRO_STD = factor(0), SPIRO_Mean = factor(0.25))) %>%
  add_row(select(ts10, -Present_opinion) %>%
            mutate(
              Step =  if_else(HK_distribution == "TRUE", "Classical HK" ,"No Identity, Constant Boundary"),
              Identity_Type = "none", SPIRO_Distribution = "none",
              Boundary_STD = 0, Conformity_STD = 0, RS = 0,
              SPIRO_STD = factor(0),SPIRO_Mean = factor(0.25))) %>%

  # Now we do changes on joined file:
  mutate(
    Step = factor(Step) %>% fct_inorder() %>% fct_rev() %>%  fct_relevel("Classical HK")) %>%

  # Selecting out obsolete variables:
  select(-(17:31)) %>%

  # Filtering only for comparable results:
  filter(Conformity %in% c(.2, .8), Conformity_STD %in% c(0, .1),
         as.numeric(as.character(SPIRO_Mean)) %in% seq(0.25, 0.85, 0.12),
         Boundary_STD %in% seq(0, 0.15, .05), Boundary >= 0.1, Boundary <= 0.3) #%>%
  #filter(Boundary %in% c(.1, .2, .3))

# Checking the final file:
tc %>% count(HK_distribution, Use_identity, Identity_Type, Step) %>% arrange(desc(n))
tc %>% count(fct_rev(Step), Boundary)
tc %>% count(fct_rev(Step))
tc %>% count(RS) %>% tail(20)
tc %>% filter(Step == "Normal SPIRO, Normal Boundary") %>% count(RS) %>%
  filter((RS %% 10) == 0) %>% arrange(n) %>% tail(20)
tc %>% filter(Step == "Normal SPIRO, Normal Boundary", Boundary %in% c(0.11, 0.12, 0.13)) %>%
  count(RS) %>% arrange((n)) %>% head()
tc %>% filter(Step == "Normal SPIRO, Normal Boundary", Boundary %in% c(0.14, 0.15, 0.16)) %>%
  count(RS) %>% arrange((n)) %>% head()
tc %>% filter(Step == "Normal SPIRO, Normal Boundary", Boundary %in% c(0.17, 0.18, 0.19)) %>%
  count(RS) %>% arrange((n)) %>% head()
tc %>% filter(Step == "Normal SPIRO, Normal Boundary", Boundary %in% c(0.21, 0.22, 0.23)) %>%
  count(RS) %>% arrange((n)) %>% head()
tc %>% filter(Step == "Normal SPIRO, Normal Boundary", Boundary %in% c(0.24, 0.25, 0.26)) %>%
  count(RS) %>% arrange((n)) %>% head()
tc %>% filter(Step == "Normal SPIRO, Normal Boundary", Boundary %in% c(0.27, 0.28, 0.29)) %>%
  count(RS) %>% arrange((n)) %>% head()



# Steps' comparison -------------------------------------------------------

# Comparison of summary statistics:
tc %>% # filter(Boundary %in% c(.1, .2, .3)) %>%
  group_by(Step) %>% get_summary_stats(ESBG)
tc %>% # filter(Boundary %in% c(.1, .2, .3)) %>%
  group_by(Step) %>% get_summary_stats(diversity)
tc %>% # filter(Boundary %in% c(.1, .2, .3)) %>%
  group_by(Step) %>% get_summary_stats(extremness)


# Comparison of summary statistics -- all steps mapped to just one run of classy HK:
tc %>% # filter(Boundary %in% c(.1, .2, .3)) %>%
  group_by(Boundary, Conformity, N, Step) %>% get_summary_stats(ESBG) %>% knitr::kable()
tc %>% # filter(Boundary %in% c(.1, .2, .3)) %>%
  group_by(Boundary, Conformity, N, Step) %>% get_summary_stats(diversity) %>% knitr::kable()
tc %>% # filter(Boundary %in% c(.1, .2, .3)) %>%
  group_by(Boundary, Conformity, N, Step) %>% get_summary_stats(extremness) %>% knitr::kable()


# Series of t-test -- whether do Steps 3 & 4 differ:
tc %>% filter(Step %in% c("Constant SPIRO, Normal Boundary", "Normal SPIRO, Normal Boundary")) %>%
  mutate(Step = factor(Step)) %>%
  t_test(formula = ESBG ~ Step, detailed = T)

results_34 = tc %>%
  filter(Step %in% c("Constant SPIRO, Normal Boundary", "Normal SPIRO, Normal Boundary"),
         !(Identity_Type == "individual" & SPIRO_STD %in% c(0.05, 0))) %>%
  mutate(Step = factor(Step)) %>%
  group_by(Boundary, Conformity, N) %>%
  t_test(formula = ESBG ~ Step, detailed = T) %>%
  arrange(estimate) #%>% knitr::kable()



# Series of t-test -- whether do Steps 1, 2, 3 & 4 differ:
tc %>% filter(Step != "Classical HK") %>%
  mutate(Step = factor(Step)) %>%
  t_test(formula = ESBG ~ Step)

results_1234 = tc %>% filter(Step != "Classical HK") %>%
  mutate(Step = factor(Step)) %>%
  group_by(Boundary, Conformity, N) %>%
  t_test(formula = ESBG ~ Step, detailed = T) %>%
  arrange(estimate) #%>% knitr::kable()



# Graphs: ESBG
tc %>%
  ggplot() +
  aes(x = ESBG) +
  facet_grid(rows = vars(Step), scales = "free_y") +
  geom_histogram(binwidth = 0.025, alpha = 0.7, fill = "steelblue") +
  theme_light()

# Graphs: diversity
tc %>%
  ggplot() +
  aes(x = diversity) +
  facet_grid(rows = vars(Step), scales = "free_y") +
  geom_histogram(binwidth = 0.025, alpha = 0.7, fill = "steelblue") +
  theme_light()

# Graphs: extremness
tc %>%
  ggplot() +
  aes(x = extremness) +
  facet_grid(rows = vars(Step), scales = "free_y") +
  geom_histogram(binwidth = 0.025, alpha = 0.7, fill = "steelblue") +
  theme_light()



# Boxplot: ESBG
tc %>%
  # filter(Boundary %in% c(.1, .2, .3)) %>%
  ggplot() +
  aes(y = Step, x = ESBG, fill = Step) +
  facet_grid(rows = vars(Conformity, N), cols = vars(Boundary), scales = "free_y") +
  geom_boxplot(alpha = 0.7) +
  guides(fill = "none") +
  theme_light()

# Boxplot: diversity
tc %>%
  # filter(Boundary %in% c(.1, .2, .3)) %>%
  ggplot() +
  aes(y = Step, x = diversity, fill = Step) +
  facet_grid(rows = vars(Conformity, N), cols = vars(Boundary), scales = "free_y") +
  geom_boxplot(alpha = 0.7) +
  guides(fill = "none") +
  theme_light()

# Boxplot: extremness
tc %>%
  # filter(Boundary %in% c(.1, .2, .3)) %>%
  ggplot() +
  aes(y = Step, x = extremness, fill = Step) +
  facet_grid(rows = vars(Conformity, N), cols = vars(Boundary), scales = "free_y") +
  geom_boxplot(alpha = 0.7) +
  guides(fill = "none") +
  theme_light()



# Data focused on differences ---------------------------------------------

tcd = tc %>%
  group_by(Conformity, N, Boundary) %>%
  mutate(
    X = if_else(Step == "Classical HK", ESBG, NA_real_),
    Difference = ESBG - mean(X, na.rm = TRUE)) %>% ungroup()

tcd %>% filter(Step != "Classical HK") %>%
  ggplot() +
  aes(x = Difference) +
  facet_grid(rows = vars(Step), scales = "free_y") +
  geom_histogram(binwidth = 0.025, alpha = 0.7, fill = "steelblue") +
  theme_light()

tcd %>% filter(Step != "Classical HK") %>%
  filter(Boundary %in% c(.1, .12, .14, .16, .18, .2, .225, .25, .3)) %>%
  ggplot() +
  aes(y = Step, x = Difference, fill = Step, col = Step) +
  facet_grid(rows = vars(Conformity, N), cols = vars(Boundary), scales = "free_y") +
  geom_vline(xintercept = 0, linewidth = 2, color = "red") +
  geom_boxplot(alpha = 0.7) +
  guides(fill = "none", col = "none") +
  labs(x = "Difference from 'Classical HK'") +
  theme_light()

tcs = tc %>% #filter(Step != "Classical HK") %>%
  group_by(Conformity, N, Boundary, Step) %>%
  summarise(across(diversity:ESBG, list(mean = mean, sd = sd))) %>%
  ungroup()

tcds = tcd %>% filter(Step != "Classical HK") %>%
  group_by(Conformity, N, Boundary, Step) %>%
  summarise(across(Difference, list(mean = mean, sd = sd))) %>%
  ungroup()


tcds %>%
  ggplot() +
  aes(y = Difference_mean, x = Boundary, group = Step, col = Step) +
  facet_grid(cols = vars(N), rows = vars(Conformity), scales = "fixed") +
  geom_hline(yintercept = 0, linewidth = 1.2, color = "black") +
  geom_line(linewidth = 2, alpha = 0.4) +
  geom_point(alpha = 0.4, size = 5) +
  # guides(fill = "none", col = "none") +
  labs(x = "Difference from 'Classical HK'") +
  theme_light()

tcds %>%
  ggplot() +
  aes(y = Difference_sd, x = Boundary, group = Step, col = Step) +
  facet_grid(cols = vars(N), rows = vars(Conformity), scales = "fixed") +
  geom_hline(yintercept = 0, linewidth = 1.2, color = "black") +
  geom_line(linewidth = 2, alpha = 0.4) +
  geom_point(alpha = 0.4, size = 5) +
  # guides(fill = "none", col = "none") +
  labs(x = "SD of difference from 'Classical HK'") +
  theme_light()

tcs %>%
  ggplot() +
  aes(y = ESBG_mean, x = Boundary, group = Step, col = Step) +
  facet_grid(cols = vars(N), rows = vars(Conformity), scales = "fixed") +
  geom_hline(yintercept = 0, linewidth = 1.2, color = "black") +
  geom_line(linewidth = 2, alpha = 0.4) +
  geom_point(alpha = 0.4, size = 5) +
  # guides(fill = "none", col = "none") +
  labs(x = "Comparison of analysis/model steps: mean") +
  theme_light()

tcs %>%
  ggplot() +
  aes(y = ESBG_sd, x = Boundary, group = Step, col = Step) +
  facet_grid(cols = vars(N), rows = vars(Conformity), scales = "fixed") +
  geom_hline(yintercept = 0, linewidth = 1.2, color = "black") +
  geom_line(linewidth = 2, alpha = 0.4) +
  geom_point(alpha = 0.4, size = 5) +
  # guides(fill = "none", col = "none") +
  labs(x = "Comparison of analysis/model steps: SD") +
  theme_light()


filter(tc, Step == "Normal SPIRO, Normal Boundary", !(Boundary %in% c(0.1, 0.2, 0.3))) %>%
  count(Boundary, RS) %>% filter(n == 1536,  RS %in% c(10, 20, 30, 40, 50, 60)) %>% count(RS) %>%
  filter(n < 18) %>% arrange(desc(n), RS) %>% head(20)

