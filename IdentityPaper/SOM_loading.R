#### Script just for reading data in for SOM

## Encoding: windows-1250
## Created: 2023-02-09 FrK
## Edited: 2023-02-13 FrK


## NOTES:
##
##



# Head --------------------------------------------------------------------

# Clearing all
rm(list=ls())

# Packages:
library(tidyverse)
library(readr)
library(dplyr)
library(tibble)
library(forcats)
library(ggplot2)
library(rstatix)
library(stringr)
library(knitr)



# Reading data in: --------------------------------------------------------

# STEP 1 ------------------------------------------------------------------

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


# For Joining script we need these data stored under different name:
ts10 = tb %>% add_row(filter(ted, file == "random position") %>% select(-file)) %>%
  filter(N %in% 100:101)

# For detective purposes we read in more smooth and detailed data:
ts10s = read_csv("ClassicalHK_smooth-table.csv", skip = 6) %>%

  # Selecting and renaming...
  select(HK_distribution = 4, 2, Use_identity = 12, 13:16,
  N = 3, Boundary = 7, 8, Conformity = 10, 11,
  39, diversity = 40, extremness = 41, ESBG = 42,
  boundary_mean = 43, boundary_sd = 44,
  conformity_mean = 45, conformity_sd = 46, 47:54) %>%

  mutate(
    Step = factor("Classical HK",
                  levels = c("Classical HK", "No Identity, Constant Boundary",
                             "No Identity, Normal Boundary",
                             "Constant SPIRO, Normal Boundary", "Normal SPIRO, Normal Boundary")),
    across(c(HK_distribution, Use_identity, SPIRO_Mean, SPIRO_STD), ~factor(.x)))



# STEP 2 ------------------------------------------------------------------

# Creating object 'raw' (tibble): Loading....
raw2 = read_csv("ClassicalHK_heterogenousParameters_RS01-05.csv", skip = 6) %>%
  add_row(read_csv("ClassicalHK_heterogenousParameters_RS06-10.csv", skip = 6))
for (i in seq(11, 56, 5)) {
  raw2 = raw2 %>%
    add_row(read_csv(paste0("ClassicalHK_heterogenousParameters_RS", i, "-", i + 4, ".csv"), skip = 6))
}

# Transforming 'raw' to clean 'ts'
ts20 = raw2 %>%
  # Selecting and renaming...
  select(HK_distribution = 4, Present_opinion = 5,
         2, Use_identity = 12,
         N = 3, Boundary = 7, 8, Conformity = 10, 11,
         39, diversity = 40, extremness = 41, ESBG = 42,
         boundary_mean = 43, boundary_sd = 44,
         conformity_mean = 45, conformity_sd = 46) %>%

  # Processing.
  mutate(
    across(.cols = c(1:2, 4), factor),
    even_N = ((N %%2) == 0))



# STEP 3 ------------------------------------------------------------------

# Creating object 'raw' (tibble): Loading....
raw3 = read_csv("GlobalIdentityHK_heterogenousParameters_RS01-05-table.csv", skip = 6) %>%
  add_row(read_csv("GlobalIdentityHK_heterogenousParameters_RS06-10-table.csv", skip = 6)) %>%
  add_row(read_csv("GlobalIdentityHK_heterogenousParameters_RS20.csv", skip = 6)) %>%
  add_row(read_csv("GlobalIdentityHK_heterogenousParameters_RS55.csv", skip = 6))
for (i in seq(11, 56, 5)) {
  raw3 = raw3 %>%
    add_row(read_csv(paste0("GlobalIdentityHK_heterogenousParameters_RS", i, "-", i + 4, "-table.csv"), skip = 6))
}

# Transforming 'raw3' to clean 't3'
ts30 = raw3 %>%
  # Selecting and renaming...
  select(HK_distribution = 4, 2, Use_identity = 12, 13:15,
         N = 3, Boundary = 7, 8, Conformity = 10, 11,
         39, diversity = 40, extremness = 41, ESBG = 42,
         boundary_mean = 43, boundary_sd = 44,
         conformity_mean = 45, conformity_sd = 46) %>%
  # Dropping duplicate observations due to experiment runs dubling:
  distinct()



# STEP 4.1 ----------------------------------------------------------------

# Creating object 'raw' (tibble): Loading....
raw4 = read_csv("Step4_indID-hetPar_RS01-05.csv", skip = 6) %>%
  add_row(read_csv("Step4_indID-hetPar_RS06-10.csv", skip = 6))
for (i in c(seq(30, 55, 5), 39, 49, "35B", "45B")) {
  raw4 = raw4 %>%
    add_row(read_csv(paste0("Step4_indID-hetPar_RS", i, ".csv"), skip = 6))
}
for (i in seq(11, 56, 5)) {
  raw4 = raw4 %>%
    add_row(read_csv(paste0("Step4_indID-hetPar_RS", i, "-", i + 4, ".csv"), skip = 6))
}


# Creating object 'raw' (tibble): Loading....
spiro = read_csv("SPIROdistribution_indID-hetPar_RS01-10.csv", skip = 6) %>%
  add_row(read_csv("SPIROdistribution_indID-hetPar_RS20.csv", skip = 6))
for (i in seq(11, 51, 10)) {
  spiro = spiro %>%
    add_row(read_csv(paste0("SPIROdistribution_indID-hetPar_RS", i, "-", i + 9, ".csv"), skip = 6))
}
spiro = select(spiro, -1) %>% distinct()


# Preparing additional data:
raw41 = read_csv("Step4.1_b011_indID-hetPar_RS18-20-table.csv", skip = 6)
for (i in c(13, 16, 19, 23, 26, 29)) {
  raw41 = read_csv(paste0("Step4.1_b0", i, "_indID-hetPar_RS01-10-table.csv"), skip = 6) %>%
    add_row(raw41)
  if (i %in% c(13, 16, 19)) {
    raw41 = read_csv(paste0("Step4.1_b0", i, "_indID-hetPar_RS01-10_B-table.csv"), skip = 6) %>%
      add_row(raw41)
  }
}
for (i in c(13, 16, 19, 23, 26, 29)) {
  for (j in seq(11, 51, 10)) {
    raw41 = read_csv(paste0("Step4.1_b0", i, "_indID-hetPar_RS", j,"-", j + 9, "-table.csv"), skip = 6) %>%
      add_row(raw41)
  }
}
for (i in c(13, 16, 19, 23, 26, 29)) {
  for (j in seq(11, 51, 10)) {
    raw41 = read_csv(paste0("Step4.1_b0", i, "_indID-hetPar_RS", j,"-", j + 9, "-table.csv"), skip = 6) %>%
      add_row(raw41)
    f = paste0("Step4.1_b0", i, "_indID-hetPar_RS", j,"-", j + 9, "_B-table.csv")
    if (file.exists(f)) raw41 = read_csv(f, skip = 6) %>% add_row(raw41)
  }
}

# Transforming 'raw4' to clean 'ts'
ts41 = left_join(raw4, spiro, by = names(raw4)[c(2:32, 34:37)]) %>%

  # Adding extra data from raw41:
  add_row(select(raw41, -`max-ticks`, -`[step]`)) %>%

  # Selecting and renaming...
  select(HK_distribution = 4, 2, Use_identity = 12, 13:16,
         N = 3, Boundary = 7, 8, Conformity = 10, 11,
         39, diversity = 40, extremness = 41, ESBG = 42,
         boundary_mean = 43, boundary_sd = 44,
         conformity_mean = 45, conformity_sd = 46, 49:56) %>%
  # Dropping duplicate observations due to experiment runs dubling:
  drop_na() %>% distinct() %>%

  # Filtering out non-needed simulations using non-needed values of SPIRO_Mean:
  filter(as.numeric(as.character(SPIRO_Mean)) %in% seq(0.25, 0.85, 0.12)) %>%

  # Processing.
  mutate(
    even_N = ((N %%2) == 0),
    across(.cols = c(1, 3:7), factor),
    spiro_mean = ((SPIRO_0.15 * .15) + (SPIRO_0.25 * .25) + (SPIRO_0.35 * .35) + (SPIRO_0.45 * .45) + (SPIRO_0.55 * .55) + (SPIRO_0.65 * .65) + (SPIRO_0.75 * .75) + (SPIRO_0.85 * .85)) / N,
    spiro_sd = sqrt((SPIRO_0.15 * ((spiro_mean - 0.15) ^ 2) + SPIRO_0.25 * ((spiro_mean - 0.25) ^ 2) + SPIRO_0.35 * ((spiro_mean - 0.35) ^ 2) + SPIRO_0.45 * ((spiro_mean - 0.45) ^ 2) + SPIRO_0.55 * ((spiro_mean - 0.55) ^ 2) + SPIRO_0.65 * ((spiro_mean - 0.65) ^ 2) + SPIRO_0.75 * ((spiro_mean - 0.75) ^ 2) + SPIRO_0.85 * ((spiro_mean - 0.85) ^ 2)) / (N - 1)) )



# Combining all together: -------------------------------------------------

# We start by joining files -- we do just changes needed for joining.
tc = ts41 %>%
  mutate(Step = "Normal SPIRO, Normal Boundary") %>%
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
              Boundary_STD = 0, Conformity_STD = 0, RS = if_else(RS < 61, RS, 0),
              SPIRO_STD = factor(0),SPIRO_Mean = factor(0.25))) %>%

  # Now we do changes on joined file:
  mutate(
    Step = factor(Step) %>% fct_inorder() %>% fct_rev() %>%  fct_relevel("Classical HK")) %>%

  # Selecting out obsolete variables:
  select(-(17:31)) %>%

  # Filtering only for comparable results:
  filter(Conformity %in% c(.2, .8), Conformity_STD %in% c(0, .1),
         as.numeric(as.character(SPIRO_Mean)) %in% seq(0.25, 0.85, 0.12),
         Boundary_STD %in% seq(0, 0.15, .05), Boundary >= 0.1, Boundary <= 0.3)



# Saving data -------------------------------------------------------------
save(tc, ts10s, file = "SOM.RData")

