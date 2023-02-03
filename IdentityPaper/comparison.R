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
              Step = "No Identity, Constant Boundary",
              Identity_Type = "none", SPIRO_Distribution = "none",
              Boundary_STD = 0, Conformity_STD = 0, RS = 0,
              SPIRO_STD = factor(0),SPIRO_Mean = factor(0.25))) %>%

  # Now we do changes on joined file:
  mutate(
    Step = factor(Step) %>% fct_inorder() %>% fct_rev()) %>%

  # Selecting out obsolete variables:
  select(-(17:31)) %>%

  # Filtering only for comparable results:
  filter(Conformity %in% c(.2, .8), Conformity_STD %in% c(0, .1),
         as.numeric(as.character(SPIRO_Mean)) %in% seq(0.25, 0.85, 0.12),
         Boundary_STD %in% seq(0, 0.15, .05), Boundary >= 0.1, Boundary <= 0.3) %>%
  filter(Boundary %in% c(.1, .2, .3))



tc %>% count(HK_distribution, Use_identity, Identity_Type, Step)
tc %>% count(fct_rev(Step), Boundary)




