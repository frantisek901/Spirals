#### Script for project with Mike, Ashley and Ashwin  II
##   (Adapted from my work for Project Group '#public-opinion-project' from SFI-CI Winter 2022
##    and the example from `nlrx` package help: ??"nlrx-package")
##
##   This script handles getting data from NetLogo model,
##   now named: groupFormation_pipeline.nlogo   
##   into R and storing them t/here for further network analyses.
##
##   NOTES:
##   A) Current idea of data structure is following:
##   Keep all needed information as the data files:
##   1) multi-layer weighted network data, where opinion distance is one layer and different types of relationships are other layers
##   2) dataframe/tibble with agents' traits and variables
##   This leads to flexible, but still lightweight datastructure allowing us in case of need store info on every tick of simmulation.
##
##   B) Current state of data stored is even better -- whole experiment is stored in one object/file!
##
##   C) Please consult help of `nlrx` package for further information and details here:
??"nlrx-package"
##

## Encoding: windows-1250
## Last edit: 2022-02-04 FranÈesko
##


## Head:
# Clear all
rm(list = ls())

# Packages
library(nlrx)
library(dplyr)
library(tibble)
library(igraph)
library(ggplot2)


#### Controling NetLogo from R with 'nlrx' package


# Creation of `nl` object: ---------------------------------------------

# Windows default NetLogo installation path (adjust to your needs!):
netlogopath = "c:/Program Files/NetLogo 6.2.2/"
## QUESTION to more experienced:
#  Does anyone know the way how to find the path to our own NetLogo instalations without setting it manually?
#  I am afraid that this manual setting forces many people avoid good things -- many ncluding me...

modelpath = "groupFormation_pipeline.nlogoo"

outpath = paste0(getwd(),"/OutputData")
# NOTE: Git doesn't cover this folder, do not forget create it manually in main project folder.

nl = nl(nlversion = "6.2.2",
        nlpath = netlogopath,
        modelpath = modelpath,
        jvmmem = 3072)  # NOTE: This reserves some memory, accomodate it according your computer.



# Attaching an experiment: ------------------------------------------------

nl@experiment = experiment(expname="firstTry",
                           outpath=outpath,
                           repetition=1,
                           tickmetrics="true",
                           idsetup="setup",
                           idgo="go",
                           runtime=50,
                           evalticks=0:50,
                           metrics = c("polarisation", "unweighted_polarisation"),  # NOTE: We might create here reporters which do not exist in the simulation!
                           metrics.turtles = list("turtles" = c("who", "Opinion-position", "P-speaking", "Speak?", "Uncertainty", "Tolerance", "Conformity", "Satisfied?", "group", "distance_to_centroid")),
                           metrics.links = list(
                             "comms" = c("[who] of end1", "[who] of end2")),
                           variables = list('number-of-agents' = list(min=20, max=100, step = 40),
                                            'agent-tolerance' = list(min=10, max=90, step = 40)),
                           constants = list("transparency" = 255))



# Attaching a simulation design: ------------------------------------------

nl@simdesign = simdesign_ff(nl=nl, nseeds=3)



# Producing results ----------------------------------

## Running a single simulation:
a = Sys.time()
results = run_nl_one(nl = nl, seed = getsim(nl, "simseeds")[1], siminputrow = 3)
Sys.time() - a
## Computational time on FranÈesko's computer: 14.6 sec.


## Running whole experiment:
a = Sys.time()
results = run_nl_all(nl = nl)
Sys.time() - a
## Computational time on FranÈesko's computer: 5.8 min. => 13 sec. per single simulation!



# Attaching results to `nl` object: ---------------------------------------

setsim(nl, "simoutput") = results



# Saving object 'nl' with whole experiment: -------------------------------

save(nl, file = "wholeExperiment.RData")



##  Please consult help of `nlrx` package for further information and details on properties of `nl` object:
??"nlrx-package"



# Create object with network data and visualize them --------

# One object for all network data pulled from 'nl':
nl.graph = nl_to_graph(nl)

# Tick 0 of the first simulation:
sim01.tick00 = nl.graph$spatial.links[[1]]
V(sim01.tick00)
E(sim01.tick00)
plot.igraph(sim01.tick00, vertex.size=8, vertex.label=NA, edge.arrow.size=0.2)

# Tick 0 of the third simulation:
sim03.tick00 = nl.graph$spatial.links[[103]]
V(sim03.tick00)
E(sim03.tick00)
plot.igraph(sim03.tick00, vertex.size=8, vertex.label=NA, edge.arrow.size=0.2)

# Tick 50 of the third simulation:
sim03.tick50 = nl.graph$spatial.links[[153]]
V(sim03.tick50)
E(sim03.tick50)
plot.igraph(sim03.tick50, vertex.size=8, vertex.label=NA, edge.arrow.size=0.2)



# Some graphs on agent variables: -----------------------------------------

df = nl@simdesign@simoutput %>% filter(`[step]` == 0) %>% select(7:9)
names(df) = c("Agents", "Negatives", "Links")

ggplot(df, aes(x = Agents, y = Negatives)) +
  geom_point() +
  theme_minimal()

ggplot(df, aes(x = Agents, y = Links, col = Negatives, size = Negatives)) +
  geom_point() +
  theme_minimal()

