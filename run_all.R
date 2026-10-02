# Run the full simulation study

library(glmmTMB)
library(mice)
library(lme4)
library(pbapply)
library(ggplot2)

source("01_mice_impute_2l_glmmTMB.R")
source("02_simulation_study.R")
source("03_plots.R")
