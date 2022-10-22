##########################################################################
###
### Age-period survival of White-tailed deer
###
### Alison Ketz
### 01/01/2021
###
##########################################################################

rm(list = ls())

setwd("~/Documents/Survival/surv_app/Case_Studies/wtd/wtd_CSL_K")

library(nimble)
library(Matrix)
library(coda)
library(lattice)
library(splines)
library(Hmisc)
library(lubridate)
library(ggplot2)
library(gridExtra)
library(xtable)
library(RColorBrewer)

####################################
### Load data
####################################

source("S5_02_load_format_data.R")


##################################
### Setting constants
##################################

source("S5_03_prelim_CSL_K.R")


##################################
### Run model
##################################

source("S5_04_run_model_CSL_K.R")

######################################
### Post processing: summaries, plots
######################################

source("S5_05_results_plots_sum_CSL_K.R")

