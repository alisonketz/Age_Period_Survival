##########################################################################################################
###
### Survival WTD - CSL_S
### Age-period survival
###
### Alison Ketz
### 10/08/2020
###
##########################################################################################################

rm(list = ls())

setwd("~/Documents/Survival/surv_app/Case_Studies/wtd/wtd_CSL_S")

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

source("S5_03_prelim_CSL_S.R")

##################################
### Run model
##################################

source("S5_04_run_model_CSL_S.R")

######################################
### Post processing: summaries, plots
######################################

source("S5_05_results_plots_sum_CSL_S.R")

