###########################################################################################
###
### Load data for survival models
###
###########################################################################################

d_fit <- read.csv("../wtd_data.csv")
n_fit <- dim(d_fit)[1]

#vector to calibrate indexes
age2date <- d_fit$left.time - d_fit$ageCapture

### Number of Age effects and Number of Period effects
nT_age <- max(d_fit$ageRight) - 1
nT_period <- max(d_fit$right.time) - 1

