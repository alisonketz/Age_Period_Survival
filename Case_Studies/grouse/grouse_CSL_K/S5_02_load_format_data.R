###########################################################################################
###
### Load data for survival models
###
###########################################################################################

d_fit <- read.csv("../grouse_data.csv")
n_fit <- dim(d_fit)[1]

#vector to calibrate indexes
age2date <- d_fit$left.time - d_fit$left.age

### Number of Age effects 
nT_age <- max(d_fit$right.age) - 1

#Number of Period effects for each 
#year seperately and combined
nT_yr1 <- length(as.numeric(names(table(d_fit[d_fit[,5]==1,1]))[1]):as.numeric(names(table(d_fit[d_fit[,5]==1,2]))[length(table(d_fit[d_fit[,5]==1,2]))]))
nT_yr2 <- length(as.numeric(names(table(d_fit[d_fit[,6]==1,1]))[1]):as.numeric(names(table(d_fit[d_fit[,6]==1,2]))[length(table(d_fit[d_fit[,6]==1,2]))]))
nT_yr3 <- length(as.numeric(names(table(d_fit[d_fit[,7]==1,1]))[1]):as.numeric(names(table(d_fit[d_fit[,7]==1,2]))[length(table(d_fit[d_fit[,7]==1,2]))]))
nT_period_total <- nT_yr1+nT_yr2+nT_yr3
nT_period_max <- max(nT_yr1, nT_yr2, nT_yr3)
nT_period <- max(d_fit$right.time)

