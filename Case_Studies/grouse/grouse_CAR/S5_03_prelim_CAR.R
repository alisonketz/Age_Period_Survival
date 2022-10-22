#######################################################
###
### Preliminary constants
###
#######################################################

n_brood <- max(d_fit$brood)


######################################################################################################################3
###
### Parameters for CAR model for Age
###
######################################################################################################################3

#create num vector
num_age <- c(1, rep(2, nT_age - 2), 1)

#create adjacency vector along both years
temp <- as.matrix(bandSparse(n = nT_age, k = c(1), symmetric = TRUE))
temp2 <- matrix(0, nT_age, nT_age)
for (i in 1:nrow(temp2)) {
  temp2[i, ] <- temp[i, ] * c(1:nT_age)
}
adj_age <- t(temp2)[which(t(temp2) != 0)]

#weight vector into single vector
weights_age <- rep(1, length(adj_age))

#number of records in the adjacency vectors
nNage <- length(adj_age)


######################################################################################################################3
###
### Parameters for CAR model for Period Effects
###
######################################################################################################################3

###
### Period effects for different years
###

#year 1
num_period_yr1 <- c(1, rep(2, nT_yr1 - 2), 1)
temp <- as.matrix(bandSparse(n = nT_yr1, k = c(1), symmetric = TRUE))
temp2 <- matrix(0, nT_yr1, nT_yr1)
for (i in 1:nrow(temp2)) {
  temp2[i, ] <- temp[i, ] * c(1:nT_yr1)
}
adj_period_yr1 <- t(temp2)[which(t(temp2) != 0)]
weights_period_yr1 <- rep(1, length(adj_period_yr1))
nNperiod_yr1 <- length(adj_period_yr1)

#year 2
num_period_yr2 <- c(1, rep(2, nT_yr2 - 2), 1)
temp <- as.matrix(bandSparse( n = nT_yr2, k = c(1), symmetric = TRUE))
temp2 <- matrix(0, nT_yr2, nT_yr2)
for (i in 1:nrow(temp2)) {
  temp2[i, ] <- temp[i, ] * c(1:nT_yr2)
}
adj_period_yr2 <- t(temp2)[which(t(temp2) != 0)]
weights_period_yr2 <- rep(1,length(adj_period_yr2))
nNperiod_yr2 <- length(adj_period_yr2)

#year 3
num_period_yr3 <- c(1, rep(2, nT_yr3 - 2), 1)
temp <- as.matrix(bandSparse(n = nT_yr3, k = c(1), symmetric = TRUE))
temp2 <- matrix(0, nT_yr3, nT_yr3)
for (i in 1:nrow(temp2)) {
  temp2[i, ] <- temp[i, ] * c(1:nT_yr3)
}
adj_period_yr3 <- t(temp2)[which(t(temp2) != 0)]
weights_period_yr3 <- rep(1, length(adj_period_yr3))
nNperiod_yr3 <- length(adj_period_yr3)


######################################################
###
### number of MCMCr iterations, Chains, and Burn-in
###
######################################################

reps <- 50000
bin <- reps * .5
n_chains <- 3
n_thin <- 1
