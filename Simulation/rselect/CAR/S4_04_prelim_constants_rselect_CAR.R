######################################################################################################################3
###
### Parameters for CAR model 
###
######################################################################################################################3

#create num vector
num_period <- c(1, rep(2, nT_period - 2), 1)

#create adjacency vector along both years
temp <- as.matrix(bandSparse(n = nT_period, k = c(1),symmetric = T))
temp2 <- matrix(0, nT_period, nT_period)
for (i in 1:nrow(temp2)) {
  temp2[i, ] <- temp[i, ] * c(1:nT_period)
}
adj_period <- t(temp2)[which(t(temp2) != 0)]

#weight vector into single vector
weights_period <- rep(1, length(adj_period))

#length of the study period in total
nperiod <- nT_period

#number of records in the adjacency vectors
nNperiod <- length(adj_period)

############################################################################
###
### Parameters for CAR model for period
###
############################################################################

#create num vector
num_age <- c(1, rep(2, nT_age - 2), 1)

#create adjacency vector along both years
temp <- as.matrix(bandSparse(n = nT_age, k = c(1), symmetric = T))
temp2 <- matrix(0, nT_age, nT_age)
for (i in 1:nrow(temp2)) {
  temp2[i, ] <- temp[i, ] * c(1:nT_age)
}
adj_age <- t(temp2)[which(t(temp2) != 0)]

#weight vector into single vector
weights_age <- rep(1, length(adj_age))

#length of the study age in total
nage <- nT_age

#number of records in the adjacency vectors
nNage <- length(adj_age)


########################################################
###
### Number of MCMC iterations, Chains, Burn-in, Thining
###
########################################################

#number of MCMCr iterations, Chains, and Burn-in
reps <- 50000
bin <- reps * .5
n_chains <- 3
n_thin <- 1