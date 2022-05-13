#######################################################
###
### Preliminary constants
###
#######################################################

n_brood <- max(d_fit$brood)

#########################################################################################################################3
###
### Function for calculating kernel convolution
###
#########################################################################################################################3

kernel_conv<- nimbleFunction(
  run = function(nT = double(0),
                 Z = double(2),
                 stauk = double(0),
                 nconst = double(0),
                 tauk = double(0),
                 nknots = double(0),
                 alphau = double(1)
  ){
    temp <- nimMatrix(value = 0, nrow = nT, ncol = nknots)
    temp1 <- nimMatrix(value = 0, nrow = nT, ncol = nknots)
    temp2 <- nimNumeric(nknots)
    KA <- nimNumeric(nT)

    for (i in 1:nT) {
      for (j in 1:nknots) {
        temp1[i, j] <- stauk * nconst * exp(-0.5 * Z[i, j]^2 * tauk)
      }
    }

    for (j in 1:nknots) {
      temp2[j] <- sum(temp1[1:nT, j])
    }

    for (i in 1:nT) {
      for (j in 1:nknots) {
        temp[i, j] <- (temp1[i, j] / temp2[j]) * alphau[j]
      }
      KA[i] <- sum(temp[i, 1:nknots])
    }
    muKA <- mean(KA[1:nT])
    KA[1:nT] <- KA[1:nT] - muKA

    returnType(double(1))
    return(KA[1:nT])
  })

Ckernel_conv <- compileNimble(kernel_conv)


#######################################################
###
### Basis calculated from the spline package
###
#######################################################

intvl_age <- 2
knots_age <- seq(1,nT_age, by = intvl_age)
knots_age <- unique(knots_age)
nknots_age <- length(knots_age)

Z_age <- matrix(0,nT_age, nknots_age)
for (i in 1:nrow(Z_age)) {
  for (j in 1:nknots_age) {
    Z_age[i, j] <- abs(i - knots_age[j]) #absolute value
  }
}

#############################################################
###
### plot of the basis functions
###
##############################################################
pdf("basis_function_age.pdf")
plot(1:nT_age,
     Z_age[, 1],
     type = "l",
     main = "Basis Function Age Effect")
for (i in 2:nknots_age) {
  lines(1:nT_age, Z_age[, i])
}
dev.off()

#############################################
###
### Spline basis matrix for Period Effects
###
#############################################

intvl_period <- 1 # interval for knots

######################
###
### Year 1 = 2014
###
######################

intvl_period_yr1 <- intvl_period # interval for knots
knots_period_yr1 <- unique(c(seq(1,
                                 nT_yr1,
                                 intvl_period_yr1),
                             nT_yr1)) #knots
nknots_period_yr1 <- length(knots_period_yr1)

Z_period_yr1 <- matrix(0,
                     nT_yr1,
                     nknots_period_yr1)
for (i in 1:nrow(Z_period_yr1)) {
  for (j in 1:nknots_period_yr1) {
    Z_period_yr1[i, j] <- abs(i - knots_period_yr1[j]) #absolute value
  }
}

######################
###
### Year 2 = 2015
###
######################

intvl_period_yr2 <- intvl_period # interval for knots

knots_period_yr2 <- unique(c(seq(1,
                               nT_yr2,
                               intvl_period_yr2),
                            nT_yr2)) #knots
nknots_period_yr2 <- length(knots_period_yr2)
Z_period_yr2 <- matrix(0,
                       nT_yr2,
                       nknots_period_yr2)
for (i in 1:nrow(Z_period_yr2)) {
  for (j in 1:nknots_period_yr2) {
    Z_period_yr2[i, j] <- abs(i - knots_period_yr2[j]) #absolute value
  }
}
######################
###
### Year 3 = 2017
###
######################

intvl_period_yr3 <- intvl_period # interval for knots
knots_period_yr3 <- unique(c(seq(1,
                                 nT_yr3,
                                 intvl_period_yr3),
                             nT_yr3)) #knots
nknots_period_yr3 <- length(knots_period_yr3)

Z_period_yr3 <- matrix(0, nT_yr3, nknots_period_yr3)
for (i in 1:nrow(Z_period_yr3)) {
  for (j in 1:nknots_period_yr3) {
    Z_period_yr3[i, j] <- abs(i - knots_period_yr3[j])
  }
}

##############################
###
### Plots of basis functions
###
##############################
pdf("basis_function_time_bs.pdf")
plot(1:nT_yr1,
    Z_period_yr1[,1],
    type = "l",
    ylim = c(-1, 1),
    main = "Basis Function Period Effects Year 1")
for (i in 2:nknots_period_yr1) {
  lines(1:nT_yr1,
        Z_period_yr1[, i])
}
plot(1:nT_yr2,
     Z_period_yr2[, 1],
     type = "l",
     ylim = c(-1, 1),
     main = "Basis Function Period Effect Year 2")
for (i in 2:nknots_period_yr2) {
  lines(1:nT_yr2,Z_period_yr2[, i])
}
plot(1:nT_yr3,
     Z_period_yr3[,1],
     type = "l",
     ylim = c(-1, 1),
     main = "Basis Function Period Effect Year 3")
for(i in 2:nknots_period_yr3){
  lines(1:nT_yr3,
        Z_period_yr3[, i])
}
dev.off()

######################################################
###
### number of MCMCr iterations, Chains, and Burn-in
###
######################################################

reps <- 50000
bin <- reps * .5
n_chains <- 3
n_thin <- 1
