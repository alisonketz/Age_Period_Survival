########################################################
###
### Parameters for model, i.e. knots, consts
###
########################################################

########################################################
###
### Function for kernel convolution 
###
########################################################

kernel_conv <- nimbleFunction(
  run = function(nT = double(0),
                 Z = double(2),
                 stauk = double(0),
                 nconst = double(0),
                 tauk = double(0),
                 nknots = double(0),
                 alphau = double(1)
  ) {
    temp <- nimMatrix(value = 0,nrow = nT, ncol = nknots)
    temp1 <- nimMatrix(value = 0,nrow = nT, ncol = nknots)
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
        temp[i,j] <- (temp1[i, j] / temp2[j]) * alphau[j]
      }
      KA[i] <- sum(temp[i, 1:nknots])
    }

    muKA <- mean(KA[1:nT])
    KA[1:nT] <- KA[1:nT] - muKA

    returnType(double(1))
    return(KA[1:nT])
  })

ckernel_conv <- compileNimble(kernel_conv)

######################################################################
###
### Determining basis function and number of knots
### based on the spline package
###
######################################################################

quant_age <- .05
knots_age <- c(1, round(quantile(right_age - 1,
                                 c(seq(quant_age, .99, by = quant_age),
                                 .99)
                                )
                        )
              )
knots_age <- unique(knots_age)
if (max(knots_age) == nT_age) {
    knots_age[nknots_age] <- knots_age[nknots_age] - 1
  }

##Basis for age effects
splinebasis <- ns(1:nT_age, knots = knots_age)

##A constraint matrix so age_effects = 0
constr_sumzero <- matrix(1, 1, nrow(splinebasis)) %*% splinebasis

##QR factorization to obtain null space of constraint
qrc <- qr(t(constr_sumzero))
Z <- qr.Q(qrc,
          complete = TRUE
          )[, (nrow(constr_sumzero) + 1):ncol(constr_sumzero)]
Z_age <- splinebasis %*% Z

#number of knows
nknots_age <- dim(Z_age)[2]
nknots_age

#plotting the basis function
pdf("basis_function_age_ns.pdf")
plot(1:nT_age,
     Z_age[,1],
     type = "l",
     ylim = c(-1, 1),
     main = "Basis Function Age Effect")
for(i in 2:nknots_age) {
  lines(1:nT_age, Z_age[, i])
}
dev.off()


#############################################################
###
### Period effects Basis Function - Kernel Convolution
###
#############################################################

intvl_period <- 1
knots_period <- seq(1, nT_period, by = intvl_period)
knots_period <- unique(knots_period)
nknots_period <- length(knots_period)

Z_period <- matrix(0, nT_period, nknots_period)
for (i in 1:nrow(Z_period)) {
  for (j in 1:nknots_period) {
    Z_period[i, j] <- abs(i - knots_period[j])
  }
}

########################################################
###
### Number of MCMC iterations, Chains, Burn-in, Thining
###
########################################################

reps <- 50000
bin <- reps * .5
n_chains <- 3
n_thin <- 1