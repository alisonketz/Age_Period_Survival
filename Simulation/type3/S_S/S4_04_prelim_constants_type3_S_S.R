########################################################
###
### Parameters for model, i.e. knots, consts
###
########################################################
########################################################
###
### Determining basis function and number of knots
### based on the spline package
###
########################################################

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

##Basis for age effect
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
     Z_age[, 1],
     type = "l",
     ylim = c(-1, 1),
     main = "Basis Function Age Effect")
for(i in 2:nknots_age) {
  lines(1:nT_age, Z_age[, i])
}
dev.off()


########################################
###
### Spline basis matrix for Period
###
##########################################

intvl_period <- 7
knots_period <- seq(1, nT_period, by = intvl_period)
knots_period <- unique(knots_period)

##Basis for period-effect
splinebasis <- bs(1:nT_period, knots = knots_period)

##A constraint matrix so period_effects = 0
constr_sumzero <- matrix(1, 1, nrow(splinebasis)) %*% splinebasis

##Get a basis for null space of constraint
qrc <- qr(t(constr_sumzero))
Z <- qr.Q(qrc,
          complete = TRUE
         )[, (nrow(constr_sumzero) + 1):ncol(constr_sumzero)]
Z_period <- splinebasis %*% Z
nknots_period <- dim(Z_period)[2]

pdf("basis_function_period_bs.pdf")
plot(1:nT_period,
     Z_period[, 1],
     type = "l",
     ylim = c(-1, 1),
     main = "Basis Function Period Effect")
for (i in 2:nknots_period) {
  lines(1:nT_period, Z_period[, i])
}
dev.off()

########################################################
###
### Number of MCMC iterations, Chains, Burn-in, Thining
###
########################################################

reps <- 50000
bin <- reps * .5
n_chains <- 3
n_thin <- 1