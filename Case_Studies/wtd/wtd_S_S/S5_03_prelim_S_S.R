#######################################################
###
### Preliminary constants
###
#######################################################
#######################################################
###
### Basis calculated from the spline package
###
#######################################################
quant_age <- .05
knots_age <- c(1,round(quantile(d_fit$ageRight-1,c(seq(quant_age,.99, by=quant_age),.99))))
knots_age <- unique(knots_age)

splinebasis <- ns(1:nT_age, knots = knots_age)

##A constraint matrix so period.effects = 0
constr_sumzero <- matrix(1, 1, nrow(splinebasis)) %*% splinebasis

##Get a basis for null space of constraint
qrc <- qr(t(constr_sumzero))
Z <- qr.Q(qrc,
          complete = TRUE)[, (nrow(constr_sumzero) + 1):ncol(constr_sumzero)]
Z_age <- splinebasis %*% Z
nknots_age <- dim(Z_age)[2]

#############################################################
###
### plot of the basis functions
###
##############################################################
pdf("basis_function_age.pdf")
plot(1:nT_age,
     Z_age[, 1],
     ylim = c(-1, 1),
     type = "l",
     main = "Basis Function Age Effect")
for (i in 2:nknots_age) {
  lines(1:nT_age, Z_age[, i])
}
dev.off()

########################################
###
### Spline basis matrix for Period
###
##########################################

intvl_period <- 12
knots_period <- c(1,
               seq(5,
                   nT_period - 1,
                   by = intvl_period))
knots_period <- unique(knots_period)
splinebasis <- bs(1:nT_period, knots = knots_period)

##A constraint matrix so period_effects = 0
constr_sumzero <- matrix(1, 1, nrow(splinebasis)) %*% splinebasis

##Get a basis for null space of constraint
qrc <- qr(t(constr_sumzero))
Z <- qr.Q(qrc,
          complete = TRUE)[, (nrow(constr_sumzero) + 1):ncol(constr_sumzero)]
Z_period <- splinebasis %*% Z
nknots_period<- dim(Z_period)[2]


######################################################
###
### number of MCMCr iterations, Chains, and Burn-in
###
######################################################

reps <- 50000
bin <- reps * .5
n_chains <- 3
n_thin <- 1
