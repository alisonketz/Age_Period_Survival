#######################################################
###
### Preliminary constants
###
#######################################################

n_brood <- max(d_fit$brood)

#######################################################
###
### Basis calculated from the spline package
###
#######################################################
quant_age <- .05
knots_age <- c(1,round(quantile(d_fit$right.age-1,c(seq(quant_age,.99, by=quant_age),.99))))
knots_age <- unique(knots_age)

splinebasis <- ns(1:nT_age, knots = knots_age)

##A constraint matrix so age.effects = 0
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

intvl_period <- 7 # interval for knots

######################
###
### Year 1 = 2015
###
######################

intvl_period_yr1<-intvl_period # interval for knots
knots_period_yr1<-unique(c(seq(1,nT_yr1,intvl_period_yr1),nT_yr1)) #knots

##Basis for time-varying cause
splinebasis <- bs(1:nT_yr1,knots=knots_period_yr1)#

##A constraint matrix so time.effects = 0
constr_sumzero <- matrix(1, 1, nrow(splinebasis)) %*% splinebasis

##Get a basis for null space of constraint
qrc <- qr(t(constr_sumzero))
Z <- qr.Q(qrc,complete=TRUE)[,(nrow(constr_sumzero)+1):ncol(constr_sumzero)]
Z_period_yr1 <- splinebasis%*%Z
nknots_period_yr1<- dim(Z_period_yr1)[2]
nknots_period_yr1
######################
###
### Year 2 = 2016
###
######################

intvl_period_yr2 <- intvl_period # interval for knots
knots_period_yr2 <- unique(c(seq(1, nT_yr2, intvl_period_yr2), nT_yr2)) #knots

##Basis for time-varying cause
splinebasis <- bs(1:nT_yr2, knots = knots_period_yr1)#

##A constraint matrix so period.effects = 0
constr_sumzero <- matrix(1, 1, nrow(splinebasis)) %*% splinebasis

##Get a basis for null space of constraint
qrc <- qr(t(constr_sumzero))
Z <- qr.Q(qrc, complete = TRUE)[, (nrow(constr_sumzero) + 1):ncol(constr_sumzero)]
Z_period_yr2 <- splinebasis%*%Z
nknots_period_yr2 <- dim(Z_period_yr2)[2]

######################
###
### Year 3 = 2017
###
######################

intvl_period_yr3 <- intvl_period # interval for knots
knots_period_yr3 <- unique(c(seq(1, nT_yr3, intvl_period_yr3), nT_yr3)) #knots

##Basis for time-varying cause
splinebasis <- bs(1:nT_yr3, knots = knots_period_yr3)#

##A constraint matrix so time.effects = 0
constr_sumzero <- matrix(1, 1, nrow(splinebasis)) %*% splinebasis

##Get a basis for null space of constraint
qrc <- qr(t(constr_sumzero))
Z <- qr.Q(qrc,complete=TRUE)[, (nrow(constr_sumzero) + 1):ncol(constr_sumzero)]
Z_period_yr3 <- splinebasis %*% Z
nknots_period_yr3 <- dim(Z_period_yr3)[2]

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
