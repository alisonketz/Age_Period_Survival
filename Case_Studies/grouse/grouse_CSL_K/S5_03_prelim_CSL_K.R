#######################################################
###
### Preliminary constants
###
#######################################################

n_brood <- max(d_fit$brood)

#######################################################
###
### Function for calculating kernel convolution
###
#######################################################

kernel_conv <- nimbleFunction(
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

#########################################################################################################################3
###
### Determining basis functions and knots 
### Using a function from the BCGAM R package
### Assuming a decreasing and convex shape for Age effects
###
#########################################################################################################################3

decconvex=function(x,t){
  n=length(x)
  k=length(t)-2
  m=k+3
  sigma=matrix(1:(m*n)*0,nrow=m,ncol=n)
  for(j in 1:k){
    i1=x<=t[j]
    sigma[j,i1] = x[i1]-t[1]
    i2=x>t[j]&x<=t[j+1]
    sigma[j,i2] = t[j]-t[1]+((t[j+1]-t[j])^3-(t[j+1]-x[i2])^3)/3/(t[j+1]-t[j])/(t[j+2]-t[j]) +(x[i2]-t[j])*(t[j+2]-t[j+1])/(t[j+2]-t[j])
    i3=x>t[j+1]&x<=t[j+2]
    sigma[j,i3] = t[j]-t[1] + (t[j+1]-t[j])^2/3/(t[j+2]-t[j]) + (t[j+2]-t[j+1])*(t[j+1]-t[j])/(t[j+2]-t[j]) +((t[j+2]-t[j+1])^3-(t[j+2]-x[i3])^3)/3/(t[j+2]-t[j+1])/(t[j+2]-t[j])
    i4=x>=t[j+2]
    sigma[j,i4] = t[j]-t[1] + (t[j+1]-t[j])^2/3/(t[j+2]-t[j]) + (t[j+2]-t[j+1])*(t[j+1]-t[j])/(t[j+2]-t[j]) +(t[j+2]-t[j+1])^2/3/(t[j+2]-t[j])
  }
  i1=x<=t[2]
  sigma[k+1,i1]=-(t[2]-x[i1])^3/3/(t[2]-t[1])^2
  i2=x>t[2]
  sigma[k+1,i2]=0
  i1=x<=t[k+1]
  sigma[k+2,i1]=x[i1]-t[1]
  i2=x>t[k+1]&x<=t[k+2]
  sigma[k+2,i2]=t[k+1]-t[1]+((t[k+2]-t[k+1])^2*(x[i2]-t[k+1])-(x[i2]-t[k+1])^3/3)/(t[k+2]-t[k+1])^2
  i3=x>t[k+2]
  sigma[k+2,i3]=t[k+1]-t[1]+((t[k+2]-t[k+1])^2*(t[k+2]-t[k+1])-(t[k+2]-t[k+1])^3/3)/(t[k+2]-t[k+1])^2
  sigma[k+3,]=x
  
  center.vector=apply(sigma,1,mean)
  
  list(sigma=-sigma, center.vector=-center.vector)
}

##############################################################
###
### Basis calculated from the BCGAM Meyer (2008) and 
### bcgam package
###
##############################################################

quant_age <- .2
knots_age <- c(1, round(quantile(d_fit$right.age - 1,
                                c(seq(quant_age,
                                      .99,
                                      by = quant_age),
                                  .99))))
knots_age_cgam <- unique(knots_age)
delta_i <- decconvex(1:nT_age, knots_age_cgam)
delta <- t(delta_i$sigma - delta_i$center.vector)
Z_age_cgam <- delta / max(delta)
nknots_age_cgam <- dim(Z_age_cgam)[2]

#######################################################
###
### Basis calculated from the spline package
###
#######################################################

quant_age <- .04
knots_age <- c(1,
               round(quantile(d_fit$right.age - 1,
                              c(seq(quant_age,
                                    .99,
                                    by = quant_age),
                                .99))))
knots_age_spline <- unique(knots_age)

splinebasis <- ns(1:nT_age, knots = knots_age_spline)

##A constraint matrix so age.effects = 0
constr_sumzero <- matrix(1, 1, nrow(splinebasis)) %*% splinebasis

##Get a basis for null space of constraint
qrc <- qr(t(constr_sumzero))
Z <- qr.Q(qrc,
          complete = TRUE)[, (nrow(constr_sumzero) + 1):ncol(constr_sumzero)]
Z_age_spline <- splinebasis %*% Z
nknots_age_spline <- dim(Z_age_spline)[2]

#############################################################
###
### plot of the basis functions
###
##############################################################

pdf("basis_function_age.pdf")
plot(1:nT_age,
     Z_age_spline[, 1],
     type = "l",
     main = "Basis Function Age Effect Spline")
for (i in 2:nknots_age_spline) {
  lines(1:nT_age, Z_age_spline[, i])
}
plot(1:nT_age,
     Z_age_cgam[, 1],
     type = "l",
     main = "Basis Function Age Effect CGAM")
for (i in 2:nknots_age_cgam) {
  lines(1:nT_age, Z_age_cgam[, i])
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
### Year 1 = 2015
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
### Year 2 = 2016
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

pdf("basis_function_period.pdf")
plot(1:nT_yr1,
    Z_period_yr1[, 1],
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
  lines(1:nT_yr2, Z_period_yr2[, i])
}
plot(1:nT_yr3,
     Z_period_yr3[, 1],
     type = "l",
     ylim = c(-1, 1),
     main = "Basis Function Period Effect Year 3")
for (i in 2:nknots_period_yr3) {
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
