#########################################################################################################################3
###
### Determining basis function and number of knots 
###
#########################################################################################################################3

#########################################################################################################################3
###
### function for kernel convolution model
###
#########################################################################################################################3

#function to calculate cumulative survival probability for likelihood
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
        temp1[i, j] <- stauk * nconst * exp(-0.5 * Z[i,j]^2 * tauk)
      }
    }
    
    for(j in 1:nknots){
      temp2[j] <- sum(temp1[1:nT, j])
    }
    
    for(i in 1:nT){
      for(j in 1:nknots){
        temp[i, j] <- (temp1[i, j]/temp2[j]) * alphau[j]
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
### Determining basis function and number of knots 
###
#########################################################################################################################3
#########################################################################################################################3

# Decreasing and Convex
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
### Basis calculated from the BCGAM Meyer (2008) and bcgam package
###
##############################################################

quant_age <- .25
knots_age_cgam <- c(1,round(quantile(right_age-1,c(seq(quant_age,.99,by=quant_age),.99))))
knots_age_cgam <- unique(knots_age_cgam)
delta_i <- decconvex(1:nT_age, knots_age_cgam)
delta <- t(delta_i$sigma - delta_i$center.vector)
Z_age_cgam  <- delta / max(delta)
nknots_age_cgam <- dim(Z_age_cgam)[2]


###############################
###
### basis functions
### using the spline package
###
###############################

quant_age <- .05
knots_age_spline <- c(1,round(quantile(df_fit$right_age-1,c(seq(quant_age,.99, by=quant_age),.99))))
knots_age_spline <- unique(knots_age_spline)
nknots_age_spline <- length(knots_age_spline)

##Basis for age hazard
splinebasis <- ns(1:nT_age, knots = knots_age_spline)#,intercept=TRUE,

##A constraint matrix so time.effects = 0
constr_sumzero <- matrix(1, 1, nrow(splinebasis)) %*% splinebasis

##Get a basis for null space of constraint
qrc <- qr(t(constr_sumzero))
Z <- qr.Q(qrc, complete = TRUE)[, (nrow(constr_sumzero) + 1):ncol(constr_sumzero)]
Z_age_spline <- splinebasis%*%Z
nknots_age_spline<- dim(Z_age_spline)[2]


#############################################################
###
### plot of the basis functions
###
##############################################################
pdf("basis_function_age_pdf")
plot(1:nT_age,
     Z_age_spline[, 1],
     ylim = c(-1, 1),
     type = "l",
     main = "Basis Function Age Effect Spline")
for (i in 2:nknots_age_spline) {
  lines(1:nT_age, Z_age_spline[, i])
}
plot(1:nT_age,
     Z_age_cgam[, 1],
     ylim = c(-1, 1),
     type = "l",
     main = "Basis Function Age Effect CGAM")
for (i in 2:nknots_age_cgam) {
  lines(1:nT_age, Z_age_cgam[, i])
}
dev.off()

#########################################################################################################################
###
### Time - kernel convolution
###
#########################################################################################################################

intvl_period <- 1
knots_period <- seq(1,
                    nT_period,
                    by = intvl_period)
knots_period <- unique(knots_period)
nknots_period <- length(knots_period)

Z_period <- matrix(0, nT_period, nknots_period)
for (i in 1:nrow(Z_period)) {
  for (j in 1:nknots_period) {
    Z_period[i, j] <- abs(i - knots_period[j]) #absolute value
  }
}



#number of MCMCr iterations, Chains, and Burn-in
reps<-50000
bin<-reps*.5
n_chains<-3
n_thin<-1