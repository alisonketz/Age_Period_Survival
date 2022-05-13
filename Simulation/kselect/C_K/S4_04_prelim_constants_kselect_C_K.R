########################################################
###
### Parameters for model, i.e. knots, consts
###
########################################################

########################################################
###
### Function for kernel convolution model
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
  ){
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
### Function to obtain constrained GAM basis function
### Basis calculated from the BCGAM Meyer (2008) and bcgam package
### Convex shape neither decreasing or increasing
###
#######################################################################

convex = function(x, t, pred.new=TRUE){
  n=length(x)
  k=length(t)-2
  m=k+2
  sigma=matrix(1:m*n,nrow=m,ncol=n)
  for(j in 1:(k-1)){
    i1=x<=t[j]
    sigma[j,i1] = 0
    i2=x>t[j]&x<=t[j+1]
    sigma[j,i2] = (x[i2]-t[j])^3 / (t[j+2]-t[j]) / (t[j+1]-t[j])/3
    i3=x>t[j+1]&x<=t[j+2]
    sigma[j,i3] = x[i3]-t[j+1]-(x[i3]-t[j+2])^3/(t[j+2]-t[j])/(t[j+2]-t[j+1])/3+(t[j+1]-t[j])^2/3/(t[j+2]-t[j])-(t[j+2]-t[j+1])^2/3/(t[j+2]-t[j])
    i4=x>t[j+2]
    sigma[j,i4]=(x[i4]-t[j+1])+(t[j+1]-t[j])^2/3/(t[j+2]-t[j])-(t[j+2]-t[j+1])^2/3/(t[j+2]-t[j])
  }
  i1=x<=t[k]
  sigma[k,i1] = 0
  i2=x>t[k]&x<=t[k+1]
  sigma[k,i2] = (x[i2]-t[k])^3 / (t[k+2]-t[k]) / (t[k+1]-t[k])/3
  i3=x>t[k+1]
  sigma[k,i3] = x[i3]-t[k+1]-(x[i3]-t[k+2])^3/(t[k+2]-t[k])/(t[k+2]-t[k+1])/3+(t[k+1]-t[k])^2/3/(t[k+2]-t[k])-(t[k+2]-t[k+1])^2/3/(t[k+2]-t[k])
  i1=x<=t[2]
  sigma[k+1,i1]=x[i1]-t[1]+(t[2]-x[i1])^3/(t[2]-t[1])^2/3
  i2=x>t[2]
  sigma[k+1,i2]=x[i2]-t[1]
  i1=x<=t[k+1]
  sigma[k+2,i1]=0
  i2=x>t[k+1]
  sigma[k+2,i2]=(x[i2]-t[k+1])^3/(t[k+2]-t[k+1])^2/3
  
  v1=1:n*0+1
  v2=x
  x.mat=cbind(v1,v2)
  
  if(pred.new==TRUE){
    list(sigma=sigma,x.mat=x.mat)}
  
  else{
    if(pred.new==FALSE){
      coef=solve(t(x.mat)%*%x.mat)%*%t(x.mat)%*%t(sigma)
      list(sigma=sigma, x.mat=x.mat, center.vector=coef)}
  }
  
}

##############################################################
###
### Basis calculated from the BCGAM Meyer (2008) and bcgam package
###
##############################################################

quant_age <- .25
knots_age <- c(1, round(quantile(right_age - 1, c(seq(quant_age, 
                                                      .99,
                                                      by = quant_age),
                                             .99))), nT_age)
knots_age <- unique(knots_age)
delta_i <- convex(1:nT_age, knots_age, pred.new = FALSE)
delta <- t(rbind(delta_i$sigma - t(delta_i$x.mat %*% delta_i$center.vector)))
delta <- delta / max(delta)
Z_age <- delta
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
     ylim = c(-1, 1),
     main = "Basis Function Age Effect")
for(i in 2:nknots_age){
  lines(1:nT_age,delta[,i])
}
dev.off()

#############################################################
###
### Period effects Basis Function - Kernel convolution
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

