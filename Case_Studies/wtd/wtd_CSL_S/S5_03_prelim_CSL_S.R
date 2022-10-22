#######################################################
###
### Preliminary constants
###
#######################################################

#########################################################################################################################3
###
### Determining basis function and number of knots 
###
#########################################################################################################################3

# Function for Basis - Convex shape neither decreasing or increasing
convex=function(x, t, pred.new=TRUE){
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
##############################################################
###
### Basis calculated from the BCGAM Meyer (2008) and bcgam package
###
##############################################################

quant_age <- .2

knots_age_cgam <- c(1, round(quantile(d_fit$ageRight - 1,
                       c(seq(quant_age, .99, by = quant_age),
                       .99))))
knots_age_cgam <- unique(knots_age_cgam)
delta_i <- convex(1:nT_age, knots_age_cgam, pred.new = FALSE)
delta <- t(rbind(delta_i$sigma - t(delta_i$x.mat %*% delta_i$center.vector) )  )
delta <- delta / max(delta)
Z_age_cgam <- delta
nknots_age_cgam<- dim(Z_age_cgam)[2]


###############################
###
### basis functions
### using the spline package
###
###############################

quant_age <- .05
knots_age_spline <- c(1,
                   round(quantile(d_fit$ageRight - 1,
                                  c(seq(quant_age, .99, by = quant_age), .99))))
knots_age_spline <- unique(knots_age_spline)

##Basis for age hazard
splinebasis <- ns(1:nT_age, knots = knots_age_spline)

##A constraint matrix so time.effects = 0
constr_sumzero <- matrix(1, 1, nrow(splinebasis)) %*% splinebasis

##Get a basis for null space of constraint
qrc <- qr(t(constr_sumzero))
Z <- qr.Q(qrc, complete = TRUE)[, (nrow(constr_sumzero) + 1):ncol(constr_sumzero)]
Z_age_spline <- splinebasis %*% Z
nknots_age_spline <- dim(Z_age_spline)[2]

#############################################################
###
### plot of the basis functions
###
##############################################################
pdf("basis_function_age.pdf")
plot(1:nT_age,
     Z_age_cgam[, 1],
     ylim = c(-1, 1),
     type = "l",
     main = "Basis Function Age Effect CGAM")
for (i in 2:nknots_age_cgam) {
  lines(1:nT_age, Z_age_cgam[, i])
}

plot(1:nT_age,
     Z_age_spline[, 1],
     ylim = c(-1, 1),
     type = "l",
     main = "Basis Function Age Effect Spline")
for (i in 2:nknots_age_spline) {
  lines(1:nT_age, Z_age_spline[, i])
}
dev.off()

#########################################
###
### Setting up the distance matrix
### for the kernel convolution on Time
###
#########################################
########################################
###
### Spline basis matrix for Period
###
##########################################

intvl_period <- 12
knots_period <- c(1,
               seq(5,
                   nT_period-1,
                   by = intvl_period))
knots_period <- unique(knots_period)
splinebasis <- bs(1:nT_period, knots = knots_period)

##A constraint matrix so period.effects = 0
constr_sumzero <- matrix(1, 1, nrow(splinebasis)) %*% splinebasis

##Get a basis for null space of constraint
qrc <- qr(t(constr_sumzero))
Z <- qr.Q(qrc,
          complete = TRUE)[, (nrow(constr_sumzero) + 1):ncol(constr_sumzero)]
Z_period <- splinebasis %*% Z
nknots_period <- dim(Z_period)[2]


######################################################
###
### number of MCMCr iterations, Chains, and Burn-in
###
######################################################

#number of MCMCr iterations, Chains, and Burn-in
reps <- 50000
bin <- reps * .5
n_chains <- 3
n_thin <- 1
