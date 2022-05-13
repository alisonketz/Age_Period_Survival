##################################################################
###
###  Running Age-Period hazard model
###
##################################################################

##################################################################
###
###  Function to calculate the
###  probability in the likelihood
###
##################################################################

state_transition <- nimbleFunction(
  run = function(records = double(0),
                 left=double(1),
                 right = double(1),
                 age_effect = double(1),
                 period_effect = double(1),
                 age2date = double(1),
                 nT_age = double(0),
                 beta0 = double(0)
  ){

    SLR <- nimNumeric(records)
    UCH <- nimMatrix(value = 0, nrow = records, ncol = nT_age)

    for (j in 1:records) {
      for (k in left[j]:(right[j] - 1)) {
        UCH[j, k] <- exp(beta0 +
                         age_effect[k] +
                         period_effect[k + age2date[j]])
      }
      SLR[j] <- exp(-sum(UCH[j, left[j]:(right[j] - 1)]))
    }
    returnType(double(1))
    return(SLR[1:records])
})

#compile state transition probability function
Cstate_transition <- compileNimble(state_transition)
 

###########################
###
### Run model
###
###########################

modelcode <- nimbleCode({

  #Priors for Age and Period effects
  beta0_temp ~ dnorm(0, .01)
  mix ~ dunif(-1, 1)
  beta0 <- beta0_temp * mix

  #####################
  ### Age Effects
  #####################

  for (k in 1:nknots_age) {
    ln_b_age[k] ~ dnorm(0, tau_age)
    b_age[k] <- exp(ln_b_age[k])
  }
  tau_age ~ dgamma(1, 1)
  
  for (t in 1:nT_age) {
    age_effect_temp[t] <- inprod(b_age[1:nknots_age],Z_age[t,1:nknots_age])
    age_effect[t] <- age_effect_temp[t] - mu_age
  }
  mu_age <- mean(age_effect_temp[1:nT_age])

  ##################################################
  ### Period effects
  ### for entire study (3 years)
  ### with separate period effects for each year,
  ### but shared hyperparams
  ##################################################### 

  mix2 ~ dunif(-1, 1)
  ln_sk_period ~ dnorm(0, sd = 1)
  sdk_period <- exp(mix2 * ln_sk_period)
  tauk_period <- 1 / sdk_period^2
  stauk_period <- sqrt(tauk_period)
  sda_period ~ T(dnorm(0, sd = 1), 0, Inf)#<- 1/sqrt(taua_period)
  taua_period <- 1 / sda_period^2
  ratioinf_period <- sdk_period / sda_period #ratio of variability
 
  ###
  ### Period effects for year 1
  ###
 
  for (i in 1:(nknots_period_yr1)) {
    alpha_period_yr1[i] ~ dnorm(0, 1)
    alphau_period_yr1[i] <- sda_period * alpha_period_yr1[i]
  }
  period_effect_yr1[1:nT_yr1] <- kernel_conv(
    nT = nT_yr1,
    Z = Z_period_yr1[1:nT_yr1, 1:nknots_period_yr1],
    stauk = stauk_period,
    nconst = nconst,
    tauk = tauk_period,
    nknots = nknots_period_yr1,
    alphau = alphau_period_yr1[1:nknots_period_yr1]
  )

  ###
  ### Period effects for year 2
  ###
 
  for (i in 1:(nknots_period_yr2)) {
    alpha_period_yr2[i] ~ dnorm(0, 1)
    alphau_period_yr2[i] <- sda_period * alpha_period_yr2[i]
  }

  period_effect_yr2[1:nT_yr2] <- kernel_conv(
    nT = nT_yr2,
    Z = Z_period_yr2[1:nT_yr2, 1:nknots_period_yr2],
    stauk = stauk_period,
    nconst = nconst,
    tauk = tauk_period,
    nknots = nknots_period_yr2,
    alphau = alphau_period_yr2[1:nknots_period_yr2]
  )

  ###
  ### Period effects for year 3
  ###

  for (i in 1:(nknots_period_yr3)) {
    alpha_period_yr3[i] ~ dnorm(0, 1)
    alphau_period_yr3[i] <- sda_period * alpha_period_yr3[i]
  }

   period_effect_yr3[1:nT_yr3] <- kernel_conv(
    nT = nT_yr3,
    Z = Z_period_yr3[1:nT_yr3, 1:nknots_period_yr3],
    stauk = stauk_period,
    nconst = nconst,
    tauk = tauk_period,
    nknots = nknots_period_yr3,
    alphau = alphau_period_yr3[1:nknots_period_yr3]
  )

  ###
  ### Combine period effects into a single vector
  ###

  period_effect[1:nT_yr1] <- period_effect_yr1[1:nT_yr1]
  period_effect[(nT_yr1 + 1):(nT_yr1 + nT_yr2)] <- period_effect_yr2[1:nT_yr2]
  period_effect[(nT_yr1 + nT_yr2 + 1):(nT_period_total)] <- period_effect_yr3[1:nT_yr3]
    
  ### Computing state transisiton probability
  SLR[1:records] <- state_transition(records = records,
                                   left = left_age[1:records],
                                   right = right_age[1:records],
                                   nT_age = nT_age,
                                   age_effect = age_effect[1:nT_age],
                                   period_effect = period_effect[1:nT_period_total],
                                   age2date = age2date[1:records],
                                   beta0 = beta0)

  for (j in 1:records) {
    censor[j] ~ dbern(SLR[j])
  }

  ##########################
  ### Derived parameters
  ##########################
  #age
  for (t in 1:nT_age) {
    llambda_age[t] <- beta0 + age_effect[t]
    UCH0_age[t] <- exp(llambda_age[t])
    S0_age[t] <- exp(-sum(UCH0_age[1:t]))
  }
  #period
  for (t in 1:nT_yr1) {
    llambda_period[t, 1] <- beta0 + period_effect_yr1[t]
    UCH0_period[t, 1] <- exp(llambda_period[t, 1])
    S0_period[t, 1] <- exp(-sum(UCH0_period[1:t, 1]))
  }
  for (t in 1:(nT_yr2)) {
    llambda_period[t, 2] <- beta0+period_effect_yr2[t]
    UCH0_period[t, 2] <- exp(llambda_period[t, 2])
    S0_period[t, 2] <- exp(-sum(UCH0_period[1:t, 2]))
  }
  for (t in 1:(nT_yr3)) {
    llambda_period[t, 3] <- beta0 + period_effect_yr3[t]
    UCH0_period[t, 3] <- exp(llambda_period[t, 3])
    S0_period[t, 3] <- exp(-sum(UCH0_period[1:t, 3]))
  }
  #combining age+period
  for (t in 1:nT_age) {#nT_age<nT_yr1
    llambda_ageperiod[t, 1] <- beta0 + period_effect_yr1[t] + age_effect[t]
    UCH0_ageperiod[t, 1] <- exp(llambda_ageperiod[t, 1])
    S0_ageperiod[t, 1] <- exp(-sum(UCH0_ageperiod[1:t, 1]))
  }
  for (t in 1:(nT_yr2)) {
    llambda_ageperiod[t, 2] <- beta0 + period_effect_yr2[t] + age_effect[t]
    UCH0_ageperiod[t, 2] <- exp(llambda_ageperiod[t, 2])
    S0_ageperiod[t, 2] <- exp(-sum(UCH0_ageperiod[1:t, 2]))
  }
  for (t in 1:(nT_yr3)) {
    llambda_ageperiod[t, 3] <- beta0 + period_effect_yr3[t] + age_effect[t]
    UCH0_ageperiod[t, 3] <- exp(llambda_ageperiod[t, 3])
    S0_ageperiod[t, 3] <- exp(-sum(UCH0_ageperiod[1:t, 3]))
  }

})#end model statement

#Data
nimData <- list(censor = d_fit$censor,
                Z_period_yr1 = Z_period_yr1,
                Z_period_yr2 = Z_period_yr2,
                Z_period_yr3 = Z_period_yr3,
                Z_age = Z_age,
                left_age = d_fit$left.age,
                right_age = d_fit$right.age,
                age2date = age2date,
                llambda_period = matrix(NA, nr = nT_period_max, nc = 3),
                UCH0_period = matrix(NA, nr = nT_period_max, nc = 3),
                S0_period = matrix(NA, nr = nT_period_max, nc = 3),
                llambda_ageperiod = matrix(NA, nr = nT_period_max, nc = 3),
                UCH0_ageperiod = matrix(NA, nr = nT_period_max, nc = 3),
                S0_ageperiod = matrix(NA, nr = nT_period_max, nc = 3),
                period_effect = rep(NA,nT_period_total)
                )

nimConsts <- list(records = n_fit,
                 nT_age = nT_age,
                 nT_period_total = nT_period_total,
                 nT_yr1 = nT_yr1,
                 nT_yr2 = nT_yr2,
                 nT_yr3 = nT_yr3,
                 nknots_age = nknots_age,
                 nknots_period_yr1 = nknots_period_yr1,
                 nknots_period_yr2 = nknots_period_yr2,
                 nknots_period_yr3 = nknots_period_yr3,
                 nconst = 1 / sqrt(2 * pi)
                 )

initsFun <- function()list(beta0_temp = rnorm(1, -5, .0001),
                          mix = 1,
                          mix2 = 1,
                          sda_period = runif(1, 0, 5),
                          ln_sk_period = rnorm(1, 0, 1),
                          alpha_period_yr1 = rep(0, nknots_period_yr1),
                          alpha_period_yr2 = rep(0, nknots_period_yr2),
                          alpha_period_yr3 = rep(0, nknots_period_yr3),
                          tau_age = runif(1, .1, 1),
                          ln_b_age = rnorm(nknots_age) * 10^-4,
                          llambda_period = matrix(0,
                                                  nr = nT_period_max,
                                                  nc = 3),
                          UCH0_period = matrix(0,
                                               nr = nT_period_max,
                                               nc = 3),
                          S0_period = matrix(0,
                                             nr = nT_period_max,
                                             nc = 3),
                          llambda_ageperiod = matrix(0,
                                                     nr = nT_period_max,
                                                     nc = 3),
                          UCH0_ageperiod = matrix(0,
                                                  nr = nT_period_max,
                                                  nc = 3),
                          S0_ageperiod = matrix(0,
                                                nr = nT_period_max,
                                                nc = 3),
                          period_effect = rep(0, nT_period_total)
                          )
nimInits <- initsFun()

Rmodel <- nimbleModel(code = modelcode,
                      constants = nimConsts,
                      data = nimData,
                      inits = initsFun()
                      )

#identify params to monitor
parameters <- c(
              "beta0",
              "beta0_temp",
              "mix",
              "mix2",
              "S0_age",
              "llambda_age",
              "age_effect",
              "period_effect_yr1",
              "period_effect_yr2",
              "period_effect_yr3",
              "S0_period",
              "llambda_period",
              "sdk_period",
              "sda_period",
              "alphau_period_yr1",
              "alphau_period_yr2",
              "alphau_period_yr3",
              "tau_age",
              "b_age",
              "S0_ageperiod",
              "llambda_ageperiod"
)
starttime <- Sys.time()
confMCMC <- configureMCMC(Rmodel,
                          monitors = parameters,
                          thin = n_thin,
                          useConjugacy = FALSE,
                          enableWAIC = TRUE)
nimMCMC <- buildMCMC(confMCMC)
Cnim <- compileNimble(Rmodel)
CnimMCMC <- compileNimble(nimMCMC,
                          project = Rmodel)
mcmcout <- runMCMC(CnimMCMC,
                   niter = reps,
                   nburnin = bin,
                   nchains = n_chains,
                   inits = initsFun,
                   samplesAsCodaMCMC = TRUE,
                   summary = TRUE,
                   WAIC = TRUE
                   )

runtime <- difftime(Sys.time(),
                    starttime,
                    units = "min")

###
### save model run
###

save(runtime, file = "runtime.Rdata")
save(mcmcout, file = "mcmcout.Rdata")

