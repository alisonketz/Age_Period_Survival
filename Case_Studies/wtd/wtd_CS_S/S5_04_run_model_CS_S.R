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
  ### Age effect cgam
  #####################

  for (k in 1:nknots_age_cgam) {
    ln_b_age_cgam[k] ~ dnorm(0, tau_age_cgam)
    b_age_cgam[k] <- exp(ln_b_age_cgam[k])
  }
  tau_age_cgam ~ dgamma(1, 1)
  for (t in 1:nT_age) {
    age_effect_temp[t] <- inprod(b_age_cgam[1:nknots_age_cgam],
                                 Z_age_cgam[t, 1:nknots_age_cgam])
    age_effect_cgam[t] <- age_effect_temp[t] - mu_age_cgam
  }
  mu_age_cgam <- mean(age_effect_temp[1:nT_age])

  #####################
  ### Age effect spline
  #####################

  for (k in 1:nknots_age_spline) {
    b_age_spline[k] ~ dnorm(0, tau_age_spline)
  }
  tau_age_spline ~ dgamma(1, 1)

  for (t in 1:nT_age) {
      age_effect_spline[t] <- inprod(b_age_spline[1:nknots_age_spline],
                                     Z_age_spline[t, 1:nknots_age_spline])
  }

  ########################
  ### Age effect combine
  ########################

  age_effect[1:nT_age] <- age_effect_cgam[1:nT_age] +
                          age_effect_spline[1:nT_age]

  ########################
  ### Period Effects
  ########################
  for (k in 1:nknots_period) {
    b_period[k] ~ dnorm(0, tau_period)
  }
  tau_period ~ dgamma(.01, .01)
 
  for (t in 1:nT_period) {
    period_effect[t] <- inprod(b_period[1:nknots_period],
                               Z_period[t, 1:nknots_period])
  }

  ### Computing state transisiton probability
  SLR[1:records] <- state_transition(records = records,
                                   left = left_age[1:records],
                                   right = right_age[1:records],
                                   nT_age = nT_age,
                                   age_effect = age_effect[1:nT_age],
                                   period_effect = period_effect[1:nT_period],
                                   age2date = age2date[1:records],
                                   beta0 = beta0)

  for (j in 1:records) {
    censor[j] ~ dbern(SLR[j])
  }

  ##########################
  ### Derived parameters
  ##########################

  for (t in 1:nT_age) {
    llambda_age[t] <- beta0 + age_effect[t]
    UCH0_age[t] <- exp(llambda_age[t])
    S0_age[t] <- exp(-sum(UCH0_age[1:t]))
  }
  for (t in 1:nT_period) {
    llambda_period[t] <- beta0 + period_effect[t]
    UCH0_period[t] <- exp(llambda_period[t])
    S0_period[t] <- exp(-sum(UCH0_period[1:t]))
  }

})#end model statement

#Data
nimData <- list(censor = d_fit$censor,
                Z_period = Z_period,
                Z_age_cgam = Z_age_cgam,
                Z_age_spline = Z_age_spline,
                left_age = d_fit$ageCapture,
                right_age = d_fit$ageRight,
                age2date = age2date
                )

nimConsts <-list(records = n_fit,
                 nT_age = nT_age,
                 nT_period = nT_period,
                 nknots_age_cgam = nknots_age_cgam,
                 nknots_age_spline = nknots_age_spline,
                 nknots_period = nknots_period
)

initsFun <- function()list(
                          tau_period = runif(1, .1, 1),
                          tau_age_cgam = runif(1, .1, 1),
                          tau_age_spline = runif(1, .1, 1),
                          beta0_temp = rnorm(1, -5, .0001),
                          mix = 1,
                          b_period = rnorm(nknots_period) * 10^-4,
                          b_age_spline = rnorm(nknots_age_spline) * 10^-4,
                          ln_b_age_cgam = runif(nknots_age_cgam, -10, -5)
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
  "S0_age",
  "llambda_age",
  "age_effect",
  "age_effect_cgam",
  "age_effect_spline",
  "tau_age_cgam",
  "tau_age_spline",
  "b_age_cgam",
  "b_age_spline",
  "b_period",
  "period_effect",
  "S0_period",
  "llambda_period",
  "tau_period"
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

