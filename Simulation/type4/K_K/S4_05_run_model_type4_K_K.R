##################################################
###
###  Function state transition probability 
###  Age-Period Model
###
###################################################

state_transition <- nimbleFunction(
  run = function(records = double(0),
                 left = double(1),
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
                         period_effect[k - age2date[j]])
      }
      SLR[j] <- exp(-sum(UCH[j, left[j]:(right[j] - 1)]))
    }
    returnType(double(1))
    return(SLR[1:records])
  })

cstate_transition <- compileNimble(state_transition)

###########################
###
### Model Specification
###
###########################

modelcode <- nimbleCode({

  ### Prior for intercept
  ### using parameter expansion for convergence
  beta0_temp ~ dnorm(0, .01)
  mix ~ dunif(-1, 1)
  beta0 <- beta0_temp * mix

  
  ###
  ### Age Effect
  ###

  mix1 ~ dunif(-1, 1)
  ln_sk_age ~ dnorm(0, sd=3)
  sdk_age <- exp(mix1 * ln_sk_age)
  tauk_age <- 1 / sdk_age^2
  stauk_age <- sqrt(tauk_age)
  sda_age ~ T(dnorm(0, sd=1), 0, Inf)
  tau_age <- 1 / sda_age^2
  for (i in 1:nknots_age) {
    alpha_age[i] ~ dnorm(0, 1)
    alphau_age[i] <- sda_age * alpha_age[i]
  }
  ratioinf_age <- sdk_age / sda_age #ratio of variability

  age_effect[1:nT_age] <- kernel_conv(
    nT = nT_age,
    Z = Z_age[1:nT_age, 1:nknots_age],
    stauk = stauk_age,
    nconst = nconst,
    tauk = tauk_age,
    nknots = nknots_age,
    alphau = alphau_age[1:nknots_age]
  )

  ###
  ### time effect for entire study (2 years)
  ###

  #Kernel convolution- time to death
  #time line

  mix2 ~ dunif(-1, 1)
  ln_sk_period ~ dnorm(0, sd = 1)
  sdk_period <- exp(mix2 * ln_sk_period)
  tauk_period <- 1 / sdk_period^2
  stauk_period <- sqrt(tauk_period)
  sda_period ~ T(dnorm(0, sd = 1), 0, Inf)
  taua_period <- 1 / sda_period^2
  for(i in 1:(nknots_period)) {
    alpha_period[i] ~ dnorm(0, 1)
    alphau_period[i] <- sda_period * alpha_period[i]
  }
  ratioinf_period <- sdk_period / sda_period #ratio of variability

  period_effect[1:nT_period] <- kernel_conv(
    nT = nT_period,
    Z = Z_period[1:nT_period, 1:nknots_period],
    stauk = stauk_period,
    nconst = nconst,
    tauk = tauk_period,
    nknots = nknots_period,
    alphau = alphau_period[1:nknots_period]
  )

  ###
  ### Computing state transition probability
  ###

  SLR[1:records] <- state_transition(records = records,
                                   left = left_age[1:records],
                                   right = right_age[1:records],
                                   nT_age = nT_age,
                                   age_effect = age_effect[1:nT_age],
                                   period_effect = period_effect[1:nT_period],
                                   age2date = age2date[1:records],
                                   beta0 = beta0)

  ###
  ### Likelihood
  ###
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

#data
nimData <- list(censor = df_fit[, 3],
                Z_age = Z_age,
                Z_period = Z_period,
                left_age = df_fit[, 1],
                right_age = df_fit[, 2],
                age2date = age2date
                )

nimConsts <- list(records = n_fit,
                 nT_age = nT_age,
                 nT_period = nT_period,
                 nknots_age = ncol(Z_age),
                 nknots_period = ncol(Z_period),
                 nconst = 1 / sqrt(2 * pi))

#specify initial values
initsFun <- function()list(
                          beta0_temp = rnorm(1, beta0_true, .0001),
                          mix = 1,
                          mix1 = 1,
                          mix2 = 1,
                          sda_period = runif(1, 0, 5),
                          ln_sk_period = rnorm(1, 0, 1),
                          alpha_period = rep(0, nknots_period),
                          sda_age = runif(1, 0, 5),
                          ln_sk_age = rnorm(1, 0, 1),
                          alpha_age = rep(0, nknots_age)
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
              "S0_age",
              "llambda_age",
              "age_effect",
              "sdk_age",
              "sda_age",
              "alpha_age",
              "period_effect",
              "S0_period",
              "llambda_period",
              "sdk_period",
              "sda_period",
              "alpha_period",
              "ratioinf_age",
              "ratioinf_period"
)

starttime <- Sys.time()
confMCMC <- configureMCMC(Rmodel,
                          monitors = parameters,
                          thin = n_thin,
                          useConjugacy = FALSE,
                          enableWAIC = TRUE)
nimMCMC <- buildMCMC(confMCMC)
Cnim <- compileNimble(Rmodel)
CnimMCMC <- compileNimble(nimMCMC, project = Rmodel)
mcmcout <- runMCMC(CnimMCMC,
                   niter = reps,
                   nburnin = bin,
                   nchains = n_chains,
                   inits = initsFun,
                   samplesAsCodaMCMC = TRUE,
                   summary = TRUE,
                   WAIC = TRUE)

runtime <- difftime(Sys.time(), starttime, units = "min")

save(runtime, file = "runtime.Rdata")
save(mcmcout, file = "mcmcout.Rdata")
