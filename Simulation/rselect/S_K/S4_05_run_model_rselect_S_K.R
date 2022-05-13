##################################################
###
###  Running Model using the Nimble package
###  Age-Period Model
###
##################################################
##################################################
###
###  Function state transition probability
###
##################################################

#function to calculate cumulative survival probability for likelihood
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
    UCH <-nimMatrix(value = 0, nrow = records, ncol = nT_age)

    for (j in 1:records) {
      for (k in left[j]:(right[j] - 1)) {
        UCH[j, k] <- exp(beta0 + age_effect[k] + period_effect[k - age2date[j]])
      }
      # total prob of surviving
      SLR[j] <- exp(-sum(UCH[j, left[j]:(right[j] - 1)]))
    }

    returnType(double(1))
    return(SLR[1:records])
})

Cstate_transition <- compileNimble(state_transition)
 
##################################################
###
### Model Statement
###
##################################################

modelcode <- nimbleCode({

  ### Prior for intercept
  ### using parameter expansion for convergence
  beta0_temp ~ dnorm(0, 0.01)
  mix ~ dunif(-1, 1)
  beta0 <- beta0_temp * mix

  #####################
  ### Age effect spline
  #####################
  
  for (k in 1:nknots_age) {
    b_age[k] ~ dnorm(0, tau_age)
  }
  tau_age ~ dgamma(0.01, 0.01)
    for (t in 1:nT_age) {
    age_effect[t] <- inprod(b_age[1:nknots_age], Z_age[t, 1:nknots_age])
  }

 #####################################
  ### Period effect kernel convolution
  #####################################

  mix2 ~ dunif(-1, 1)
  ln_sk_period ~ dnorm(0, sd = 1)
  sdk_period <- exp(mix2 * ln_sk_period)
  tauk_period <- 1 / sdk_period^2
  stauk_period <- sqrt(tauk_period)
  sda_period ~ T(dnorm(0, sd = 1), 0, Inf)
  taua_period <- 1 / sda_period^2
  for (i in 1:(nknots_period)) {
    alpha_period[i] ~ dnorm(0, 1)
    alphau_period[i] <- sda_period * alpha_period[i]
  }
  ratioinf_period <- sdk_period / sda_period #ratio of variability/smoothing

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

###Data
nimData <- list(censor = df_fit[, 3],
                Z_period = Z_period,
                Z_age = Z_age,
                left_age = df_fit[, 1],
                right_age = df_fit[, 2],
                age2date = age2date
                )

###Constants
nimConsts <- list(records = n_fit,
                 nT_age = nT_age,
                 nT_period = nT_period,
                 nknots_age = nknots_age,
                 nknots_period = nknots_period,
                 nconst = 1 / sqrt(2 * pi))

### Initial values
initsFun <- function()list(
                          beta0_temp = rnorm(1, beta0_true, .0001),
                          mix = 1,
                          mix2 = 1,
                          sda_period = runif(1, 0, 5),
                          ln_sk_period = rnorm(1, 0, 1),
                          alpha_period = rep(0, nknots_period),
                          tau_age = runif(1),
                          b_age = rnorm(nknots_age) * 10^-4
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
              "tau_age",
              "period_effect",
              "S0_period",
              "llambda_period",
              "sdk_period",
              "sda_period",
              "alpha_period",
              "ratioinf_period"
)

starttime <- Sys.time()
confMCMC <- configureMCMC(Rmodel,
                          monitors = parameters,
                          thin = n_thin,
                          useConjugacy = FALSE,
                          enableWAIC = TRUE)
confMCMC$removeSamplers(c("b_age"))
confMCMC$addSampler(target = "b_age", type = "RW_block")
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
