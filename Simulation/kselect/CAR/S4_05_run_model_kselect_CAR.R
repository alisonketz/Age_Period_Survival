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

  # Priors
  ###
  ### ICAR Age Effect
  ###
  
  age_effect[1:nage] ~ dcar_normal(adj_age[1:nNage],
                                   weights = weights_age[1:nNage],
                                   num_age[1:nage],
                                   tau_age)
  tau_age ~ dgamma(1, 1)
  mu_age <- mean(age_effect[1:nage])

  #ICAR time effect
  period_effect[1:nperiod] ~ dcar_normal(adj_period[1:nNperiod],
                                       weights = weights_period[1:nNperiod],
                                       num_period[1:nperiod],
                                       tau_period,
                                       zero_mean = 1)
  tau_period ~ dgamma(1, 1)

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
    llambda_age[t] <- age_effect[t]
    UCH0_age[t] <- exp(llambda_age[t])
    S0_age[t] <- exp(-sum(UCH0_age[1:t]))
  }
  for (t in 1:nT_period) {
    llambda_period[t] <- mu_age + period_effect[t]
    UCH0_period[t] <- exp(llambda_period[t])
    S0_period[t] <- exp(-sum(UCH0_period[1:t]))
  }

})#end model statement

#data
nimData <- list(censor = df_fit[, 3],
                left_age = df_fit[, 1],
                right_age = df_fit[, 2],
                age2date = age2date
                )

nimConsts <- list(records = n_fit,
                 nT_age = nT_age,
                 nT_period = nT_period,
                 nNperiod = nNperiod,
                 adj_period = adj_period,
                 weights_period = weights_period,
                 nperiod = nperiod,
                 num_period = num_period,
                 nNage = nNage,
                 adj_age = adj_age,
                 weights_age = weights_age,
                 nage = nage,
                 num_age = num_age
)

#specify initial values

initsFun <- function()list(period_effect = rep(.1, nT_period),
                          age_effect = rep(-6, nT_age),
                          tau_age = runif(1, 5, 10),
                          tau_period = runif(1, 2, 5)
                          )
nimInits <- initsFun()


Rmodel <- nimbleModel(code = modelcode,
                      constants = nimConsts,
                      data = nimData,
                      inits = initsFun()
                      )

#identify params to monitor
parameters <- c("S0_age",
                # "llambda_age",
                "age_effect",
                "tau_age",
                "period_effect",
                "S0_period",
                # "llambda_period",
                "tau_period",
                "mu_age"
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