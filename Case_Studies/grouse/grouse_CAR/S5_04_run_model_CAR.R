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

  #####################
  ### Age Effects
  #####################

  age_effect[1:nT_age] ~ dcar_normal(adj_age[1:nNage],
                                     weights = weights_age[1:nNage],
                                     num_age[1:nT_age],
                                     tau_age)
  tau_age ~ dgamma(1, 1)
  mu_age <- mean(age_effect[1:nT_age])

  #####################
  ### Period Effects
  #####################

  period_effect_yr1[1:nT_yr1] ~ dcar_normal(adj_period_yr1[1:nNperiod_yr1],
                                       weights = weights_period_yr1[1:nNperiod_yr1],
                                       num_period_yr1[1:nT_yr1],
                                       tau_period,
                                       zero_mean=1)

  period_effect_yr2[1:nT_yr2] ~ dcar_normal(adj_period_yr2[1:nNperiod_yr2],
                                       weights=weights_period_yr2[1:nNperiod_yr2],
                                       num_period_yr2[1:nT_yr2],
                                       tau_period,
                                       zero_mean=1)
  period_effect_yr3[1:nT_yr3] ~ dcar_normal(adj_period_yr3[1:nNperiod_yr3],
                                       weights=weights_period_yr3[1:nNperiod_yr3],
                                       num_period_yr3[1:nT_yr3],
                                       tau_period,
                                       zero_mean=1)
  tau_period ~ dgamma(1, 1)

  period_effect[1:nT_period_total] <- c(period_effect_yr1[1:nT_yr1],
                            period_effect_yr2[1:nT_yr2],
                            period_effect_yr3[1:nT_yr3])

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
    llambda_period[t, 1] <- mu_age + period_effect_yr1[t]
    UCH0_period[t, 1] <- exp(llambda_period[t, 1])
    S0_period[t, 1] <- exp(-sum(UCH0_period[1:t, 1]))
  }
  for (t in 1:(nT_yr2)) {
    llambda_period[t, 2] <- mu_age + period_effect_yr2[t]
    UCH0_period[t, 2] <- exp(llambda_period[t, 2])
    S0_period[t, 2] <- exp(-sum(UCH0_period[1:t, 2]))
  }
  for (t in 1:(nT_yr3)) {
    llambda_period[t, 3] <- mu_age + period_effect_yr3[t]
    UCH0_period[t, 3] <- exp(llambda_period[t, 3])
    S0_period[t, 3] <- exp(-sum(UCH0_period[1:t, 3]))
  }
  #combining age+period
  for (t in 1:nT_age) {#nT_age<nT_yr1
    llambda_ageperiod[t, 1] <- period_effect_yr1[t] + age_effect[t]
    UCH0_ageperiod[t, 1] <- exp(llambda_ageperiod[t, 1])
    S0_ageperiod[t, 1] <- exp(-sum(UCH0_ageperiod[1:t, 1]))
  }
  for (t in 1:(nT_yr2)) {
    llambda_ageperiod[t, 2] <- period_effect_yr2[t] + age_effect[t]
    UCH0_ageperiod[t, 2] <- exp(llambda_ageperiod[t, 2])
    S0_ageperiod[t, 2] <- exp(-sum(UCH0_ageperiod[1:t, 2]))
  }
  for (t in 1:(nT_yr3)) {
    llambda_ageperiod[t, 3] <- period_effect_yr3[t] + age_effect[t]
    UCH0_ageperiod[t, 3] <- exp(llambda_ageperiod[t, 3])
    S0_ageperiod[t, 3] <- exp(-sum(UCH0_ageperiod[1:t, 3]))
  }

})#end model statement

#Data
nimData <- list(censor = d_fit$censor,
                left_age = d_fit$left.age,
                right_age = d_fit$right.age,
                age2date = age2date,
                adj_age = adj_age,
                weights_age = weights_age,
                num_age = num_age,
                adj_period_yr1 = adj_period_yr1,
                weights_period_yr1 = weights_period_yr1,
                num_period_yr1 = num_period_yr1,
                adj_period_yr2 = adj_period_yr2,
                weights_period_yr2 = weights_period_yr2,
                num_period_yr2 = num_period_yr2,
                adj_period_yr3 = adj_period_yr3,
                weights_period_yr3 = weights_period_yr3,
                num_period_yr3 = num_period_yr3,
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
                 nNage=nNage,
                 nNperiod_yr1 = nNperiod_yr1,
                 nNperiod_yr2 = nNperiod_yr2,
                 nNperiod_yr3 = nNperiod_yr3,
                 nT_yr1 = nT_yr1,
                 nT_yr2 = nT_yr2,
                 nT_yr3 = nT_yr3,
                 beta0 = 0
                 )

initsFun <- function()list(
                          period_effect_yr1 = rep(0, nT_yr1),
                          period_effect_yr2 = rep(0, nT_yr2),
                          period_effect_yr3 = rep(0, nT_yr3),
                          age_effect = rep(-6, nT_age),
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
              "S0_age",
              "llambda_age",
              "age_effect",
              "period_effect_yr1",
              "period_effect_yr2",
              "period_effect_yr3",
              "S0_period",
              "llambda_period",
              "tau_age",
              "tau_period",
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

