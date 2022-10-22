ageperiod_surv_sim_data_type3 <- function() {

  # sample size of individuals
  n <- 1000

  # intercept of the log hazard
  beta0 <- -4

  ### maximum age interval
  nT_age <- 130

  ### maximum period interval
  nT_period <- 130

  #####################################################
  ### Age at entry. np staggered entry
  #####################################################

  ### For the type3 simlation, all individuals enter at age 1
  left_age <- rep(1, n)

  ###########################################
  ### Staggered entry period effects
  ###########################################

  ### Staggered entry during the first 50 intervals
  interval_ent <- 50
  left_period <- nimble::rcat(n, prob = rep(1 / interval_ent, interval_ent))

  ### The maximum age that an individual can be alive and in the study
  ### must be calculated to restrict the right age during generation
  maxtimes <- nT_period - left_period

  ########################################
  ### Age effects for the log hazard
  ### Baseline Weibull hazard of Age Effects
  ########################################

  lam <- .95
  age <- seq(1, nT_age - 1, by = 1)
  age_effect <- 3 * lam * age^(-(1 - lam))

  #Inducing wiggles during the early age intervals
  anthro <- .1 * cos(1/8 * pi * age)
  anthro[1:50] <- anthro[1:50] * seq(0, 1, length = length(1:50))
  anthro[51:103] <- anthro[51:103] * seq(1, 0, length = length(51:103))
  anthro[103:(nT_age-1)] <- 0
  age_effect <- age_effect + anthro
  age_effect <- age_effect - mean(age_effect)


  ########################################
  ### Period effects for the log hazard
  ########################################

  period <- seq(1, nT_period - 1, by = 1)
  period_effect <- .5 * sin(5/120 * pi * period + 5)
  period_effect <- period_effect - mean(period_effect)


  ########################################
  ### Log hazard
  ########################################

  hazard <- matrix(NA, n, nT_age)

  for (i in 1:n) {
    hazard[i,] <- c(rep(0, left_age[i] - 1), beta0 +
                      age_effect[left_age[i]:(left_age[i] + maxtimes[i] - 1)] +
                      period_effect[left_period[i]:(left_period[i] + maxtimes[i] - 1)],
                      rep(0, nT_age - (left_age[i] - 1 + maxtimes[i])))
  }

  ########################################
  ### Probability of mortality
  ########################################

  test_stat <- matrix(0, n, nT_age)
  for (i in 1:n) {
    for (j in left_age[i]:(left_age[i] + maxtimes[i] - 1)) {
      if (j == left_age[i]) {
        test_stat[i, j] <- (1 - exp(-exp(hazard[i, j])))
      } else {
        test_stat[i, j] <-
           (1 - exp(-exp(hazard[i, j]))) *
           exp(-sum(exp(hazard[i, left_age[i]:(j - 1)])))
      }
    }
    test_stat[i, left_age[i] + maxtimes[i]] <-
            exp(-sum(exp(hazard[i, left_age[i]:(left_age[i]
            + maxtimes[i] - 1)])))
  }

  ########################################
  ### Calculating right censoring
  ########################################
  fail_int <- rep(0, n)
  right_age <- rep(0, n)
  rt_censor <- rep(0, n)
  for (i in 1:n) {
    fail_int[i] <- which(rmultinom(1, 1, test_stat[i, 1:nT_age]) == 1)
    right_age[i] <- ifelse(fail_int[i] >= left_age[i] + maxtimes[i],
                          left_age[i] + maxtimes[i],
                          fail_int[i] + 1)
    rt_censor[i] <- ifelse(fail_int[i] < (left_age[i] + maxtimes[i]), 0, 1)
  }

  right_period <- right_age - left_age + left_period

  ########################################
  ### Return values
  ########################################

  return(list(n = n,
              nT_age = nT_age,
              nT_period = nT_period,
              beta0 = beta0,
              age_effect = age_effect,
              period_effect = period_effect,
              hazard = hazard,
              test_stat = test_stat,
              right_age = right_age,
              left_age = left_age,
              right_period = right_period,
              left_period = left_period,
              rt_censor = rt_censor,
              prop_right_cens = sum(right_age == nT_age)/length(right_age)
              )
  )
}
