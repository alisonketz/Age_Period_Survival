ageperiod_surv_sim_data_type4 <- function() {
  
  # sample size of individuals from each year
  n_yr1 <- 1000
  n_yr2 <- 1000
  n <- n_yr1 + n_yr2
  
  # intercept
  beta0 <- -5

  ### maximum age interval
  nT_age <- 400

  ### maximum period interval
  nT_period <- 104

  #####################################################
  ### Age at entry with staggered entry,
  ### drawn from a stable age distribution
  ### based on the age effect hazard
  ### for 2 'years', and for the age that
  ### individuals that go from birth to censoring
  #####################################################

  #Age effects on log hazard
  #baseline is a quadratic function
  age <- seq(1, nT_age - 1, by = 1) 
  age_effect <- 5e-05 * age^2 - 0.02 * age + 1 

  #Inducing wiggles during the early age intervals
  anthro <- .5 * cos(1 / 26 * pi * age)
  anthro[1:19] <- 0
  anthro[20:38] <- anthro[20:38] * seq(0, 1, length = length(20:38))
  anthro[50:170] <- anthro[50:170] * seq(1, 0, length = length(50:170))
  anthro[171:(nT_age - 1)] <- 0 
  age_effect <- age_effect + anthro
  age_effect <- age_effect - mean(age_effect)

  hazard_se <- -exp((beta0 + age_effect))
  stat_se <- rep(NA, nT_age - nT_period)
  stat_se[1] <- (1 - exp(hazard_se[1]))
  for (j in 2:(nT_age - nT_period)) {
    stat_se[j] <- (1 - exp(hazard_se[j])) * exp(sum(hazard_se[1:(j - 1)]))
  }
  left_age <- nimble::rcat(n, stat_se)

  ###########################################
  # Staggered entry period effects
  ###########################################

  weeks_entry <- 8
  left_yr1 <- nimble::rcat(n_yr1, prob = rep(1 / weeks_entry, weeks_entry))
  left_yr2 <- nimble::rcat(n_yr2, prob = rep(1 / weeks_entry, weeks_entry))
  left_period <- c(left_yr1, left_yr2)

  # no staggered entry
  maxtimes <- nT_period - left_period

  ########################################
  # Age effects for the log hazard
  ########################################

  age <- seq(1, nT_age - 1, by = 1)
  age_effect <- 5e-05 * age^2 - 0.02 * age + 1
  age_effect <- age_effect + anthro
  age_effect <- age_effect - mean(age_effect)
 
  ########################################
  ### Period effects for the log hazard 
  ########################################

  period <- seq(1, nT_period - 1, by = 1)
  period_effect <- 1 * sin(2/52 * pi * (period) + 1)
  period_effect <- period_effect - mean(period_effect)

  ########################################
  ### Log hazard
  ########################################
  hazard <- matrix(NA, n, nT_age)
  for (i in 1:n) {
    hazard[i,] <- c(rep(0,left_age[i] - 1), beta0+
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
      if (j==left_age[i]) {
        test_stat[i, j] <- (1 - exp(-exp(hazard[i, j])))
      } else {
        test_stat[i, j] <-
          (1 - exp(-exp(hazard[i, j])))*exp(-sum(exp(hazard[i, left_age[i]:(j - 1)])))
      }
    }
    test_stat[i, left_age[i] + maxtimes[i]] <-
          exp(-sum(exp(hazard[i, left_age[i]:(left_age[i] + maxtimes[i] - 1)])))
  }

  ########################################
  ### Calculating right censoring
  ########################################
  fail_int <- rep(0, n)
  right_age <- rep(0, n)
  rt_censor <- rep(0, n)
  for(i in 1:n) {
    fail_int[i] <- which(rmultinom(1, 1, test_stat[i,1:nT_age]) == 1)
    right_age[i] <- ifelse(fail_int[i] >= left_age[i] + maxtimes[i],
                           left_age[i] + maxtimes[i], fail_int[i] + 1)
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
              prop_right_cens = sum(right_age == nT_age) / length(right_age)
              )
         )
}
