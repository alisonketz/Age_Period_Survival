out <- mcmcout$samples
fit_sum <- mcmcout$summary$all.chains


pdf("traceplots_type3_CAR.pdf")
traceplot(out[, "mu_age"], ylab = "mu_age")
traceplot(out[, "tau_age"], ylab = "tau_age")
traceplot(out[, "tau_period"], ylab = "tau_period")
traceplot(out[, "age_effect[1]"], ylab = "age_effect[1]")
traceplot(out[, "age_effect[3]"], ylab = "age_effect[3]")
traceplot(out[, "age_effect[10]"], ylab = "age_effect[10]")
traceplot(out[, "age_effect[51]"], ylab = "age_effect[51]")
traceplot(out[, "period_effect[1]"], ylab = "period_effect[1]")
traceplot(out[, "period_effect[3]"], ylab = "period_effect[3]")
traceplot(out[, "period_effect[10]"], ylab = "period_effect[10]")
traceplot(out[, "period_effect[51]"], ylab = "period_effect[51]")
dev.off()

# The Colorblind palette with grey:
cbPalette <- c("#0072B2",
               "#D55E00",
               "#CC79A7",
               "#999999",
               "#E69F00",
               "#56B4E9")

study_end <- nT_age
S0_age_indx <- grep("S0_age", rownames(fit_sum))
S0_age <- fit_sum[S0_age_indx, ]
intervals <- 1:study_end

survival <- c(S0_age[1:study_end,1])
survival_lower <- c(S0_age[1:study_end,4])
survival_upper <- c(S0_age[1:study_end,5])
survival_true <- c(S0_age_true)

out_s_age <- data.frame(intervals,
                        survival,
                        survival_lower,
                        survival_upper,
                        survival_true)

p_s_age <- ggplot(data = out_s_age, aes(x = intervals)) +
  geom_ribbon(aes(ymin = survival_lower, ymax = survival_upper),
              alpha = .3,
              linetype = 0) +
  geom_line(aes(x = intervals, y = survival), linetype = 1, size = 2) +
  geom_line(aes(x = intervals,y=survival_true),
            linetype = 2,
            size = 1,
            color = cbPalette[5]) +
  ylim(0,1) +
  ggtitle("Survival - Age") +
  xlab("Age") +
  ylab("Survival Probability") +
  theme_bw()

ggsave("S0_age.pdf",p_s_age)

##################################
###
### Log hazard - Age effect cgam
###
##################################

te_indx <- grep("age_effect", rownames(fit_sum))[1:(nT_age)]

age_effect_mean <- fit_sum[te_indx, 1]
age_effect_lower <- fit_sum[te_indx, 4]
age_effect_upper <- fit_sum[te_indx, 5]

intervals <- 1:nT_age

out_age_effect <- data.frame(intervals,
                            age_effect_mean,
                            age_effect_lower,
                            age_effect_upper,
                            truth = ageperiod_out$age_effect[1:(nT_age)] +
                                    beta0_true
                            )

age_effect_plot <- ggplot(data = out_age_effect, aes(x = intervals)) +
  geom_line(aes(x = intervals, y = age_effect_mean),
            size = 1) +
  geom_ribbon(aes(ymin = age_effect_lower, ymax = age_effect_upper),
              alpha = .2,
              linetype = 0) +
  geom_line(aes(x = intervals, y = truth),
            size = 1,
            col = "red") +
  ggtitle("Age Effect Posterior") +
  xlab("Age") +
  ylab("Effect Size") +
  theme_bw()

ggsave("age_effect.pdf", age_effect_plot)

##################################
###
### Survival - Period
###
##################################

study_end <- nT_period
S0_period_indx <- grep("S0_period", rownames(fit_sum))
S0_period <- fit_sum[S0_period_indx, ]
intervals <- 1:study_end

survival <- S0_period[1:study_end,1]
survival_lower <- S0_period[1:study_end,4]
survival_upper <- S0_period[1:study_end,5]
survival_true <- c(S0_period_true)

out_s_period <- data.frame(intervals,
                        survival,
                        survival_lower,
                        survival_upper,
                        survival_true)

p_s_period <- ggplot(data = out_s_period, aes(x = intervals)) +
  geom_ribbon(aes(ymin = survival_lower, ymax = survival_upper),
              alpha = .3,
              linetype = 0) +
  geom_line(aes(x = intervals, y = survival),
            linetype=1,
            size=2)+
  geom_line(aes(x = intervals, y = survival_true),
            linetype = 2,
            size = 1,
            color = cbPalette[5])+
  ylim(0,1) +
  ggtitle("Survival - Time") +
  xlab("Time") +
  ylab("Survival Probability") +
  theme_bw()

ggsave("S0_period.pdf", p_s_period)

##################################
###
### Log Hazard - Period axis
###
##################################

te_indx <- grep("period_effect", rownames(fit_sum))

period_effect_mean <- fit_sum[te_indx, 1]
period_effect_lower <- fit_sum[te_indx, 4]
period_effect_upper <- fit_sum[te_indx, 5]

intervals <- 1:nT_period

out_period_effect <- data.frame(intervals,
                                period_effect_mean,
                                period_effect_lower,
                                period_effect_upper,
                                truth = ageperiod_out$period_effect[1:nT_period]
                                )

period_effect_plot <- ggplot(data = out_period_effect, aes(x = intervals)) +
  geom_line(aes(x = intervals, y = period_effect_mean), size = 1) +
  geom_ribbon(aes(ymin = period_effect_lower, ymax = period_effect_upper),
              alpha = .2,
              linetype = 0) +
  geom_line(aes(x = intervals,y = truth),
            size = 1,
            col = "red") +
            ggtitle("Time Effect Posterior") +
            xlab("Time") +
            ylab("Effect Size") +
            theme_bw()

ggsave("period_effect.pdf", period_effect_plot)


###############################################
###
### Save results
###
##############################################

waic <- mcmcout$WAIC

sink("results_type3_CAR.txt")
print(fit_sum[c(grep("mu", rownames(fit_sum)),
                grep("tau", rownames(fit_sum)),
                grep("beta", rownames(fit_sum))), ])
cat("beta0_true: ", beta0_true, "\n")
print(gelman.diag(out[, c(grep("mu", rownames(fit_sum)),
                          grep("tau", rownames(fit_sum)))]))
print(effectiveSize(out[, c(grep("mu", rownames(fit_sum)),
                           grep("tau", rownames(fit_sum)))]))
cat("runtime", runtime, "\n")
cat("reps: ", reps, "\n")
cat("waic: ", waic$lppd)
sink()

###############################################
###
### Save results objects
###
###############################################

save(waic, file = "waic.Rdata")

ess <- effectiveSize(out[,c(grep("beta", rownames(fit_sum)),
                           grep("sd", rownames(fit_sum)),
                           grep("tau", rownames(fit_sum)))])
save(ess, file = "ess.Rdata")
gd <- gelman.diag(out[,c(grep("beta", rownames(fit_sum)),
                        grep("sd", rownames(fit_sum)),
                        grep("tau", rownames(fit_sum)))])
save(gd, file = "gd.Rdata")
save(fit_sum, file = "fit_sum.Rdata")
save(reps, file = "reps.Rdata")
save(nT_age, file = "nknots_age.Rdata")
save(nT_period, file = "nknots_period.Rdata")