###
### Results - plotting and summary statistics
###

# load("mcmcout_Rdata")
# load("runtime.Rdata")

out <- mcmcout$samples
fit_sum <- mcmcout$summary$all.chains

pdf("traceplots.pdf")
traceplot(out[, "beta0"], ylab = "beta0")
traceplot(out[, "tau_age"], ylab = "tau_age")
traceplot(out[, "tau_period"], ylab = "tau_period")

traceplot(out[, "b_age[1]"], ylab = "b_age[1]")
traceplot(out[, "b_age[3]"], ylab = "b_age[3]")

traceplot(out[, "b_period[1]"], ylab = "b_period[1]")
traceplot(out[, "b_period[3]"], ylab = "b_period[3]")

traceplot(out[, "age_effect[1]"], ylab = "age_effect[1]")
traceplot(out[, "age_effect[3]"], ylab = "age_effect[3]")
traceplot(out[, "age_effect[10]"], ylab = "age_effect[10]")
traceplot(out[, "age_effect[51]"], ylab = "age_effect[51]")

traceplot(out[, "period_effect[1]"], ylab = "period_effect[1]")
traceplot(out[, "period_effect[3]"], ylab = "period_effect[3]")
traceplot(out[, "period_effect[10]"], ylab = "period_effect[10]")
traceplot(out[, "period_effect[51]"], ylab = "period_effect[51]")
dev.off()


###################################################################
###
### print plot S0_age
###
####################################################################

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
Weeks <- 1:study_end

survival <- S0_age[(1):(study_end), 1]
survival_lower <- S0_age[(1):(study_end), 4]
survival_upper <- S0_age[(1):(study_end), 5]

out_s_age <- data.frame(Weeks,
                        survival,
                        survival_lower,
                        survival_upper)

p_s_age <- ggplot(data = out_s_age,aes(x = Weeks)) +
  geom_ribbon(aes(ymin = survival_lower,
                  ymax = survival_upper),
              alpha = .3,
              linetype = 0) +
  geom_line(aes(x = Weeks,
                y = survival),
            linetype = 1,
            size = 2) +
  ylim(0,1)+
  ggtitle("Survival - Age") +
  xlab("Age") +
  ylab("Survival Probability") +
  theme_bw()

p_s_age
ggsave("S0_age.pdf", p_s_age)

##################################
###
### Log hazard - Age Effects
###
##################################

te_indx <- grep("age_effect", rownames(fit_sum))
age_effect_mean <- fit_sum[te_indx, 1]
age_effect_lower <- fit_sum[te_indx, 4]
age_effect_upper <- fit_sum[te_indx, 5]
weeks <- 1:nT_age
out_age_effect <- data.frame(weeks,
                            age_effect_mean,
                            age_effect_lower,
                            age_effect_upper)

age_effect_plot=ggplot(data =out_age_effect,aes(x = weeks))+
  geom_line(aes(x = weeks,y=age_effect_mean),size=1)+
  geom_ribbon(aes(ymin = age_effect_lower,
                  ymax = age_effect_upper),
              alpha = .2,
              linetype = 0) +
  ggtitle("Age Effect Posterior") +
  xlab("Age") + ylab("Effect Size") +
  theme_bw() +
  scale_color_brewer("Year",palette = 'Set1') +
  scale_fill_brewer("Year",palette = 'Set1') +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

age_effect_plot

ggsave("age_effect.pdf",age_effect_plot)

##################################
###
### survival - time axis
###
##################################

study_end <- nT_period

S0_period_indx <- grep("S0_period", rownames(fit_sum))
S0_period <- fit_sum[S0_period_indx, ]
Weeks <- 1:study_end

survival <- S0_period[1:study_end, 1]
survival_lower <- S0_period[1:study_end, 4]
survival_upper <- S0_period[1:study_end, 5]

out_s_period <- data.frame(Weeks,
                           survival,
                           survival_lower,
                           survival_upper)

p_s_period <- ggplot(data = out_s_period, aes(x = Weeks)) +
  geom_ribbon(aes(ymin = survival_lower,
                  ymax = survival_upper),
              alpha = .3,
              linetype = 0) +
  geom_line(aes(x = Weeks,
                y = survival),
            linetype = 1,
            size = 2) +
  ylim(0,1) +
  ggtitle("Survival - Period") +
  xlab("Period") + ylab("Survival Probability") +
  theme_bw()

p_s_period
ggsave("S0_period.pdf",p_s_period)

##################################
###
### Log Hazard - Period Effects
###
##################################

te_indx <- grep("period_effect", rownames(fit_sum))

period_effect_mean <- fit_sum[te_indx, 1]
period_effect_lower <- fit_sum[te_indx, 4]
period_effect_upper <- fit_sum[te_indx, 5]
weeks <- 1:nT_period
out_period_effect <- data.frame(weeks,
                               period_effect_mean,
                               period_effect_lower,
                               period_effect_upper)

period_effect_plot <- ggplot(data =out_period_effect, aes(x = weeks)) +
  geom_line(aes(x = weeks,
                y = period_effect_mean),
            size = 1) +
  geom_ribbon(aes(ymin = period_effect_lower,
                  ymax = period_effect_upper),
              alpha = .2,
              linetype = 0) +
  ggtitle("Period Effect Posterior") +
  xlab("Period") + ylab("Effect Size") +
  theme_bw() +
  scale_color_brewer("Year", palette = 'Set1') +
  scale_fill_brewer("Year", palette = 'Set1') +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

period_effect_plot

ggsave("period_effect.pdf", period_effect_plot)

###############################################
###
### Save results
###
##############################################

waic <- mcmcout$WAIC
save(waic,file = "waic.Rdata")
ess <- effectiveSize(out[,c(grep("beta0",rownames(fit_sum)),
                            grep("sd",rownames(fit_sum)),
                            grep("tau",rownames(fit_sum)))])
save(ess,file="ess.Rdata")
gd <- gelman.diag(out[,c(grep("beta0",rownames(fit_sum)),
                         grep("sd",rownames(fit_sum)),
                         grep("tau",rownames(fit_sum)))])
save(gd, file = "gd.Rdata")
save(fit_sum, file = "fit_sum.Rdata")
save(reps, file = "reps.Rdata")

###############################################
###
### Print results
###
##############################################

sink("results_wtd_C_S.txt")
print(fit_sum[c(grep("beta0",rownames(fit_sum)),grep("tau",rownames(fit_sum)),grep("sd",rownames(fit_sum))),])
print(gd)
print(ess)
cat("Runtime: ",runtime,"\n")
cat("Reps: ",reps,"\n")
print(waic)
sink()
