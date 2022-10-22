###
### Results - plotting and summary statistics
###

# load("mcmcout_Rdata")
# load("runperiod_Rdata")

out <- mcmcout$samples
fit_sum <- mcmcout$summary$all.chains

pdf("traceplots.pdf")
traceplot(out[, "beta0"], ylab = "beta0")
traceplot(out[, "tau_age_cgam"], ylab = "tau_age_cgam")
traceplot(out[, "tau_age_spline"], ylab = "tau_age_spline")
traceplot(out[, "sdk_period"], ylab = "sdk_period")
traceplot(out[, "sda_period"], ylab = "sda_period")

traceplot(out[, "b_age_cgam[1]"], ylab = "b_age_cgam[1]")
traceplot(out[, "b_age_cgam[3]"], ylab = "b_age_cgam[3]")
traceplot(out[, "b_age_spline[1]"], ylab = "b_age_spline[1]")
traceplot(out[, "b_age_spline[3]"], ylab = "b_age_spline[3]")

traceplot(out[, "age_effect[1]"], ylab = "age_effect[1]")
traceplot(out[, "age_effect[3]"], ylab = "age_effect[3]")
traceplot(out[, "age_effect[51]"], ylab = "age_effect[51]")

traceplot(out[, "period_effect_yr1[1]"], ylab = "period_effect_yr1[1]")
traceplot(out[, "period_effect_yr1[3]"], ylab = "period_effect_yr1[3]")
traceplot(out[, "period_effect_yr1[51]"], ylab = "period_effect_yr1[51]")

traceplot(out[, "period_effect_yr2[1]"], ylab = "period_effect_yr2[1]")
traceplot(out[, "period_effect_yr2[3]"], ylab = "period_effect_yr2[3]")
traceplot(out[, "period_effect_yr2[51]"], ylab = "period_effect_yr2[51]")

traceplot(out[, "period_effect_yr3[1]"], ylab = "period_effect_yr3[1]")
traceplot(out[, "period_effect_yr3[3]"], ylab = "period_effect_yr3[3]")
traceplot(out[, "period_effect_yr3[51]"], ylab = "period_effect_yr3[51]")
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
S0_age_indx <- grep("S0_age", rownames(fit_sum))[1:nT_age]
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
  ylim(0, 1)+
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

############
### cgam
############

te_indx <- grep("age_effect_cgam", rownames(fit_sum))
age_effect_mean <- fit_sum[te_indx, 1]
age_effect_lower <- fit_sum[te_indx, 4]
age_effect_upper <- fit_sum[te_indx, 5]
weeks <- 1:nT_age
out_age_effect <- data.frame(weeks,
                            age_effect_mean,
                            age_effect_lower,
                            age_effect_upper)

age_effect_cgam_plot <- ggplot(data = out_age_effect,
                          aes(x = weeks)) +
  geom_line(aes(x = weeks, y = age_effect_mean),
            size = 1) +
  geom_ribbon(aes(ymin = age_effect_lower,
                  ymax = age_effect_upper),
              alpha = .2,
              linetype = 0) +
  ggtitle("Age Effect Posterior CGAM") +
  xlab("Age") + ylab("Effect Size") +
  theme_bw() +
  scale_color_brewer("Year", palette = 'Set1') +
  scale_fill_brewer("Year", palette = 'Set1') +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

age_effect_cgam_plot

ggsave("age_effect_cgam.pdf", age_effect_cgam_plot)

############
### spline
############

te_indx <- grep("age_effect_spline", rownames(fit_sum))
age_effect_mean <- fit_sum[te_indx, 1]
age_effect_lower <- fit_sum[te_indx, 4]
age_effect_upper <- fit_sum[te_indx, 5]
weeks <- 1:nT_age
out_age_effect <- data.frame(weeks,
                            age_effect_mean,
                            age_effect_lower,
                            age_effect_upper)

age_effect_spline_plot <- ggplot(data = out_age_effect,
                          aes(x = weeks)) +
  geom_line(aes(x = weeks, y = age_effect_mean),
            size = 1) +
  geom_ribbon(aes(ymin = age_effect_lower,
                  ymax = age_effect_upper),
              alpha = .2,
              linetype = 0) +
  ggtitle("Age Effect Posterior Spline") +
  xlab("Age") + ylab("Effect Size") +
  theme_bw() +
  scale_color_brewer("Year", palette = 'Set1') +
  scale_fill_brewer("Year", palette = 'Set1') +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

age_effect_spline_plot

ggsave("age_effect_spline.pdf", age_effect_spline_plot)

############
### cgam
############

te_indx <- grep("age_effect", rownames(fit_sum))[1:nT_age]
age_effect_mean <- fit_sum[te_indx, 1]
age_effect_lower <- fit_sum[te_indx, 4]
age_effect_upper <- fit_sum[te_indx, 5]
weeks <- 1:nT_age
out_age_effect <- data.frame(weeks,
                            age_effect_mean,
                            age_effect_lower,
                            age_effect_upper)

age_effect_plot <- ggplot(data = out_age_effect,
                          aes(x = weeks)) +
  geom_line(aes(x = weeks, y = age_effect_mean),
            size = 1) +
  geom_ribbon(aes(ymin = age_effect_lower,
                  ymax = age_effect_upper),
              alpha = .2,
              linetype = 0) +
  ggtitle("Age Effect Posterior") +
  xlab("Age") + ylab("Effect Size") +
  theme_bw() +
  scale_color_brewer("Year", palette = 'Set1') +
  scale_fill_brewer("Year", palette = 'Set1') +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

age_effect_plot

ggsave("age_effect.pdf", age_effect_plot)

##################################
###
### survival - time axis
###
##################################

S0_period_indx <- grep("S0_period", rownames(fit_sum))
S0_period <- fit_sum[S0_period_indx, ]


###
### year 1
###

study_end <- nT_yr1
S0_period_yr1 <- S0_period[1:nT_yr1, ]
Weeks <- 1:study_end

survival_yr1 <- S0_period_yr1[1:study_end, 1]
survival_lower_yr1 <- S0_period_yr1[1:study_end, 4]
survival_upper_yr1 <- S0_period_yr1[1:study_end, 5]

out_s_period_yr1 <- data.frame(Weeks,
                               survival_yr1,
                               survival_lower_yr1,
                               survival_upper_yr1)

p_s_period_yr1 <- ggplot(data =out_s_period_yr1,
                         aes(x = Weeks)) +
  geom_ribbon(aes(ymin = survival_lower_yr1,
                  ymax = survival_upper_yr1),
              alpha = .3,
              linetype = 0) +
  geom_line(aes(x = Weeks,
                y = survival_yr1),
            linetype = 1, size = 2) +
  ylim(0, 1) +
  ggtitle("Survival - Period Year 1")+xlab("Time (Days)") +
  ylab("Survival Probability") +
  theme_bw()

p_s_period_yr1
ggsave("S0_time_yr1.pdf",p_s_period_yr1)

###
### year 2
###
study_end <- nT_yr2
S0_period_yr2 <- S0_period[(nT_period_max+1):(nT_period_max+nT_yr2),]
Weeks <- 1:study_end

survival_yr2 <- S0_period_yr2[1:study_end, 1]
survival_lower_yr2 <- S0_period_yr2[1:study_end, 4]
survival_upper_yr2 <- S0_period_yr2[1:study_end, 5]

out_s_period_yr2=data.frame(Weeks,survival_yr2,survival_lower_yr2,survival_upper_yr2)

p_s_period_yr2 = ggplot(data =out_s_period_yr2,aes(x = Weeks))+
  geom_ribbon(aes(ymin=survival_lower_yr2,ymax=survival_upper_yr2),alpha=.3,linetype=0)+
  geom_line(aes(x = Weeks,y=survival_yr2),linetype=1,size=2)+
  ylim(0, 1)+
  ggtitle("Survival - Period Year 1")+xlab("Time (Days)")+ylab("Survival Probability")+
  theme_bw()

p_s_period_yr2
ggsave("S0_time_yr2.pdf",p_s_period_yr2)

###
### year 3
###

study_end <- nT_yr3
S0_period_yr3 <- S0_period[(2 * nT_period_max + 1):(2 * nT_period_max+nT_yr3), ]
Weeks <- 1:study_end

survival_yr3 <- S0_period_yr3[1:study_end, 1]
survival_lower_yr3 <- S0_period_yr3[1:study_end, 4]
survival_upper_yr3 <- S0_period_yr3[1:study_end, 5]

out_s_period_yr3 <- data.frame(Weeks,survival_yr3,survival_lower_yr3,survival_upper_yr3)

p_s_period_yr3 = ggplot(data =out_s_period_yr3,aes(x = Weeks))+
  geom_ribbon(aes(ymin=survival_lower_yr3,ymax=survival_upper_yr3),alpha=.3,linetype=0)+
  geom_line(aes(x = Weeks,y=survival_yr3),linetype=1,size=2)+
  ylim(0, 1)+
  ggtitle("Survival - Period Year 1")+xlab("Time (Days)")+ylab("Survival Probability")+
  theme_bw()

p_s_period_yr3
ggsave("S0_time_yr3.pdf",p_s_period_yr3)

##################################
###
### log hazard - time axis
###
##################################

te_indx_yr1 <- grep("period_effect_yr1", rownames(fit_sum))
period_effect_mean <- fit_sum[te_indx_yr1, 1]
period_effect_lower <- fit_sum[te_indx_yr1, 4]
period_effect_upper <- fit_sum[te_indx_yr1, 5]
weeks <- 1:(nT_yr1)

out_period_effect <- data.frame(weeks,
                                period_effect_mean,
                                period_effect_lower,
                                period_effect_upper)

period_effect_plot_yr1 <- ggplot(data = out_period_effect,
                                 aes(x = weeks)) +
  geom_line(aes(x = weeks, y = period_effect_mean),
            size = 1) +
  geom_ribbon(aes(ymin = period_effect_lower,
                  ymax = period_effect_upper),
              alpha = .2,
              linetype = 0) +
  ggtitle("Period Effect Posterior Year 1") +
  xlab("Period") + ylab("Effect Size") +
  theme_bw() +
  scale_color_brewer("Year",palette = 'Set1') +
  scale_fill_brewer("Year",palette = 'Set1') +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

period_effect_plot_yr1

ggsave("period_effect_yr1.pdf",period_effect_plot_yr1)

##############
### Year 2
##############

te_indx_yr2 <- grep("period_effect_yr2", rownames(fit_sum))
period_effect_mean <- fit_sum[te_indx_yr2, 1]
period_effect_lower <- fit_sum[te_indx_yr2, 4]
period_effect_upper <- fit_sum[te_indx_yr2, 5]
weeks <- 1:(nT_yr2)

out_period_effect <- data.frame(weeks,period_effect_mean,period_effect_lower,period_effect_upper)

period_effect_plot_yr2 <- ggplot(data = out_period_effect,
                                 aes(x = weeks)) +
  geom_line(aes(x = weeks,
                y = period_effect_mean),
            size = 1) +
  geom_ribbon(aes(ymin = period_effect_lower,
                  ymax = period_effect_upper),
              alpha = .2,
              linetype = 0) +
  ggtitle("Period Effect Posterior Year 1") +
  xlab("Time") + ylab("Effect Size") +
  theme_bw() +
  scale_color_brewer("Year", palette = 'Set1') +
  scale_fill_brewer("Year", palette = 'Set1') +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

period_effect_plot_yr2
ggsave("period_effect_yr2.pdf",period_effect_plot_yr2)

##############
### Year 3
##############

te_indx_yr3 <- grep("period_effect_yr3", rownames(fit_sum))
period_effect_mean  <- fit_sum[te_indx_yr3, 1]
period_effect_lower <- fit_sum[te_indx_yr3, 4]
period_effect_upper <- fit_sum[te_indx_yr3, 5]
weeks <- 1:(nT_yr3)
out_period_effect <- data.frame(weeks,
                                period_effect_mean,
                                period_effect_lower,
                                period_effect_upper)

period_effect_plot_yr3 <- ggplot(data = out_period_effect,
                                 aes(x = weeks)) +
  geom_line(aes(x = weeks,
                y=period_effect_mean),
            size = 1) +
  geom_ribbon(aes(ymin = period_effect_lower,
                  ymax = period_effect_upper),
              alpha = .2,
              linetype = 0) +
  ggtitle("Time Effect Posterior Year 1") +
  xlab("Time") + ylab("Effect Size") +
  theme_bw() +
  scale_color_brewer("Year",palette = 'Set1') +
  scale_fill_brewer("Year",palette = 'Set1') +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

period_effect_plot_yr3
ggsave("period_effect_yr3.pdf",period_effect_plot_yr3)

##################################
###
### survival - ageperiod axis
###
##################################

S0_ageperiod_indx <- grep("S0_ageperiod", rownames(fit_sum))
S0_ageperiod <- fit_sum[S0_ageperiod_indx,]

###
### year 1
###

study_end <- nT_yr1
S0_ageperiod_yr1 <- S0_ageperiod[1:nT_yr1, ]
Weeks <- 1:study_end

survival_yr1 <- S0_ageperiod_yr1[1:study_end, 1]
survival_lower_yr1 <- S0_ageperiod_yr1[1:study_end, 4]
survival_upper_yr1 <- S0_ageperiod_yr1[1:study_end, 5]
out_s_ageperiod_yr1 <- data.frame(Weeks,
                                  survival_yr1,
                                  survival_lower_yr1,
                                  survival_upper_yr1)
p_s_ageperiod_yr1 = ggplot(data = out_s_ageperiod_yr1,
                           aes(x = Weeks)) +
  geom_ribbon(aes(ymin = survival_lower_yr1,
                  ymax = survival_upper_yr1),
              alpha = .3,
              linetype = 0) +
  geom_line(aes(x = Weeks,y=survival_yr1),
            linetype = 1,
            size = 2) +
  ylim(0, 1)+
  ggtitle("Survival - Age+Period Year 1") +
  xlab("Period (Days)") + ylab("Survival Probability") +
  theme_bw()

p_s_ageperiod_yr1
ggsave("S0_ageperiod_yr1.pdf", p_s_ageperiod_yr1)

###
### year 2
###

study_end <- nT_yr2
S0_ageperiod_yr2 <- S0_ageperiod[(nT_period_max+1):(nT_period_max+nT_yr2),]
Weeks <- 1:study_end

survival_yr2 <- S0_ageperiod_yr2[1:study_end, 1]
survival_lower_yr2 <- S0_ageperiod_yr2[1:study_end, 4]
survival_upper_yr2 <- S0_ageperiod_yr2[1:study_end, 5]

out_s_ageperiod_yr2 <- data.frame(Weeks,
                                  survival_yr2,
                                  survival_lower_yr2,
                                  survival_upper_yr2)
p_s_ageperiod_yr2 <- ggplot(data = out_s_ageperiod_yr2,
                            aes(x = Weeks)) +
  geom_ribbon(aes(ymin = survival_lower_yr2,
                  ymax = survival_upper_yr2),
              alpha = .3,
              linetype = 0) +
  geom_line(aes(x = Weeks,y=survival_yr2),
            linetype = 1,
            size = 2) +
  ylim(0, 1) +
  ggtitle("Survival Age+Period Year 2") +
  xlab("Ageperiod (Days)") + ylab("Survival Probability") +
  theme_bw()

p_s_ageperiod_yr2
ggsave("S0_ageperiod_yr2.pdf",p_s_ageperiod_yr2)

###
### year 3
###

study_end <- nT_yr3
S0_ageperiod_yr3 <- S0_ageperiod[(2 * nT_period_max + 1):(2 * nT_period_max + nT_yr3),]
Weeks <- 1:study_end

survival_yr3 <- S0_ageperiod_yr3[1:study_end, 1]
survival_lower_yr3 <- S0_ageperiod_yr3[1:study_end, 4]
survival_upper_yr3 <- S0_ageperiod_yr3[1:study_end, 5]

out_s_ageperiod_yr3 <- data.frame(Weeks,
                                  survival_yr3,
                                  survival_lower_yr3,
                                  survival_upper_yr3)

p_s_ageperiod_yr3 <- ggplot(data = out_s_ageperiod_yr3,
                            aes(x = Weeks)) +
  geom_ribbon(aes(ymin = survival_lower_yr3,
                  ymax = survival_upper_yr3),
              alpha = .3,
              linetype = 0) +
  geom_line(aes(x = Weeks,y=survival_yr3),
            linetype = 1,
            size = 2) +
  ylim(0, 1) +
  ggtitle("Survival Age+Period Year 3") +
  xlab("ageperiod (Days)") + ylab("Survival Probability")+
  theme_bw()

p_s_ageperiod_yr3
ggsave("S0_ageperiod_yr3.pdf",p_s_ageperiod_yr3)

###############################################
###
### Save results
###
##############################################

waic <- mcmcout$WAIC
save(waic,file = "waic.Rdata")
ess <- effectiveSize(out[,c(grep("beta", rownames(fit_sum)),
                            grep("sd", rownames(fit_sum)),
                            grep("tau", rownames(fit_sum)))])
save(ess,file="ess.Rdata")
gd <- gelman.diag(out[,c(grep("beta", rownames(fit_sum)),
                         grep("sd", rownames(fit_sum)),
                         grep("tau", rownames(fit_sum)))],
                  multivariate = FALSE)
save(gd, file = "gd.Rdata")
save(fit_sum, file = "fit_sum.Rdata")
save(reps, file = "reps.Rdata")



###############################################
###
### Print results
###
##############################################

sink("results_grouse_CS_K.txt")
print(fit_sum[c(grep("beta", rownames(fit_sum)),
                grep("sd", rownames(fit_sum)),
                grep("tau", rownames(fit_sum))), ])
print(gd)
print(ess)
cat("Runtime: ", runtime, "\n")
cat("Reps: ", reps, "\n")
print(mcmcout$WAIC)
sink()

