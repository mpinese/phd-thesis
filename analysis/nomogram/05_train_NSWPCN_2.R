load("03_NSWPCN_subset.rda")

library(survival)

x = data[,c("Patient.Sex", "History.Diagnosis.AgeAt.Cent", "Path.LocationBody", "Path.Size.Cent", "Path.Ca199.Preop", "Molec.S100A2.DCThresh", "Molec.S100A4.DCThresh")]
colnames(x) = c("SexM", "AgeCent", "LocBody", "SizeCent", "Ca199", "A2", "A4")
x$SexM = x$Sex == "M"
x$Ca199 = x$Ca199 > 100

x.treat = data[,"Treat.MarginPositive", drop = FALSE]

y = Surv(as.numeric(data$History.Death.Date - data$History.Diagnosis.Date), data$History.DSDeath.Event)

# Note no surgery dates, though for almost all pts there were only a few days difference.

temp = NA
temp = ls()
rm(list = temp[!(temp %in% c("x", "y"))])


sel = !is.na(y[,1]) & !is.na(y[,2]) & !is.na(x$A2) & !is.na(x$A4) & !is.na(x$LocBody)
x = x[sel,]
y = y[sel,]
rm(sel)


# Remove CA-19-9 measurements as they're mostly missing
x = x[,colnames(x) != "Ca199"]


data = as.data.frame(cbind(Time = y[,1], DSD = y[,2], x))
rm(x, y)
data$DSD = data$DSD == 1





# Investigate functional form with the Cox model.
fit.cph.NoAge = coxph(Surv(Time, DSD) ~ SexM + LocBody + SizeCent + A2 + A4, data = data)
fit.cph.NoSize = coxph(Surv(Time, DSD) ~ SexM + AgeCent + LocBody + A2 + A4, data = data)

scatter.smooth(data$AgeCent, resid(fit.cph.NoAge, type = "martingale"), main = "Functional form: Age (Centered)", xlab = "AgeCent", ylab = "Martingale residual")
scatter.smooth(data$SizeCent, resid(fit.cph.NoSize, type = "martingale"), main = "Functional form: Size (Centered)", xlab = "SizeCent", ylab = "Martingale residual")


# It looks like age has a minor nonlinear component, 
# leading to a quadratic-like U shape.  The size relationship
# appears to have a knee, close to size == 0, around which
# the relationship is approximately linear.

# Model age as:  AgeCent + AgeCent^2
# Model size as: SizeCent + SizeCent*I(SizeCent < 0)

data$SizeSmall = data$SizeCent * (data$SizeCent < 0)


# Check for violations of CPH
fit.cph = coxph(Surv(Time, DSD) ~ SexM + AgeCent + I(AgeCent^2) + LocBody + SizeCent + SizeSmall + A2 + A4, data = data)
cox.zph(fit.cph)

# Looks like there's a violation of CPH with gender.  Not unexpected.
# First check whether there is any evidence of gender interaction
anova(coxph(Surv(Time, DSD) ~ SexM*(AgeCent + I(AgeCent^2) + LocBody + SizeCent + SizeSmall + A2 + A4), data = data))
# Nope, good.  We're not interested in gender effects so just stratify.

fit.cph = coxph(Surv(Time, DSD) ~ strata(SexM) + AgeCent + I(AgeCent^2) + LocBody + SizeCent + SizeSmall + A2 + A4, data = data)
cox.zph(fit.cph)

# Looks good.  Slight snifter with age but I'm not particularly concerned.
# Split into age groups and do KM plots to verify
temp.age = cut(data$AgeCent, 4)
temp = survfit(Surv(Time, DSD) ~ temp.age, data)
library(ggplot2)
ggplot(data.frame(surv = temp$surv, time = temp$time, age = rep(names(temp$strata), temp$strata)), aes(y = log(-log(surv)), x = log(time), col = age)) + geom_line()

# Not perfect but it'll do.


# Check for outliers
plot(resid(fit.cph, type = "deviance"))
# Nothing obvious.

# Check for influential observations
temp = resid(fit.cph, type = "dfbetas")
library(reshape2)
colnames(temp) = names(fit.cph$coefficients)
temp = melt(temp)
colnames(temp) = c("Patient", "Coefficient", "dfbetas")
ggplot(temp, aes(y = dfbetas, x = Patient, col = Coefficient)) + geom_point()

# There are quite a few (~ 10) rather influential observations.  These
# could do with some checking, but first collapse down the model -- 
# there's little point doing dfbeta fucking about based on coefficients
# that will never get fit in the end anyway.


library(glmulti)
nobs.coxph <<- function(obj, ...) sum(obj$y[,2])
# Note: Exhaustive search at level 2 is only feasible for at most 5 variables
#fit.cph.as = glmulti(Surv(Time, DSD) ~ strata(SexM) + AgeCent + I(AgeCent^2) + LocBody + SizeCent + SizeSmall + A2 + A4, data = data, marginality = TRUE, method = "h", fitfunction = "coxph", crit = "bic", level = 2)
fit.cph.as = glmulti(Surv(Time, DSD) ~ strata(SexM) + AgeCent + I(AgeCent^2) + LocBody + SizeCent + SizeSmall + A2 + A4, data = data, marginality = TRUE, method = "g", fitfunction = "coxph", crit = "bic", level = 2)
# After 990 generations:
# Best model: Surv(Time,DSD)~1+strata(SexM)+SizeCent+A2+A4
# Crit= 1796.11096470108
# Mean crit= 1825.07531885092
# Improvements in best and average IC have bebingo en below the specified goals.
# Algorithm is declared to have converged.
# Completed.

rm(nobs.coxph)
library(MASS)
stepAIC(fit.cph, k = log(nrow(data)))

# Consensus, excellent.


# Now generate the restricted fit and examine the dfbetas.
fit.cph = coxph(Surv(Time, DSD) ~ strata(SexM) + SizeCent + A2 + A4, data = data)
cox.zph(fit.cph)	# Still good
plot(resid(fit.cph, type = "deviance"))
temp = resid(fit.cph, type = "dfbetas")
colnames(temp) = names(fit.cph$coefficients)
temp = melt(temp)
colnames(temp) = c("Patient", "Coefficient", "dfbetas")
ggplot(temp, aes(y = abs(dfbetas), x = Patient, col = Coefficient)) + geom_point()

# Tighter, as expected.  11 points have |dfbetas| > 0.2, 20 > 0.15, and 48 > 0.1.

sort(apply(abs(resid(fit.cph, type = "dfbetas")), 1, max), decreasing = TRUE)



# So where are we now?
#   1. On the basis of pre-operative assessability and data availability, variables 
#      were filtered down to Sex, AgeCent, LocBody, SizeCent, A2, A4.
#   2. Functional forms for the continuous variates AgeCent and SizeCent indicated
#      a possible slight quadratic effect on AgeCent, and a knee on SizeCent.  These
#      were modelled by incorporating additional terms.
#   3. Analysis of a full model fit (with additional nonlinear terms included) indicated
#      violation of PH for gender.  This was dealt with by stratification.  A slight
#      PH violation by age was deemed unimportant.
#   4. Variable selection by BIC (both stepwise and genetic all-subset) settled on a 
#      final model of Surv(Time,DSD)~1+strata(SexM)+SizeCent+A2+A4.  This model
#      was refit by coxph.
#   5. PH was verified on the final model.  Deviance residuals showed no egregious outliers.
#      dfBetaS indicated a number of influential observations, which require checking.











# Start with a demonstrative full fit with a flexible generalized gamma model
library(flexsurv)

fit.gg = flexsurvreg(Surv(Time, DSD) ~ SexM + SizeCent + A2 + A4,
	anc = list(
		sigma = ~ SexM,
		Q = ~ SexM),
	data = data, dist = "gengamma")

AIC(fit.gg)
BIC(fit.gg)


# It'll be tricky to assess the goodness of fit of this model.  Best just roll
# with it and leave that to the validation cohorts.

plot(fit.gg)
plot(fit.gg, type = "hazard")
plot(fit.gg, type = "cumhaz")

# Still, fit hardly looks great.  Difficult to judge though due to the 
# effective stratification by Sex and the averaging performed for the 
# plots -- there are no data points at the half-male half-female individual
# being plotted here, so it's not surprising that the fit isn't perfect.

# Better to plot male and female curves at various combinations of the
# other factors.  Or even just median.

temp.data = expand.grid(A4 = c(FALSE, TRUE), A2 = c(FALSE, TRUE), SexM = c(FALSE, TRUE), SizeCent = 0)
temp.preds = summary(fit.gg, newdata = temp.data, type = "survival", t = seq(0, 365*5, 30))
temp.survfit = survfit(Surv(Time, DSD) ~ SexM + A2 + A4, data)
plot(temp.survfit, col = 1:length(temp.preds), fun = "cloglog")
for (i in 1:length(temp.preds))
{
	lines(temp.preds[[i]]$time, log(-log(temp.preds[[i]]$est)), col = i)
}

#plot(temp.survfit)
plot(temp.survfit, col = 1:length(temp.preds))
for (i in 1:length(temp.preds))
{

	lines(temp.preds[[i]]$time, temp.preds[[i]]$est, col = i)
}
