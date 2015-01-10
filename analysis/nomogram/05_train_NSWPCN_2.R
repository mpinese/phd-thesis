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



# There's going to be an awful lot of model manipulation and black
# magic going on.  Create a holdout validation set for model comparison
# once the chips are down.
set.seed(20150110)
sel.val = sample.int(nrow(data), floor(nrow(data)/5))
sel.val = 1:nrow(data) %in% sel.val
mean(sel.val)
data.val = data[sel.val,,drop = FALSE]
data = data[!sel.val,,drop = FALSE]


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
# temp.age = cut(data$AgeCent, 4)
# temp = survfit(Surv(Time, DSD) ~ temp.age, data)
# library(ggplot2)
# ggplot(data.frame(surv = temp$surv, time = temp$time, age = rep(names(temp$strata), temp$strata)), aes(y = log(-log(surv)), x = log(time), col = age)) + geom_line()
# Not perfect but it'll do.


# Check for outliers
plot(resid(fit.cph, type = "deviance"))
abline(h = c(-2.5, 2.5))
# Nothing obvious.
plot(resid(fit.cph, type = "deviance")[order(data$SexM, data$A2, data$A4)], col = (4*data$SexM + 2*data$A2 + data$A4)[order(data$SexM, data$A2, data$A4)], pch = 16)
abline(h = c(-2.5, 2.5))
boxplot(resid(fit.cph, type = "deviance") ~ data$SexM + data$A2 + data$A4, varwidth = TRUE)
abline(h = 0)
boxplot(resid(fit.cph, type = "martingale") ~ data$SexM + data$A2 + data$A4, varwidth = TRUE)
abline(h = 0)

# Check for influential observations
temp = resid(fit.cph, type = "dfbetas")
library(reshape2)
colnames(temp) = names(fit.cph$coefficients)
temp = melt(temp)
colnames(temp) = c("Patient", "Coefficient", "dfbetas")
ggplot(temp, aes(y = dfbetas, x = Patient, col = Coefficient)) + geom_point()

# There is quite a number of rather influential observations.  These
# could do with some checking, but first collapse down the model -- 
# there's little point doing dfbeta fucking about based on coefficients
# that will never get fit in the end anyway.


library(glmulti)
nobs.coxph <<- function(obj, ...) sum(obj$y[,2])
# Note: Exhaustive search at level 2 is only feasible for at most 5 variables
#fit.cph.as = glmulti(Surv(Time, DSD) ~ strata(SexM) + AgeCent + I(AgeCent^2) + LocBody + SizeCent + SizeSmall + A2 + A4, data = data, marginality = TRUE, method = "h", fitfunction = "coxph", crit = "bic", level = 2)
set.seed(20150110)
fit.cph.as = glmulti(Surv(Time, DSD) ~ strata(SexM) + AgeCent + I(AgeCent^2) + LocBody + SizeCent + SizeSmall + A2 + A4, data = data, marginality = TRUE, method = "g", fitfunction = "coxph", crit = "bic", level = 2)
# After 830 generations:
# Best model: Surv(Time,DSD)~1+strata(SexM)+SizeCent+A2+A4
# Crit= 1367.16344569113
# Mean crit= 1401.37248769175
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

# The classic threshold for concern is 2/sqrt(n).  Unsure about the interpretation
# of n here -- all samples, or just uncensored?  Thankfully similar for NSWPCN.
# Use all n as a more stringent cutoff.
2/sqrt(nrow(data))
ggplot(temp, aes(y = abs(dfbetas), x = Patient, col = Coefficient)) + geom_point() + geom_hline(yintercept = 2/sqrt(nrow(data)))

sort(apply(abs(resid(fit.cph, type = "dfbetas")), 1, max), decreasing = TRUE)
sum(apply(abs(resid(fit.cph, type = "dfbetas")), 1, max) > 2/sqrt(nrow(data)))

# 19 potentially problematic points.




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





library(flexsurv)

fit.gg = flexsurvreg(Surv(Time, DSD) ~ SexM + SizeCent + A2 + A4,
	anc = list(
		sigma = ~ SexM,
		Q = ~ SexM),
	data = data, dist = "gengamma")

fit.gf = flexsurvreg(Surv(Time, DSD) ~ SexM + SizeCent + A2 + A4,
	anc = list(
		sigma = ~ SexM,
		Q = ~ SexM,
		P = ~ SexM),
	data = data, dist = "genf")

fit.gg$loglik
fit.gf$loglik
pchisq(2*(fit.gf$loglik - fit.gg$loglik), 2, lower.tail = FALSE)

fit.gg
AIC(fit.gg)
BIC(fit.gg)

# Plot fit stratified by sex, separate curves for A2, A4 status, at
# median (approx.) Size.

temp.grid = expand.grid(A4 = c(FALSE, TRUE), A2 = c(FALSE, TRUE), SexM = c(FALSE, TRUE), SizeCent = 0)
temp.grid$ID = sprintf("SexM=%s, A2=% -5s, A4=% -5s", temp.grid$SexM, temp.grid$A2, temp.grid$A4)
temp.preds = summary(fit.gg, newdata = temp.grid, type = "survival", t = seq(0, 365*5, 30))
temp.preds2 = do.call(rbind, temp.preds)
temp.preds2$group = rep(gsub(".*ID=", "", names(temp.preds)), each = nrow(temp.preds[[1]]))
temp.preds.cox = survfit(fit.cph, newdata = temp.grid)
temp.survfit = survfit(Surv(Time, DSD) ~ SexM + A2 + A4, data)
temp.data = data.frame(time = temp.survfit$time, surv = temp.survfit$surv, upper = temp.survfit$lower, lower = temp.survfit$upper, group = rep(names(temp.survfit$strata), temp.survfit$strata), model = "KM")
temp.data = rbind(temp.data, data.frame(time = temp.preds2$time, surv = temp.preds2$est, upper = temp.preds2$ucl, lower = temp.preds2$lcl, group = temp.preds2$group, model = "GG"))
temp.data = rbind(temp.data, data.frame(time = temp.preds.cox$time, surv = temp.preds.cox$surv, upper = temp.preds.cox$upper, lower = temp.preds.cox$lower, group = rep(temp.grid$ID, temp.preds.cox$strata), model = "CPH"))
temp.data$Sex = c("Male", "Female")[grepl("SexM=FALSE", temp.data$group)+1]
temp.data$A2 = c("A2-", "A2+")[grepl("A2=TRUE", temp.data$group)+1]
temp.data$A4 = c("A4-", "A4+")[grepl("A4=TRUE", temp.data$group)+1]

ggplot(temp.data, aes(x = log(time), y = log(-log(surv)), ymin = log(-log(lower)), ymax = log(-log(upper)), colour = model, fill = model)) + 
	geom_ribbon(alpha = 0.25, colour = NA) + 
	geom_line() + 
	xlim(4, 7) + ylim(-4, 2) + 
	facet_grid(A2 ~ A4 ~ Sex)
ggplot(temp.data, aes(x = time, y = surv, ymin = lower, ymax = upper, colour = model, fill = model)) + 
	geom_ribbon(alpha = 0.25, colour = NA) + 
	geom_line() + xlim(0, 2000) + ylim(0, 1) + 
	facet_grid(A2 ~ A4 ~ Sex)

# Nice plots.  Some deviation though not significant.  Most concerning
# is the A2- A4- female group, survival of which is underestimated by
# the flexsurv model.  To approach this in a modelling sense would
# require interaction terms between Sex and A2, A4.  
# Overfitting seems likely considering the very few data available for
# the A2+/A4- group.  Perhaps just add a single "DoubleNegFemale" term.

fit.gg2 = flexsurvreg(Surv(Time, DSD) ~ SexM + SizeCent + A2 + A4 + I(SexM == FALSE & A2 == FALSE & A4 == FALSE),
	anc = list(
		sigma = ~ SexM,
		Q = ~ SexM),
	data = data, dist = "gengamma")

AIC(fit.gg)
AIC(fit.gg2)
AIC(fit.gg) - AIC(fit.gg2)
# Equivocal on AIC.  BIC would favour gg then.

pchisq(-2*(fit.gg$loglik - fit.gg2$loglik), 1, lower.tail = FALSE)
# Not good evidence on LRT
temp.preds = summary(fit.gg2, newdata = temp.grid, type = "survival", t = seq(0, 365*5, 30))
temp.preds2 = do.call(rbind, temp.preds)
temp.preds2$group = rep(gsub(".*ID=", "", names(temp.preds)), each = nrow(temp.preds[[1]]))
temp.data = rbind(temp.data, data.frame(time = temp.preds2$time, surv = temp.preds2$est, upper = temp.preds2$ucl, lower = temp.preds2$lcl, group = temp.preds2$group, model = "GG2", Sex = NA, A2 = NA, A4 = NA))
temp.data$Sex = c("Male", "Female")[grepl("SexM=FALSE", temp.data$group)+1]
temp.data$A2 = c("A2-", "A2+")[grepl("A2=TRUE", temp.data$group)+1]
temp.data$A4 = c("A4-", "A4+")[grepl("A4=TRUE", temp.data$group)+1]

ggplot(temp.data, aes(x = log(time), y = log(-log(surv)), ymin = log(-log(lower)), ymax = log(-log(upper)), colour = model, fill = model)) + 
	geom_ribbon(alpha = 0.25, colour = NA) + 
	geom_line() + 
	xlim(4, 7) + ylim(-4, 2) + 
	facet_grid(A2 ~ A4 ~ Sex)
ggplot(temp.data, aes(x = time, y = surv, ymin = lower, ymax = upper, colour = model, fill = model)) + 
	geom_ribbon(alpha = 0.25, colour = NA) + 
	geom_line() + xlim(0, 2000) + ylim(0, 1) + 
	facet_grid(A2 ~ A4 ~ Sex)

# That fixes it pretty well.  Despite the modelling indicators
# favouring the simpler model, my gut feeling is to go with the
# more complex one.  Validation will sort out any tweaking.


temp.data$lower[temp.data$model != "KM"] = NA
temp.data$upper[temp.data$model != "KM"] = NA
ggplot(temp.data, aes(x = log(time), y = log(-log(surv)), ymin = log(-log(lower)), ymax = log(-log(upper)), colour = model, fill = model)) + 
	geom_ribbon(alpha = 0.25, colour = NA) + 
	geom_line() + 
	xlim(4, 7) + ylim(-4, 2) + 
	facet_grid(A2 ~ A4 ~ Sex)
ggplot(temp.data, aes(x = time, y = surv, ymin = lower, ymax = upper, colour = model, fill = model)) + 
	geom_ribbon(alpha = 0.25, colour = NA) + 
	geom_line() + xlim(0, 2000) + ylim(0, 1) + 
	facet_grid(A2 ~ A4 ~ Sex)




# Down to just two models, gg and gg2.  Select one based on 
# results for the validation set.  What metric to use?
library(pec)
pec_times = sort(unique(data.val$Time))
pec_preds_gg = t(sapply(summary(fit.gg, newdata = data.val, type = "survival", t = pec_times), function(x) x$est))
pec_preds_gg2 = t(sapply(summary(fit.gg2, newdata = data.val, type = "survival", t = pec_times), function(x) x$est))
pec(list(gg = pec_preds_gg, gg2 = pec_preds_gg2),
	Surv(Time, DSD) ~ strata(SexM) + SizeCent + A2 + A4, 
	data = data.val,
	times = pec_times,
	exact = TRUE,
	cens.model = "cox")

