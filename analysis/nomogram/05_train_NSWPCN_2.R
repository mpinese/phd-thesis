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
library(ggplot2)
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



# Throw in an RSF for fun
library(randomForestSRC)
set.seed(20150111)
fit.rsf = rfsrc(
	Surv(Time, DSD) ~ SexM + AgeCent + LocBody + SizeCent + A2 + A4,
	data = data,
	mtry = 1,
	splitrule = "logrankscore",
	nsplit = 2, 
	ntree = 1000)



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




# Down to just two useful prediction models, gg and gg2.  
# In addition, consider the CPH, despite statistical issues,
# because it's so common.  Select one based on 
# results for the validation set.  What metric to use?
# Use IBS.  Sadly pec doesn't just support plug-in IBS.  DIY
# time.

calcIBS = function(surv, pred, pred_times, max_time)
{
	stopifnot(nrow(surv) == nrow(pred) && length(pred_times) == ncol(pred))

	n = nrow(surv)
	marg_survfit = survfit(surv ~ 1)
	marg_censfit = survfit(Surv(surv[,1], !surv[,2]) ~ 1)
	marg_surv_func = approxfun(marg_survfit$time, marg_survfit$surv, method = "constant", yleft = 1, yright = 0, rule = 2:1, f = 0)
	marg_cens_func = approxfun(marg_censfit$time, marg_censfit$surv, method = "constant", yleft = 1, yright = 0, rule = 2:1, f = 0)

	pred_funcs = apply(pred, 1, function(pat_preds) approxfun(pred_times, pat_preds, yleft = 1, yright = min(pat_preds), rule = 2))

	indiv_patient_bsc = function(pat_i, tstars)
	{
		observed_time = surv[pat_i, 1]
		observed_event = surv[pat_i, 2]
		pred_func = pred_funcs[[pat_i]]
		category = 1*(observed_time <= tstars & observed_event) + 2*(observed_time > tstars) + 3*(observed_time <= tstars & !observed_event)
		bsc = rep(NA, length(tstars))
		bsc[category == 1] = pred_func(tstars[category == 1])^2 / marg_cens_func(observed_time)
		bsc[category == 2] = (1 - pred_func(tstars[category == 2]))^2 / marg_cens_func(tstars[category == 2])
		bsc[category == 3] = 0
		bsc
	}

	bsc_func = function(tstars) { rowMeans(sapply(1:n, function(pat_i) indiv_patient_bsc(pat_i, tstars))) }

	weight_func = function(tstars) { (1 - marg_surv_func(tstars)) / (1 - marg_surv_func(max_time)) }

	# Be slack and do trapezoidal int. with a fine grid.  It should be possible 
	# to calulate the int. exactly but I cbfed.
	int_grid = seq(0, max_time, length.out = 1e3)
	bsc_vals = bsc_func(int_grid)
	weight_vals = weight_func(int_grid)
	int_vals = bsc_vals * weight_vals
	ibsc = (2*sum(int_vals) - int_vals[1] - int_vals[length(int_vals)]) * (diff(range(int_grid))) / (2*length(int_vals))

	return(list(bsc = bsc_vals, weights = weight_vals, eval_times = int_grid, ibsc = ibsc))
}


# Survival predictions for the models
ibs_times = sort(unique(data.val$Time))
ibs_preds_gg = as.matrix(t(sapply(summary(fit.gg, newdata = data.val, type = "survival", t = ibs_times), function(x) x$est)))
ibs_preds_gg2 = as.matrix(t(sapply(summary(fit.gg2, newdata = data.val, type = "survival", t = ibs_times), function(x) x$est)))
temp_cox_preds = survfit(fit.cph, newdata = data.val)
ibs_preds_cph = simplify2array(tapply(1:length(temp_cox_preds$time), rep(names(temp_cox_preds$strata), temp_cox_preds$strata), function(strat_i) { 
	approx(x = temp_cox_preds$time[strat_i], y = temp_cox_preds$surv[strat_i], xout = ibs_times, method = "constant", yleft = 1, rule = 2, f = 0)$y } ))
ibs_preds_cph = t(ibs_preds_cph[,rownames(data.val)])
temp_rsf_preds = predict(fit.rsf, newdata = data.val)
ibs_preds_rsf = t(apply(temp_rsf_preds$survival, 1, function(survs) approx(temp_rsf_preds$time.interest, survs, xout = ibs_times, method = "constant", yleft = 1, rule = 2, f = 0)$y))
# Patients (from data.val) are in rows, times (from ibs_times) in columns.

# Add a no-information KM predictor
temp_km0 = survfit(Surv(Time, DSD) ~ 1, data)
ibs_preds_km0 = t(matrix(rep(approx(temp_km0$time, temp_km0$surv, xout = ibs_times, method = "constant", yleft = 1, rule = 2, f = 0)$y, times = nrow(data.val)), ncol = nrow(data.val)))

# Point estimates.  The normalization constant simplifies comparison between different time scales.
calcIBS(Surv(data.val$Time, data.val$DSD), ibs_preds_gg, ibs_times, max(data.val$Time))$ibs / max(data.val$Time)
calcIBS(Surv(data.val$Time, data.val$DSD), ibs_preds_gg2, ibs_times, max(data.val$Time))$ibs / max(data.val$Time)
calcIBS(Surv(data.val$Time, data.val$DSD), ibs_preds_cph, ibs_times, max(data.val$Time))$ibs / max(data.val$Time)
calcIBS(Surv(data.val$Time, data.val$DSD), ibs_preds_rsf, ibs_times, max(data.val$Time))$ibs / max(data.val$Time)
calcIBS(Surv(data.val$Time, data.val$DSD), ibs_preds_km0, ibs_times, max(data.val$Time))$ibs / max(data.val$Time)

# Bootstrap to get an idea of distributions of the normalized IBSCs
library(plyr)
set.seed(20150111)
bsc_boots = laply(1:5e2, function(i) {
	if (i %% 5e1 == 0)	{ message(i) }
	boot_samp = sample.int(nrow(data.val), replace = TRUE)
	gg = calcIBS(Surv(data.val$Time, data.val$DSD)[boot_samp,], ibs_preds_gg[boot_samp,], ibs_times, max(data.val$Time))$bsc
	gg2 = calcIBS(Surv(data.val$Time, data.val$DSD)[boot_samp,], ibs_preds_gg2[boot_samp,], ibs_times, max(data.val$Time))$bsc
	cph = calcIBS(Surv(data.val$Time, data.val$DSD)[boot_samp,], ibs_preds_cph[boot_samp,], ibs_times, max(data.val$Time))$bsc
	rsf = calcIBS(Surv(data.val$Time, data.val$DSD)[boot_samp,], ibs_preds_rsf[boot_samp,], ibs_times, max(data.val$Time))$bsc
	km0 = calcIBS(Surv(data.val$Time, data.val$DSD)[boot_samp,], ibs_preds_km0[boot_samp,], ibs_times, max(data.val$Time))$bsc
	rbind(gg, gg2, cph, rsf, km0)
})
ibs_eval_times = calcIBS(Surv(data.val$Time, data.val$DSD), ibs_preds_gg, ibs_times, max(data.val$Time))$eval_times
plot(0 ~ 0, type = "n", xlim = c(0, max(data.val$Time)), ylim = c(0, 1))
for (i in 1:dim(bsc_boots)[1])
{
	lines(bsc_boots[i, 1, ] ~ ibs_eval_times, col = rgb(0, 0, 1, 10/dim(bsc_boots)[1]), lwd = 3)
	lines(bsc_boots[i, 2, ] ~ ibs_eval_times, col = rgb(1, 0, 0, 10/dim(bsc_boots)[1]), lwd = 3)
	lines(bsc_boots[i, 3, ] ~ ibs_eval_times, col = rgb(0, 1, 0, 10/dim(bsc_boots)[1]), lwd = 3)
#	lines(bsc_boots[i, 4, ] ~ ibs_eval_times, col = rgb(0, 0, 0, 10/dim(bsc_boots)[1]), lwd = 3)
	lines(bsc_boots[i, 5, ] ~ ibs_eval_times, col = rgb(0, 0, 0, 10/dim(bsc_boots)[1]), lwd = 3)
}
temp = melt(aaply(bsc_boots, 2:3, median))
colnames(temp) = c("Model", "Time", "BSC")
ggplot(temp, aes(x = Time, y = BSC, colour = Model)) + geom_line()

plot(0 ~ 0, type = "n", xlim = c(0, max(data.val$Time)), ylim = c(-0.05, 0.1), xlab = "Time", ylab = "Improvement over KM0")
for (i in 1:dim(bsc_boots)[1])
{
	lines(bsc_boots[i, 5, ] - bsc_boots[i, 1, ] ~ ibs_eval_times, col = rgb(0, 0, 1, 40/dim(bsc_boots)[1]), lwd = 1)
	lines(bsc_boots[i, 5, ] - bsc_boots[i, 2, ] ~ ibs_eval_times, col = rgb(1, 0, 0, 40/dim(bsc_boots)[1]), lwd = 1)
	lines(bsc_boots[i, 5, ] - bsc_boots[i, 3, ] ~ ibs_eval_times, col = rgb(0, 1, 0, 40/dim(bsc_boots)[1]), lwd = 1)
	lines(bsc_boots[i, 5, ] - bsc_boots[i, 4, ] ~ ibs_eval_times, col = rgb(0, 0, 0, 40/dim(bsc_boots)[1]), lwd = 1)
}

temp = sapply(list(gg = ibs_preds_gg, gg2 = ibs_preds_gg2, cph = ibs_preds_cph, rsf = ibs_preds_rsf, km0 = ibs_preds_km0), function(preds) calcIBS(Surv(data.val$Time, data.val$DSD), preds, ibs_times, max(data.val$Time))$bsc)
temp = melt(temp)
colnames(temp) = c("Time", "Model", "BS")
ggplot(temp, aes(x = Time, y = BS, colour = Model)) + geom_line()

temp = melt(aaply(bsc_boots, 2:3, quantile, probs = c(0.05, 0.5, 0.95)))
colnames(temp) = c("Model", "Time", "Quantile", "Value")
temp$Quantile = paste("Q", gsub("%", "", temp$Quantile), sep = "")
temp = dcast(temp, Model + Time ~ Quantile, value.var = "Value")
ggplot(temp, aes(x = Time, y = Q50, ymin = Q5, ymax = Q95, colour = Model, fill = Model)) + geom_line() + geom_ribbon(alpha = 0.2, colour = NA)

bsc_boots_diff = aaply(bsc_boots, 2, function(x) x - bsc_boots[,5,])[1:4,,]
temp = melt(aaply(bsc_boots_diff, c(1,3), quantile, probs = c(0.05, 0.5, 0.95)))
colnames(temp) = c("Model", "Time", "Quantile", "Value")
temp$Quantile = paste("Q", gsub("%", "", temp$Quantile), sep = "")
temp = dcast(temp, Model + Time ~ Quantile, value.var = "Value")
ggplot(temp, aes(x = Time, y = Q50, ymin = Q5, ymax = Q95, colour = Model, fill = Model)) + geom_line() + geom_ribbon(alpha = 0.2, colour = NA)
ggplot(temp, aes(x = Time, y = Q50, ymin = Q5, ymax = Q95, colour = Model, fill = Model)) + geom_line() + geom_ribbon(alpha = 0.2, colour = NA) + xlim(0, 500)
ggplot(temp, aes(x = Time, y = Q50, colour = Model)) + geom_line()

set.seed(20150111)
ibsc_boots = t(sapply(1:5e2, function(i) {
	if (i %% 5e1 == 0)	{ message(i) }
	boot_samp = sample.int(nrow(data.val), replace = TRUE)
	gg = calcIBS(Surv(data.val$Time, data.val$DSD)[boot_samp,], ibs_preds_gg[boot_samp,], ibs_times, max(data.val$Time[boot_samp]))$ibs / max(data.val$Time[boot_samp])
	gg2 = calcIBS(Surv(data.val$Time, data.val$DSD)[boot_samp,], ibs_preds_gg2[boot_samp,], ibs_times, max(data.val$Time[boot_samp]))$ibs / max(data.val$Time[boot_samp])
	cph = calcIBS(Surv(data.val$Time, data.val$DSD)[boot_samp,], ibs_preds_cph[boot_samp,], ibs_times, max(data.val$Time[boot_samp]))$ibs / max(data.val$Time[boot_samp])
	rsf = calcIBS(Surv(data.val$Time, data.val$DSD)[boot_samp,], ibs_preds_rsf[boot_samp,], ibs_times, max(data.val$Time[boot_samp]))$ibs / max(data.val$Time[boot_samp])
	km0 = calcIBS(Surv(data.val$Time, data.val$DSD)[boot_samp,], ibs_preds_km0[boot_samp,], ibs_times, max(data.val$Time[boot_samp]))$ibs / max(data.val$Time[boot_samp])
	c(gg, gg2, cph, rsf, km0)
}))
colnames(ibsc_boots) = c("gg", "gg2", "cph", "rsf", "km0")

boxplot(ibsc_boots)
plot(density(ibsc_boots[,1]), col = "blue", lwd = 2)
lines(density(ibsc_boots[,2]), col = "red", lwd = 2)
lines(density(ibsc_boots[,3]), col = "green", lwd = 2)
lines(density(ibsc_boots[,4]), col = "purple", lwd = 2)
lines(density(ibsc_boots[,5]), col = "black", lwd = 2)

plot(density(ibsc_boots[,5] - ibsc_boots[,1]), col = "blue", lwd = 2, ylim = c(0, 200))
lines(density(ibsc_boots[,5] - ibsc_boots[,2]), col = "red", lwd = 2)
lines(density(ibsc_boots[,5] - ibsc_boots[,3]), col = "green", lwd = 2)
lines(density(ibsc_boots[,5] - ibsc_boots[,4]), col = "purple", lwd = 2)



# All models perform equivalently on the validation set.  Select
# the simplest: gg.

