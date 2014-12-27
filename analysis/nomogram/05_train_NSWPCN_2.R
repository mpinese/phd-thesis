load("03_NSWPCN_subset.rda")

pdf("05_train_NSWPCN_2.pdf")

options(echo = TRUE, warn = 1)

library(survival)
library(rms)
library(randomForestSRC)

data.all = cbind(time = data.y[,1], event = data.y[,2], data.x.all, data.x.extra)
dd <<- datadist(data.all)
options(datadist = "dd")

#####################################################################
## RANDOM FOREST 

set.seed(1234)
rsf.molec_preop = rfsrc(
	Surv(time, event) ~ 
		Patient.Sex + 
		History.Diagnosis.AgeAt.Cent + 
		Path.LocationBody + 
		Path.Size.Cent + 
		Molec.S100A4.DCThresh + 
		Molec.S100A2.DCThresh, 
	data = data.all,
	mtry = 1,
	splitrule = "logrankscore",
	nsplit = 2, 
	ntree = 1000)


#####################################################################
## KAPLAN MEIER

km.molec_preop = survfit(
	Surv(time, event) ~ 
		Molec.S100A4.DCThresh + 
		Molec.S100A2.DCThresh, 
	data = data.all)
km.molec_preop.strata = rep(gsub(" ", "", names(km.molec_preop$strata)), km.molec_preop$strata)


#####################################################################
## WEIBULL TESTING

temp.start = 1
for (temp.i in 1:length(km.molec_preop$strata))
{
	temp.times = km.molec_preop$time[temp.start:(temp.start+km.molec_preop$strata[temp.i])]
	temp.surv = km.molec_preop$surv[temp.start:(temp.start+km.molec_preop$strata[temp.i])]
	if (temp.i == 1) { 
		plot(log(temp.times), log(-log(temp.surv)), type = "s", col = temp.i, xlab = "log(t)", ylab = "log(-log(hat(S)(t)))", ylim = c(-5, 2)) 
	} else {
		lines(log(temp.times), log(-log(temp.surv)), type = "s", col = temp.i) 
	}
	temp.start = temp.start + km.molec_preop$strata[temp.i] + 1
}
# Not great fits to Weibull.  :/


#####################################################################
## WEIBULL

# weibull.molec_preop = psm(
# 	Surv(time, event) ~ 
# 		Patient.Sex + 
# 		History.Diagnosis.AgeAt.Cent + 
# 		Path.LocationBody + 
# 		Path.Size.Cent + 
# 		Molec.S100A4.DCThresh + 
# 		Molec.S100A2.DCThresh, 
# 	data = data.all, dist = "weibull", x = TRUE, y = TRUE)
# anova(weibull.molec_preop)
# fastbw(weibull.molec_preop)
# validate(weibull.molec_preop, B = 1000)


#####################################################################
## COX PROPORTIONAL HAZARD

temp.cox = cph(
	Surv(time, event) ~ 1,
	data = data.all, x = TRUE, y = TRUE)
for (i in 1:ncol(data.all)) { 
	scatter.smooth(data.all[,i], resid(temp.cox, type = "martingale"), xlab = colnames(data.all)[i], ylab = "Martingale residual", main = sprintf("Cox form: %s", colnames(data.all)[i]))
	abline(h = 0) 
}
# Looks like a threshold effect for size.  Introduce a new
# variable to code this.
data.all$Path.Size.Cent.Thresh = pmin(data.all$Path.Size.Cent, 0)
dd <<- datadist(data.all)
options(datadist = "dd")

temp.cox = cph(
	Surv(time, event) ~ 
		Patient.Sex + 
		History.Diagnosis.AgeAt.Cent + 
		Path.LocationBody + 
		Path.Size.Cent.Thresh + 
		Molec.S100A4.DCThresh + 
		Molec.S100A2.DCThresh, 
	data = data.all, x = TRUE, y = TRUE, surv = FALSE)
cox.zph(temp.cox)
plot(cox.zph(temp.cox))

# zph analysis indicated violations for Patient.Sex and Molec.S100A4.DCThresh.
# Stratify on these to get around the problem
cox.molec_preop = cph(
	Surv(time, event) ~ 
		strat(Patient.Sex) + 
		History.Diagnosis.AgeAt.Cent + 
		Path.LocationBody + 
		Path.Size.Cent.Thresh + 
		strat(Molec.S100A4.DCThresh) + 
		Molec.S100A2.DCThresh, 
	data = data.all, x = TRUE, y = TRUE, surv = FALSE)
cox.zph(cox.molec_preop)			# All good now
plot(survfit(cox.molec_preop), col = 1:4, lwd = 2)

library(survivalMPL)
for (i in 1:4)
{
	temp.data = data.all[xor(data.all$Patient.Sex == "M", i %% 2) & xor(data.all$Molec.S100A4.DCThresh, i <= 2),]
	temp.mpl = coxph_mpl(
	Surv(time, event) ~ 
		History.Diagnosis.AgeAt.Cent + 
		Path.LocationBody + 
		Path.Size.Cent.Thresh + 
		Molec.S100A2.DCThresh, 
	data = temp.data, control = coxph_mpl.control(n.obs = sum(temp.data$event), basis = "gaussian"))
	plot(temp.mpl)
}

temp.dfbetas <- residuals(cox.molec_preop, type = "dfbetas")
for (j in 1:ncol(temp.dfbetas)) {
	plot(temp.dfbetas[,j], ylab = names(coef(cox.molec_preop))[j], ylim = c(-max(c(abs(temp.dfbetas[,j]), 1)), max(c(abs(temp.dfbetas[,j]), 1))))
	abline(h = c(-1, 0, 1), col = c("red", "grey", "red"))
}
# Some outliers that will need checking if they persist.  First
# wait for the updated data from DC -- maybe they'll be resolved.

plot(survfit(cox.molec_preop))
temp.haz = basehaz(cox.molec_preop, centered = TRUE)
temp.strata = sort(unique(temp.haz$strata))
plot(exp(-hazard) ~ time, temp.haz, type = "n", xlab = "Time", ylab = "Survival", ylim = c(0, 1))
for (i in 1:length(temp.strata))
{
	lines(temp.haz$time[temp.haz$strata == temp.strata[i]], exp(-temp.haz$hazard[temp.haz$strata == temp.strata[i]]), type = "s", col = i, lwd = 2)
}
legend("topright", legend = temp.strata, lty = "solid", col = 1:length(temp.strata), lwd = 2, inset = 0.05)

temp.haz = basehaz(cox.molec_preop, centered = FALSE)
temp.strata = sort(unique(temp.haz$strata))
plot(hazard ~ time, temp.haz, type = "n", xlab = "Time", ylab = "Cumulative Hazard", main = "Base hazard estimate")
for (i in 1:length(temp.strata))
{
	lines(temp.haz$time[temp.haz$strata == temp.strata[i]], temp.haz$hazard[temp.haz$strata == temp.strata[i]], type = "s", col = i, lwd = 2)
}
legend("bottomright", legend = temp.strata, lty = "solid", col = 1:length(temp.strata), lwd = 2, inset = 0.05)




# Weibull / Exponential
plot(log(hazard) ~ log(time), temp.haz, type = "n", xlab = "log(Time)", ylab = "log(Cumulative Hazard)", main = "Base hazard estimate: Weibull / Exponential Transform")
for (i in 1:length(temp.strata))
{
	lines(log(temp.haz$time[temp.haz$strata == temp.strata[i]]), log(temp.haz$hazard[temp.haz$strata == temp.strata[i]]), type = "s", col = i, lwd = 2)
}
legend("bottomright", legend = temp.strata, lty = "solid", col = 1:length(temp.strata), lwd = 2, inset = 0.05)

# Gompertz
plot(log(hazard) ~ time, temp.haz, type = "n", xlab = "Time", ylab = "log(Cumulative Hazard)", main = "Base hazard estimate: Gompertz Transform")
for (i in 1:length(temp.strata))
{
	lines(temp.haz$time[temp.haz$strata == temp.strata[i]], log(temp.haz$hazard[temp.haz$strata == temp.strata[i]]), type = "s", col = i, lwd = 2)
}
legend("bottomright", legend = temp.strata, lty = "solid", col = 1:length(temp.strata), lwd = 2, inset = 0.05)

# Log-logistic
logit = function(x) { log(x) - log(1-x) }
plot(logit(exp(-hazard)) ~ log(time), temp.haz, type = "n", xlab = "log(Time)", ylab = "logit(exp(-Cumulative Hazard))", main = "Base hazard estimate: Log-Logistic Transform")
for (i in 1:length(temp.strata))
{
	lines(log(temp.haz$time[temp.haz$strata == temp.strata[i]]), logit(exp(-temp.haz$hazard[temp.haz$strata == temp.strata[i]])), type = "s", col = i, lwd = 2)
}
legend("bottomright", legend = temp.strata, lty = "solid", col = 1:length(temp.strata), lwd = 2, inset = 0.05)






library(mfp)
plot(hazard ~ time, temp.haz, type = "n", xlab = "Time", ylab = "Cumulative Hazard", main = "Base hazard estimate")
temp.time = seq(0, max(temp.haz$time), 1)
for (i in 1:length(temp.strata))
{
	lines(temp.haz$time[temp.haz$strata == temp.strata[i]], temp.haz$hazard[temp.haz$strata == temp.strata[i]], type = "s", col = i, lwd = 2)
	temp.haz.fit = mfp(log(hazard) ~ fp(time) + 0, data = temp.haz[temp.haz$strata == temp.strata[i],])
	lines(temp.time, exp(predict(temp.haz.fit, newdata = data.frame(time = temp.time))), col = i, lty = "dotted", lwd = 2)
}
#legend("bottomright", legend = temp.strata, lty = "solid", col = 1:length(temp.strata), lwd = 2, inset = 0.05)


library(bshazard)
temp.lp = predict(cox.molec_preop)
for (i in 1:length(temp.strata))
{
	temp.bshaz = bshazard(Surv(time, event) ~ temp.lp[temp.haz$strata == temp.strata[i]], data.all[temp.haz$strata == temp.strata[i],])
	plot(temp.bshaz)
}


save.image("05_NSWPCN_fits_2.rda")
