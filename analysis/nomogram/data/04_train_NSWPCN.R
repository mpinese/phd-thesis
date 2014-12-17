load("03_NSWPCN_subset.rda")

pdf("04_train_NSWPCN.pdf")


library(survival)


# What's going on with the diagnosis date prognostic?  I suspect something
# else is correlated with it, and that other thing is the true prognostic.
temp = coxph(data.y ~ pspline(History.Diagnosis.Date.Cent), data = data.x.all)
summary(temp)
termplot(temp, se = TRUE, rug = TRUE)
temp = coxph(data.y ~ History.Diagnosis.Date.Cent, data = data.x.all)
summary(temp)
termplot(temp, se = TRUE, rug = TRUE)
# A decrease in hazard as time goes on.  Better surgery?  Chemo?  Cohort
# selection?

# Could it be margins?
scatter.smooth(data.x.all$History.Diagnosis.Date.Cent, data.x.all$Treat.MarginPositive)
# Plausible.
anova(coxph(data.y ~ Treat.MarginPositive + History.Diagnosis.Date.Cent, data.x.all))
# Pretty close, but not everything.

# Earlier diagnosis?
scatter.smooth(data.x.all$History.Diagnosis.Date.Cent, data.x.all$History.Diagnosis.AgeAt.Cent)
# No.
scatter.smooth(data.x.all$History.Diagnosis.Date.Cent, data.x.all$Path.Size.Cent)
# No.

# Surgery on different locations?
scatter.smooth(data.x.all$History.Diagnosis.Date.Cent, data.x.all$Path.LocationBody)
# No.

# Changes in subtype (cohort selection rules)?
anova(lm(History.Diagnosis.Date.Cent ~ Path.Subtype, data.x.all))
# Looks promising
anova(coxph(data.y ~ Treat.MarginPositive + Path.Subtype + History.Diagnosis.Date.Cent, data.x.all))
# Fit issues but I think this could be it.
anova(coxph(data.y ~ Treat.MarginPositive + I(Path.Subtype == "Tubular") + History.Diagnosis.Date.Cent, data.x.all))
# Bam.
anova(coxph(data.y ~ I(Path.Subtype == "Tubular") + Treat.MarginPositive + History.Diagnosis.Date.Cent, data.x.all))

temp1 = sort(data.x.all$History.Diagnosis.Date.Cent)
temp2 = data.x.all$Path.Subtype[order(data.x.all$History.Diagnosis.Date.Cent)]
plot(temp1, cumsum(temp2 == "Tubular") / sum(temp2 == "Tubular"), type = "s", col = "blue", xlab = "Centered Time of Diagnosis", ylab = "Fraction Collected", main = "Time dependency of subtype collection")
lines(temp1, cumsum(temp2 == "Adenosquamous") / sum(temp2 == "Adenosquamous"), type = "s", col = "red")
lines(temp1, cumsum(temp2 == "NotSpecified") / sum(temp2 == "NotSpecified"), type = "s", col = "green")
legend("bottomright", inset = 0.05, legend = c("Tubular", "Adenosquamous", "NotSpecified"), col = c("blue", "green", "red"), lty = "solid")
plot(temp1, cumsum(temp2 == "Tubular"), type = "s", col = "blue", xlab = "Centered Time of Diagnosis", ylab = "Total Collected", main = "Time dependency of subtype collection")
lines(temp1, cumsum(temp2 == "Adenosquamous"), type = "s", col = "red")
lines(temp1, cumsum(temp2 == "NotSpecified"), type = "s", col = "green")
legend("topleft", inset = 0.05, legend = c("Tubular", "Adenosquamous", "NotSpecified"), col = c("blue", "green", "red"), lty = "solid")
# So early cancers were often adenosquamous.  Then rates were approximately stable for
# all three main types, but by around 2003 (centered time == 0), no more adenosquamous cancers
# were collected.
survdiff(data.y ~ data.x.all$Path.Subtype)


# Practically speaking this means that I can remove History.Diagnosis.Date as 
# a confounder, as it is captured by subtype (which in turns appears to be 
# largely dependent on A2/A4).

data.x.molec_preop = data.x.molec_preop[,colnames(data.x.molec_preop) != "History.Diagnosis.Date.Cent"]
data.x.molec_postop = data.x.molec_postop[,colnames(data.x.molec_postop) != "History.Diagnosis.Date.Cent"]
data.x.conv_preop = data.x.conv_preop[,colnames(data.x.conv_preop) != "History.Diagnosis.Date.Cent"]
data.x.conv_postop = data.x.conv_postop[,colnames(data.x.conv_postop) != "History.Diagnosis.Date.Cent"]
data.x.all = data.x.all[,colnames(data.x.all) != "History.Diagnosis.Date.Cent"]


library(rJava)
.jinit(parameters = "-Xmx4g")
library(glmulti)
nobs.coxph <<- function(object, ...) { object$n }


set.seed(1234)
glmulti.molec_preop.bic = glmulti(
	data.y ~ 
		Patient.Sex + 
		History.Diagnosis.AgeAt.Cent + 
#		ns(History.Diagnosis.AgeAt.Cent, 4) + 		# EDA indicates that linear terms are sufficient
		Path.LocationBody + 
		Path.Size.Cent + 
#		ns(Path.Size.Cent, 4) + 					# EDA indicates that linear terms are sufficient
		Stage.pT.Simplified +
		Molec.S100A4.DCThresh + 
		Molec.S100A2.DCThresh, 
	data = data.x.molec_preop,
	level = 2,
	method = "g",
	crit = "bic", 
	includeobjects = TRUE, 
	plotty = FALSE, 
	report = TRUE,
	marginality = TRUE,
	fitfunction = "coxph")
set.seed(1234)
glmulti.molec_preop.aicc = glmulti(
	data.y ~ 
		Patient.Sex + 
		History.Diagnosis.AgeAt.Cent + 
#		ns(History.Diagnosis.AgeAt.Cent, 4) + 		# EDA indicates that linear terms are sufficient
		Path.LocationBody + 
		Path.Size.Cent + 
#		ns(Path.Size.Cent, 4) + 					# EDA indicates that linear terms are sufficient
		Stage.pT.Simplified +
		Molec.S100A4.DCThresh + 
		Molec.S100A2.DCThresh, 
	data = data.x.molec_preop,
	level = 2,
	method = "g",
	crit = "aicc", 
	includeobjects = TRUE, 
	plotty = FALSE, 
	report = TRUE,
	marginality = TRUE,
	fitfunction = "coxph")


set.seed(1234)
glmulti.conv_preop.bic = glmulti(
	data.y ~ 
		Patient.Sex + 
		History.Diagnosis.AgeAt.Cent + 
#		ns(History.Diagnosis.AgeAt.Cent, 4) + 
		Path.LocationBody + 
		Path.Size.Cent + 
#		ns(Path.Size.Cent, 4) + 
		Stage.pT.Simplified,
	data = data.x.conv_preop,
	level = 2,
	method = "g",
	crit = "bic", 
	includeobjects = TRUE, 
	plotty = FALSE, 
	report = TRUE,
	marginality = TRUE,
	fitfunction = "coxph")
set.seed(1234)
glmulti.conv_preop.aicc = glmulti(
	data.y ~ 
		Patient.Sex + 
		History.Diagnosis.AgeAt.Cent + 
#		ns(History.Diagnosis.AgeAt.Cent, 4) + 
		Path.LocationBody + 
		Path.Size.Cent + 
#		ns(Path.Size.Cent, 4) + 
		Stage.pT.Simplified,
	data = data.x.conv_preop,
	level = 2,
	method = "g",
	crit = "aicc", 
	includeobjects = TRUE, 
	plotty = FALSE, 
	report = TRUE,
	marginality = TRUE,
	fitfunction = "coxph")


set.seed(1234)
glmulti.molec_postop.bic = glmulti(
	data.y ~ 
		Patient.Sex + 
		History.Diagnosis.AgeAt.Cent + 
#		ns(History.Diagnosis.AgeAt.Cent, 4) + 
		Path.LocationBody + 
		Path.Size.Cent + 
#		ns(Path.Size.Cent, 4) + 
		Stage.pT.Simplified + 
		Molec.S100A4.DCThresh + 
		Molec.S100A2.DCThresh + 
		Treat.MarginPositive + 
#		Path.Subtype + 					# This field needs data review if it's to be kept.
		Path.Differentiation + 
		Path.LN.InvolvedFraction + 
#		ns(Path.LN.InvolvedFraction, 4) + 
		Path.Invasion.Vascular + 
		Path.Invasion.Perineural + 
		Stage.pN,
	data = data.x.molec_postop,
	level = 2,
	method = "g",
	crit = "bic", 
	includeobjects = TRUE, 
	plotty = FALSE, 
	report = TRUE,
	marginality = TRUE,
	fitfunction = "coxph")
set.seed(1234)
glmulti.molec_postop.aicc = glmulti(
	data.y ~ 
		Patient.Sex + 
		History.Diagnosis.AgeAt.Cent + 
#		ns(History.Diagnosis.AgeAt.Cent, 4) + 
		Path.LocationBody + 
		Path.Size.Cent + 
#		ns(Path.Size.Cent, 4) + 
		Stage.pT.Simplified + 
		Molec.S100A4.DCThresh + 
		Molec.S100A2.DCThresh + 
		Treat.MarginPositive + 
#		Path.Subtype + 					# This field needs data review if it's to be kept.
		Path.Differentiation + 
		Path.LN.InvolvedFraction + 
#		ns(Path.LN.InvolvedFraction, 4) + 
		Path.Invasion.Vascular + 
		Path.Invasion.Perineural + 
		Stage.pN,
	data = data.x.molec_postop,
	level = 2,
	method = "g",
	crit = "aicc", 
	includeobjects = TRUE, 
	plotty = FALSE, 
	report = TRUE,
	marginality = TRUE,
	fitfunction = "coxph")


set.seed(1234)
glmulti.conv_postop.bic = glmulti(
	data.y ~ 
		Patient.Sex + 
		History.Diagnosis.AgeAt.Cent + 
#		ns(History.Diagnosis.AgeAt.Cent, 4) + 
		Path.LocationBody + 
		Path.Size.Cent + 
#		ns(Path.Size.Cent, 4) + 
		Stage.pT.Simplified + 
		Treat.MarginPositive + 
#		Path.Subtype +  				# This field needs data review if it's to be kept.
		Path.Differentiation + 
		Path.LN.InvolvedFraction + 
#		ns(Path.LN.InvolvedFraction, 4) + 
		Path.Invasion.Vascular + 
		Path.Invasion.Perineural + 
		Stage.pN,
	data = data.x.conv_postop,
	level = 2,
	method = "g",
	crit = "bic", 
	includeobjects = TRUE, 
	plotty = FALSE, 
	report = TRUE,
	marginality = TRUE,
	fitfunction = "coxph")
set.seed(1234)
glmulti.conv_postop.aicc = glmulti(
	data.y ~ 
		Patient.Sex + 
		History.Diagnosis.AgeAt.Cent + 
#		ns(History.Diagnosis.AgeAt.Cent, 4) + 
		Path.LocationBody + 
		Path.Size.Cent + 
#		ns(Path.Size.Cent, 4) + 
		Stage.pT.Simplified + 
		Treat.MarginPositive + 
#		Path.Subtype +  				# This field needs data review if it's to be kept.
		Path.Differentiation + 
		Path.LN.InvolvedFraction + 
#		ns(Path.LN.InvolvedFraction, 4) + 
		Path.Invasion.Vascular + 
		Path.Invasion.Perineural + 
		Stage.pN,
	data = data.x.conv_postop,
	level = 2,
	method = "g",
	crit = "aicc", 
	includeobjects = TRUE, 
	plotty = FALSE, 
	report = TRUE,
	marginality = TRUE,
	fitfunction = "coxph")


glmulti.conv_preop.bic@formulas[1:10]
plot(glmulti.conv_preop.bic, "p")
plot(glmulti.conv_preop.bic, "r")
plot(glmulti.conv_preop.bic, "s")
plot(glmulti.conv_preop.bic, "w")

glmulti.conv_postop.bic@formulas[1:10]
plot(glmulti.conv_postop.bic, "p")
plot(glmulti.conv_postop.bic, "r")
plot(glmulti.conv_postop.bic, "s")
plot(glmulti.conv_postop.bic, "w")

glmulti.molec_preop.bic@formulas[1:10]
plot(glmulti.molec_preop.bic, "p")
plot(glmulti.molec_preop.bic, "r")
plot(glmulti.molec_preop.bic, "s")
plot(glmulti.molec_preop.bic, "w")

glmulti.molec_postop.bic@formulas[1:10]
plot(glmulti.molec_postop.bic, "p")
plot(glmulti.molec_postop.bic, "r")
plot(glmulti.molec_postop.bic, "s")
plot(glmulti.molec_postop.bic, "w")



glmulti.conv_preop.aic@formulas[1:10]
plot(glmulti.conv_preop.aic, "p")
plot(glmulti.conv_preop.aic, "r")
plot(glmulti.conv_preop.aic, "s")
plot(glmulti.conv_preop.aic, "w")

glmulti.conv_postop.aic@formulas[1:10]
plot(glmulti.conv_postop.aic, "p")
plot(glmulti.conv_postop.aic, "r")
plot(glmulti.conv_postop.aic, "s")
plot(glmulti.conv_postop.aic, "w")

glmulti.molec_preop.aic@formulas[1:10]
plot(glmulti.molec_preop.aic, "p")
plot(glmulti.molec_preop.aic, "r")
plot(glmulti.molec_preop.aic, "s")
plot(glmulti.molec_preop.aic, "w")

glmulti.molec_postop.aic@formulas[1:10]
plot(glmulti.molec_postop.aic, "p")
plot(glmulti.molec_postop.aic, "r")
plot(glmulti.molec_postop.aic, "s")
plot(glmulti.molec_postop.aic, "w")


rm(nobs.coxph, temp1, temp2, temp)

save.image("04_NSWPCN_fits.rda")
