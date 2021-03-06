args = commandArgs(TRUE)
mode = 0
if (length(args) > 0)
{
	mode = as.numeric(args[1]) 
}

load("03_NSWPCN_subset.rda")

options(echo = TRUE, warn = 1)

library(survival)

if (mode == 0 || mode == 9)
{
pdf("04_train_NSWPCN.pdf")

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
rm(temp, temp1, temp2)
}

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

if (mode == 0 || mode == 1) {
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
	fitfunction = "coxph",
	name = "glmulti.molec_preop.bic")
saveRDS(glmulti.molec_preop.bic, file = "04_molec_preop_bic.rds")
}

if (mode == 0 || mode == 2) {
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
	fitfunction = "coxph",
	name = "glmulti.molec_preop.aicc")
saveRDS(glmulti.molec_preop.aicc, file = "04_molec_preop_aicc.rds")
}

if (mode == 0 || mode == 3) {
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
	fitfunction = "coxph",
	name = "glmulti.conv_preop.bic")
saveRDS(glmulti.conv_preop.bic, file = "04_conv_preop_bic.rds")
}

if (mode == 0 || mode == 4) {
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
	fitfunction = "coxph",
	name = "glmulti.conv_preop.aicc")
saveRDS(glmulti.conv_preop.aicc, file = "04_conv_preop_aicc.rds")
}

if (mode == 0 || mode == 5) { 
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
	fitfunction = "coxph",
	name = "glmulti.molec_postop.bic")
saveRDS(glmulti.molec_postop.bic, file = "04_molec_postop_bic.rds")
}

if (mode == 0 || mode == 6) {
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
	fitfunction = "coxph",
	name = "glmulti.molec_postop.aicc")
saveRDS(glmulti.molec_postop.aicc, file = "04_molec_postop_aicc.rds")
}

if (mode == 0 || mode == 7) {
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
	fitfunction = "coxph",
	name = "glmulti.conv_postop.bic")
saveRDS(glmulti.conv_postop.bic, file = "04_conv_postop_bic.rds")
}

if (mode == 0 || mode == 8) {
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
	fitfunction = "coxph",
	name = "glmulti.conv_postop.aicc")
saveRDS(glmulti.conv_postop.aicc, file = "04_conv_postop_aicc.rds")
}

if (mode == 0 || mode == 10)
{
	library(randomForestSRC)

	temp.data = cbind(time = data.y[,1], event = data.y[,2], data.x.molec_preop)
	# temp.K = 10

	# temp.tasks = expand.grid(mtry = c(1, 2, 3), splitrule = c("logrank", "logrankscore"), nsplit = c(0, 1, 2, 5))
	# set.seed(1234)
	# temp.folds = sample(rep(1:temp.K, ceiling(nrow(temp.data) / temp.K))[1:nrow(temp.data)])

	# preds = lapply(1:temp.K, function(fold_i) {
	# 	message(sprintf("Fold %d", fold_i))
	# 	train = temp.data[fold_i != temp.folds,]
	# 	test = temp.data[fold_i == temp.folds,]

	# 	fold_preds = sapply(1:nrow(temp.tasks), function(task_i) {
	# 		message(sprintf("  Task %d", task_i))
	# 		fit = rfsrc(
	# 			Surv(time, event) ~ 
	# 				Patient.Sex + 
	# 				History.Diagnosis.AgeAt.Cent + 
	# 				Path.LocationBody + 
	# 				Path.Size.Cent + 
	# 				Stage.pT.Simplified + 
	# 				Molec.S100A4.DCThresh + 
	# 				Molec.S100A2.DCThresh, 
	# 			data = train,
	# 			mtry = temp.tasks$mtry[task_i],
	# 			splitrule = temp.tasks$splitrule[task_i],
	# 			nsplit = temp.tasks$nsplit[task_i])
	# 		preds = predict(fit, test)$predicted
	# 		preds
	# 	})

	# 	fold_preds
	# })

	# preds2 = do.call(rbind, preds)
	# temp.y2 = do.call(rbind, lapply(1:temp.K, function(fold_i) data.y[fold_i == temp.folds,]))
	# temp.y2 = Surv(temp.y2[,1], temp.y2[,2])
	# temp.results = cbind(temp.tasks, C = apply(preds2, 2, function(p1) survConcordance(temp.y2 ~ p1)$concordance))
	# temp.results[order(temp.results$C),]

	# From the above: mtry = 1, splitrule = "logrankscore", nsplit = 2 gives good performance.

	rsf.molec_preop = rfsrc(
		Surv(time, event) ~ 
			Patient.Sex + 
			History.Diagnosis.AgeAt.Cent + 
			Path.LocationBody + 
			Path.Size.Cent + 
			Stage.pT.Simplified + 
			Molec.S100A4.DCThresh + 
			Molec.S100A2.DCThresh, 
		data = temp.data,
		mtry = 1,
		splitrule = "logrankscore",
		nsplit = 2)
	saveRDS(rsf.molec_preop, file = "04_rsf_molec_preop.rds")
}

if (mode == 0 || mode == 9) {
glmulti.molec_preop.bic = readRDS("04_molec_preop_bic.rds")
glmulti.molec_preop.aicc = readRDS("04_molec_preop_aicc.rds")
glmulti.molec_postop.bic = readRDS("04_molec_postop_bic.rds")
glmulti.molec_postop.aicc = readRDS("04_molec_postop_aicc.rds")
glmulti.conv_preop.bic = readRDS("04_conv_preop_bic.rds")
glmulti.conv_preop.aicc = readRDS("04_conv_preop_aicc.rds")
glmulti.conv_postop.bic = readRDS("04_conv_postop_bic.rds")
glmulti.conv_postop.aicc = readRDS("04_conv_postop_aicc.rds")
rsf.molec_preop = readRDS("04_rsf_molec_preop.rds")

if ("nobs.coxph" %in% ls()) { rm(nobs.coxph) }

glmulti.conv_preop.bic@formulas[1:10]
plot(glmulti.conv_preop.bic, "p")
#plot(glmulti.conv_preop.bic, "r")
plot(glmulti.conv_preop.bic, "s")
plot(glmulti.conv_preop.bic, "w")

glmulti.conv_postop.bic@formulas[1:10]
plot(glmulti.conv_postop.bic, "p")
#plot(glmulti.conv_postop.bic, "r")
plot(glmulti.conv_postop.bic, "s")
plot(glmulti.conv_postop.bic, "w")

glmulti.molec_preop.bic@formulas[1:10]
plot(glmulti.molec_preop.bic, "p")
#plot(glmulti.molec_preop.bic, "r")
plot(glmulti.molec_preop.bic, "s")
plot(glmulti.molec_preop.bic, "w")

glmulti.molec_postop.bic@formulas[1:10]
plot(glmulti.molec_postop.bic, "p")
#plot(glmulti.molec_postop.bic, "r")
plot(glmulti.molec_postop.bic, "s")
plot(glmulti.molec_postop.bic, "w")

glmulti.conv_preop.aicc@formulas[1:10]
plot(glmulti.conv_preop.aicc, "p")
#plot(glmulti.conv_preop.aicc, "r")
plot(glmulti.conv_preop.aicc, "s")
plot(glmulti.conv_preop.aicc, "w")

glmulti.conv_postop.aicc@formulas[1:10]
plot(glmulti.conv_postop.aicc, "p")
#plot(glmulti.conv_postop.aicc, "r")
plot(glmulti.conv_postop.aicc, "s")
plot(glmulti.conv_postop.aicc, "w")

glmulti.molec_preop.aicc@formulas[1:10]
plot(glmulti.molec_preop.aicc, "p")
#plot(glmulti.molec_preop.aicc, "r")
plot(glmulti.molec_preop.aicc, "s")
plot(glmulti.molec_preop.aicc, "w")

glmulti.molec_postop.aicc@formulas[1:10]
plot(glmulti.molec_postop.aicc, "p")
#plot(glmulti.molec_postop.aicc, "r")
plot(glmulti.molec_postop.aicc, "s")
plot(glmulti.molec_postop.aicc, "w")

save.image("04_NSWPCN_fits.rda")
}
