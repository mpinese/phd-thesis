# Construct a *preoperative* function based on the Brennan nomogram.  The
# preoperative nature will mean that most prognostic components will need
# to be marginalized out.

# Variable				Preoperative?	Available?		Marginals
# Age					Yes				Yes				Linear.  90 => 0, 30 => 8.  Therefore f(x) = (x-90)*-2/15 = -2/15x + 12
# Sex					Yes				Yes				Male risk delta 3
# Portal Vein			NO								14.4% YES, risk delta 10, marginal 1.4
# Splenectomy			NO								9.9% YES, risk delta 62, marginal 6.1
# Margin of resection	NO								20.7% POS, risk delta 4, marginal 0.8
# Head.vs.Other			Yes				Yes				Head risk delta 51
# Differentiation		NO 								14.2% Well, risk delta 0, marginal 0
# 														56.4% Mod, risk delta 14, marginal 7.9
# 														29.5% Poor, risk delta 35, marginal 10.3.  Overall marginal 18.2
# Posterior.margin		NO								86.0% POS, risk delta 22, marginal 18.9
# Numb.pos.nodes		NO								Mean 2.1, approx marginal 15
# Numb.neg.nodes		NO								Mean 16.9, approx marginal 9
# Back.pain				Yes				NO				13.7% YES, risk delta 15, marginal 2.0
# T.stage				Yes				Yes
# Weight Loss			Yes				NO 				53.7% YES, risk delta 3,  marginal 1.6
# Max.path.axis			Yes				Yes

# So the preoperative MSKCC score would be:
#   S = 1.4 + 6.1 + 0.8 + 18.2 + 18.9 + 15 + 9 + 15*Back.pain + 3*Weight.Loss + -2/15*Age + 12 + 3*I(Sex == "M") + 51*I(Head.vs.Other == "Head") + T.stage + Max.path.axis
#     = 81.4 + 15*Back.pain + 3*Weight.Loss + -2/15*Age + 3*I(Sex == "M") + 51*I(Head.vs.Other == "Head") + fT(T.stage) + fS(Max.path.axis)
#
# fT(T.stage) == 36*I(T.stage == "T1") + 10*I(T.stage == "T3") + 63*I(T.stage == "T4")
# fS(Max.path.axis) is a spline.  I've digitized the points on the Brennan paper:
#


nomogram.mskcc = list(
	inputs = list(
	History.Diagnosis.AgeAt = list(
		margins = data.frame(value = 65, fraction = 1),
		scorefunc = function(x) { x = x; -2/15*pmin(pmax(x, 0), 90) + 12 }),
	Patient.Sex = list(
		margins = data.frame(value = c("M", "F"), fraction = c(0.501, 1-0.501)),
		scorefunc = function(x) { 3*I(x == "M") }),
	Portal.Vein = list(
		margins = data.frame(value = c(TRUE, FALSE), fraction = c(0.144, 1-0.144)),
		scorefunc = function(x) { 10*I(x == TRUE) }),
	Splenectomy = list(
		margins = data.frame(value = c(TRUE, FALSE), fraction = c(0.099, 1-0.099)),
		scorefunc = function(x) { 62*I(x == TRUE) }),
	Treat.MarginPositive = list(
		margins = data.frame(value = c(TRUE, FALSE), fraction = c(0.207, 1-0.207)),
		scorefunc = function(x) { 4*I(x == TRUE) }),
	Path.LocationBody = list(
		margins = data.frame(value = c(FALSE, TRUE), fraction = c(0.894, 1-0.894)),
		scorefunc = function(x) { 51*I(x == TRUE) }),
	Path.Differentiation = list(
		margins = data.frame(value = c("1", "2", "3", "4"), fraction = c(0.142, 0.564, 1-0.142-0.564, 0)),
		scorefunc = function(x) { 14*I(x == "2") + 35*I(x == "3") + 35*I(x == "4") }),		# Undifferentiated (4) not covered by the MSKCC nomogram; here assign the same score as poorly differentiated (3)
	Posterior.Margin = list(
		margins = data.frame(value = c(TRUE, FALSE), fraction = c(0.86, 1-0.86)),
		scorefunc = function(x) { 22*I(x == TRUE) }),
	Path.LN.Involved = list(
		margins = data.frame(value = 2.1, fraction = 1),
		scorefunc = function(x) { 
			x = pmin(40, pmax(x, 0))
			fitfun = splinefun(c(0, 1, 2, 3, 4, 10, 15, 20, 25, 30, 35, 40), c(0, 14.56, 24.64, 30.28, 33.00, 39.05, 43.89, 48.83, 53.77, 58.61, 63.55, 68.49), method = "natural")
			fitfun(x)
		}),
	Path.LN.Negative = list(
		margins = data.frame(value = 16.9, fraction = 1),
		scorefunc = function(x) { (pmin(pmax(x, 0), 90)-90)*-11/90 }),
	Back.pain = list(
		margins = data.frame(value = c(TRUE, FALSE), fraction = c(0.137, 1-0.137)),
		scorefunc = function(x) { 15*I(x == TRUE) }),
	Stage.pT.Simplified = list(
		margins = data.frame(value = c("T1", "T2", "T34"), fraction = c(0.037, 0.119, 1-0.037-0.119)),
		scorefunc = function(x) { 36*I(x == "T1") + 11*I(x == "T34") }),
		# The following matches the original Brennan nomogram, but was not used as there are too few T4
		# tumours in either the NSWPCN *or* the MSKCC cohorts -- how the T4 coefficient was ever estimated,
		# I'll never know.  The T34 coefficient of 11 was arrived at as (0.828*10+(1-0.037-0.119-0.828)*63)/(1-0.037-0.119),
		# being a frequency-weighted average of the T3 and T4 coefficients.
		# margins = data.frame(value = c("T1", "T2", "T3", "T4"), fraction = c(0.037, 0.119, 0.828, 1-0.037-0.119-0.828)),
		# scorefunc = function(x) { 36*I(x == "T1") + 10*I(x == "T3") + 63*I(x == "T4") }),
	Weight.loss = list(
		margins = data.frame(value = c(TRUE, FALSE), fraction = c(0.537, 1-0.537)),
		scorefunc = function(x) { 3*I(x == TRUE) }),
	Path.Size = list(
		margins = data.frame(),
		scorefunc = function(x) {
			x = pmin(16, pmax(x, 0))
			fitfun = splinefun(c(0, 1, 2, 3, 4, 6, 8, 10, 12, 14, 16), c(0, 29.74, 59.48, 86.70, 100, 97.29, 90.03, 82.77, 75.51, 68.25, 61.10), method = "natural")
			fitfun(x)
		}) ),
	outputs = list(
		DSS12mo = function(s) {
			x = pmax(50, pmin(350, s))
			fitfun = splinefun(c(79.0323, 115.02, 165.524, 197.278, 221.774, 242.339, 261.089, 279.839, 299.194, 323.992, 337.298), c(0.94, 0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1, 0.06))
			y = fitfun(x)
			pmax(0, pmin(1, y))
		},
		DSS24mo = function(s) {
			x = pmax(50, pmin(350, s))
			fitfun = splinefun(c(71.1694, 97.7823, 129.536, 153.73, 174.294, 193.347, 211.794, 231.452, 255.645, 303.125), c(0.86, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1, 0.01))
			y = fitfun(x)
			pmax(0, pmin(1, y))
		},
		DSS36mo = function(s) {
			x = pmax(50, pmin(350, s))
			fitfun = splinefun(c(69.3548, 101.109, 125.302, 145.867, 164.919, 183.367, 202.722, 226.915, 274.093), c(0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1, 0.01))
			y = fitfun(x)
			pmax(0, pmin(1, y))
		}) 
	)


applyNomogram = function(nomogram, data)
{
	scores = rowSums(sapply(names(nomogram$inputs), function(input) {
		if (input %in% colnames(data)) {
			return(nomogram$inputs[[input]]$scorefunc(data[,input]))
		}
		warning(sprintf("Marginalizing missing variable: %s", input))
		margin_score = sum(nomogram$inputs[[input]]$scorefunc(nomogram$inputs[[input]]$margins$value) * nomogram$inputs[[input]]$margins$fraction)
		return(rep(margin_score, nrow(data)))
	}))

	outputs = sapply(nomogram$outputs, function(f) f(scores))
	cbind(Score = scores, outputs)
}


load("03_NSWPCN_subset.rda")
mskcc.scores.nswpcn.postop = applyNomogram(nomogram.mskcc, data.x.all)
mskcc.scores.nswpcn.preop = applyNomogram(nomogram.mskcc, cbind(data.x.conv_preop, data.x.all$Path.Size, data.x.all$History.Diagnosis.AgeAt))

data.glasgow = readRDS("05_Glasgow.rds")
data.glasgow$Path.LN.Negative = data.glasgow$Path.LN.Inspected - data.glasgow$Path.LN.Involved
data.glasgow$History.Diagnosis.AgeAt = data.glasgow$History.Diagnosis.AgeAt.Cent + 68
data.glasgow$Path.Size = data.glasgow$Path.Size.Cent + 30
mskcc.scores.glasgow.postop = applyNomogram(nomogram.mskcc, data.glasgow)
mskcc.scores.glasgow.preop = applyNomogram(nomogram.mskcc, data.glasgow[,c("History.Diagnosis.AgeAt", "Patient.Sex", "Path.LocationBody", "Stage.pT.Simplified", "Path.Size")])


library(survival)
coxph(data.y ~ mskcc.scores.nswpcn.preop[,1])
coxph(data.y ~ mskcc.scores.nswpcn.postop[,1])
coxph(Surv(data.glasgow$History.Death.EventTimeDays, data.glasgow$History.DSDeath.Event) ~ mskcc.scores.glasgow.preop[,1])
coxph(Surv(data.glasgow$History.Death.EventTimeDays, data.glasgow$History.DSDeath.Event) ~ mskcc.scores.glasgow.postop[,1])


pdf("07_test_MSKCC_on_NSWPCN_Glasgow.pdf")

library(glmulti)
library(randomForestSRC)
library(survival)

library(timeROC)
mskcc.cdroc.nswpcn.preop = timeROC(data.y[,1]/365.25*12, data.y[,2], mskcc.scores.nswpcn.preop[,1], cause = 1, times = seq(1, 36, 1), iid = TRUE)
mskcc.cdroc.nswpcn.postop = timeROC(data.y[,1]/365.25*12, data.y[,2], mskcc.scores.nswpcn.postop[,1], cause = 1, times = seq(1, 36, 1), iid = TRUE)
mskcc.cdroc.glasgow.preop = timeROC(data.glasgow$History.Death.EventTimeDays/365.25*12, as.numeric(as.character(data.glasgow$History.Death.Cause)), mskcc.scores.glasgow.preop[,1], cause = 1, times = seq(1, 36, 1), iid = TRUE)
mskcc.cdroc.glasgow.postop = timeROC(data.glasgow$History.Death.EventTimeDays/365.25*12, as.numeric(as.character(data.glasgow$History.Death.Cause)), mskcc.scores.glasgow.postop[,1], cause = 1, times = seq(1, 36, 1), iid = TRUE)
plotAUCcurve(mskcc.cdroc.nswpcn.preop, conf.int = FALSE, add = FALSE, col = "red")
plotAUCcurve(mskcc.cdroc.nswpcn.postop, conf.int = FALSE, add = TRUE, col = "blue")
plotAUCcurve(mskcc.cdroc.glasgow.preop, conf.int = FALSE, add = TRUE, col = "orange")
plotAUCcurve(mskcc.cdroc.glasgow.postop, conf.int = FALSE, add = TRUE, col = "purple")
legend("topright", legend = c("NSWPCN Preop", "NSWPCN Postop", "Glasgow Preop", "Glasgow Postop"), col = c("red", "blue", "orange", "purple"), lty = "solid")

library(risksetROC)
risksetROC(data.y[,1]/365.25*12, status = data.y[,2], marker = mskcc.scores.nswpcn.preop[,1], predict.time = 12)
risksetROC(data.y[,1]/365.25*12, status = data.y[,2], marker = mskcc.scores.nswpcn.postop[,1], predict.time = 12)
risksetAUC(data.y[,1]/365.25*12, status = data.y[,2], marker = mskcc.scores.nswpcn.preop[,1], tmax = 36)
risksetAUC(data.y[,1]/365.25*12, status = data.y[,2], marker = mskcc.scores.nswpcn.postop[,1], tmax = 36)

risksetROC(data.glasgow$History.Death.EventTimeDays/365.25*12, status = data.glasgow$History.DSDeath.Event, marker = mskcc.scores.glasgow.preop[,1], predict.time = 12)
risksetROC(data.glasgow$History.Death.EventTimeDays/365.25*12, status = data.glasgow$History.DSDeath.Event, marker = mskcc.scores.glasgow.postop[,1], predict.time = 12)
risksetAUC(data.glasgow$History.Death.EventTimeDays/365.25*12, status = data.glasgow$History.DSDeath.Event, marker = mskcc.scores.glasgow.preop[,1], tmax = 36)
risksetAUC(data.glasgow$History.Death.EventTimeDays/365.25*12, status = data.glasgow$History.DSDeath.Event, marker = mskcc.scores.glasgow.postop[,1], tmax = 36)
