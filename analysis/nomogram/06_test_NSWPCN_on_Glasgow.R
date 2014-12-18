load("04_NSWPCN_fits.rda")
data.glasgow = readRDS("05_Glasgow.rds")

pdf("06_test_NSWPCN_on_Glasgow.pdf")

rm(args, data, data.x.confounders, data.x.conv_preop, data.x.conv_postop, 
	data.x.management, data.x.molec_preop, data.x.molec_postop, data.x.postop.conventional,
	data.x.preop.conventional, data.x.preop.molecular, mode, origcols, tsub)

library(glmulti)
library(randomForestSRC)
library(survival)

model0 = coxph(data.y ~ 1)
model1 = glmulti.conv_preop.bic
model2 = glmulti.molec_preop.bic
model3 = glmulti.conv_postop.bic
model4 = glmulti.molec_postop.bic
model5 = rsf.molec_preop

preds0 = rep(0, nrow(data.glasgow))
preds1 = predict(model1, select = 1, newdata = data.glasgow)$averages[1,]
preds2 = predict(model2, select = 1, newdata = data.glasgow)$averages[1,]
preds3 = predict(model3, select = 1, newdata = data.glasgow)$averages[1,]
preds4 = predict(model4, select = 1, newdata = data.glasgow)$averages[1,]
preds5 = predict(model5, newdata = data.glasgow)$predicted


library(timeROC)
cdroc0 = timeROC(data.glasgow$History.Death.EventTimeDays/365.25*12, as.numeric(as.character(data.glasgow$History.Death.Cause)), preds0, cause = 1, times = seq(1, 36, 1), iid = TRUE)
cdroc1 = timeROC(data.glasgow$History.Death.EventTimeDays/365.25*12, as.numeric(as.character(data.glasgow$History.Death.Cause)), preds1, cause = 1, times = seq(1, 36, 1), iid = TRUE)
cdroc2 = timeROC(data.glasgow$History.Death.EventTimeDays/365.25*12, as.numeric(as.character(data.glasgow$History.Death.Cause)), preds2, cause = 1, times = seq(1, 36, 1), iid = TRUE)
cdroc3 = timeROC(data.glasgow$History.Death.EventTimeDays/365.25*12, as.numeric(as.character(data.glasgow$History.Death.Cause)), preds3, cause = 1, times = seq(1, 36, 1), iid = TRUE)
cdroc4 = timeROC(data.glasgow$History.Death.EventTimeDays/365.25*12, as.numeric(as.character(data.glasgow$History.Death.Cause)), preds4, cause = 1, times = seq(1, 36, 1), iid = TRUE)
cdroc5 = timeROC(data.glasgow$History.Death.EventTimeDays/365.25*12, as.numeric(as.character(data.glasgow$History.Death.Cause)), preds5, cause = 1, times = seq(1, 36, 1), iid = TRUE)
#plotAUCcurve(cdroc0, conf.int = TRUE)
plotAUCcurve(cdroc1, conf.int = FALSE, add = FALSE, col = "red")		# Preop-conv
plotAUCcurve(cdroc2, conf.int = FALSE, add = TRUE, col = "blue")		# Preop-molec
#plotAUCcurve(cdroc3, conf.int = FALSE, add = TRUE, col = "orange")		# Postop-conv
#plotAUCcurve(cdroc4, conf.int = FALSE, add = TRUE, col = "green")		# Postop-molec
plotAUCcurve(cdroc5, conf.int = FALSE, add = TRUE, col = "purple")		# Preop-molec (RSF)
legend("topright", legend = c("Preop-Conv", "Preop-Molec", "Preop-Molec (RSF)"), col = c("red", "blue", "purple"), lty = "solid")
plotAUCcurveDiff(cdroc2, cdroc1, conf.int = TRUE)
plotAUCcurveDiff(cdroc5, cdroc2, conf.int = TRUE)


library(risksetROC)
risksetROC(data.glasgow$History.Death.EventTimeDays/365.25*12, status = data.glasgow$History.Death.Cause == "1", marker = preds1, predict.time = 12)
risksetROC(data.glasgow$History.Death.EventTimeDays/365.25*12, status = data.glasgow$History.Death.Cause == "1", marker = preds2, predict.time = 12)
risksetROC(data.glasgow$History.Death.EventTimeDays/365.25*12, status = data.glasgow$History.Death.Cause == "1", marker = preds5, predict.time = 12)
risksetAUC(data.glasgow$History.Death.EventTimeDays/365.25*12, status = data.glasgow$History.Death.Cause == "1", marker = preds1, tmax = 36)
risksetAUC(data.glasgow$History.Death.EventTimeDays/365.25*12, status = data.glasgow$History.Death.Cause == "1", marker = preds2, tmax = 36)
risksetAUC(data.glasgow$History.Death.EventTimeDays/365.25*12, status = data.glasgow$History.Death.Cause == "1", marker = preds5, tmax = 36)


library(pec)


