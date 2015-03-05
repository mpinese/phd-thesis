load("../biosurv/data/07_data_for_SIS.rda")

x = x.diag_dsd
y = y.diag_dsd
samps = samps.diag_dsd

temp = NA
temp = ls()
rm(list = temp[!(temp %in% c("x", "y", "samps"))])

library(survival)
temp.pct = seq(0.001, 0.999, 0.001)
temp.reg = survreg(y ~ 1)
temp.reg
temp.ptime = predict(temp.reg, type = "quantile", p = temp.pct, newdata = data.frame(dummy = 1))
plot(survfit(y ~ 1), main = "Marginal survival function")
lines(temp.ptime, 1 - temp.pct)
temp.reg = survreg(Surv(y[,1], !y[,2]) ~ 1)
temp.reg
temp.ptime = predict(temp.reg, type = "quantile", p = temp.pct, newdata = data.frame(dummy = 1))
plot(survfit(Surv(y[,1], !y[,2]) ~ 1), main = "Marginal censoring function")
lines(temp.ptime, 1 - temp.pct)


library(messina)
library(doMC)
registerDoMC(32)

messina.coxcoef.1 = messinaSurv(x, y, obj_min = 1, obj_func = "coxcoef", seed = 20150301)
messina.tau.7 = messinaSurv(x, y, obj_min = 0.7, obj_func = "tau", seed = 20150301)
messina.reltau.8 = messinaSurv(x, y, obj_min = 0.8, obj_func = "reltau", seed = 20150301)

save.image("01_surv_apgi_calcs.rda")

