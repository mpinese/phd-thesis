load("03_NSWPCN_subset.rda")

#pdf("05_train_NSWPCN_2.pdf")

options(echo = TRUE, warn = 1)

library(flexsurv)

data.all = cbind(time = data.y[,1], event = data.y[,2], data.x.all, data.x.extra)
data.all$Path.Size.Cent.Thresh = pmin(data.all$Path.Size.Cent, 0)

data.all$event[data.all$time > 2000] = 0
data.all$time = pmin(2000, data.all$time)


# Optimised by AIC, starting from model:
# Patient.Sex + History.Diagnosis.AgeAt.Cent + Path.LocationBody + Path.Size.Cent.Thresh + Molec.S100A4.DCThresh + Molec.S100A2.DCThresh
fit1 = flexsurvspline(Surv(time, event) ~ Patient.Sex + Path.Size.Cent.Thresh + Molec.S100A4.DCThresh + Molec.S100A2.DCThresh
	+ gamma2(Patient.Sex), 
	data = data.all, k = 1, scale = "odds")
plot(fit1)


# Optimised by AIC, starting from model:
# Patient.Sex + History.Diagnosis.AgeAt.Cent + Path.LocationBody + Path.Size.Cent.Thresh + Molec.S100A4.DCThresh + Molec.S100A2.DCThresh
# lnorm chosen based on Q != ~0 result from gengamma fit.
fit2 = flexsurvreg(Surv(time, event) ~ Path.Size.Cent.Thresh + Molec.S100A4.DCThresh + Molec.S100A2.DCThresh
	+ sdlog(Patient.Sex),
	data = data.all, dist = "lnorm")
plot(fit2)


summary(fit1, newdata = data.all[1,,drop=FALSE])
summary(fit2, newdata = data.all[2,,drop=FALSE])

#save.image("05_NSWPCN_fits_2.rda")
