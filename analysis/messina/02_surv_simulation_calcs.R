library(messina)




generator = function(delta_x, noise_sd, surv_dists, censor_dist, balance, n, p, ...)
{
	n1 = round(balance * n)
	n2 = n - n1
	true_class = rep(c(0, 1), c(n1, n2))
	p = p1 + p2

	x = matrix(NA, nrow = p, ncol = n)
	for (i in 1:p1)	{ x[i,] = true_class*delta_x + rnorm(n, sd = noise_sd) }
	x[(p1+1):nrow(x),] = rnorm(n*p2)
	x = (x - rowMeans(x)) / apply(x, 1, sd) * median(apply(x[1:p1,], 1, sd)) + rowMeans(x)
	true_x = rep(c(TRUE, FALSE), c(p1, p2))
	
	time_event = rep(NA, n)
	time_event[true_class == 0] = surv_dists[["0"]](n1)
	time_event[true_class == 1] = surv_dists[["1"]](n2)
	time_cens = censor_dist(n)

	time_observed = pmin(time_event, time_cens)
	event_observed = time_event <= time_cens

	y = Surv(time_observed, event_observed)

	fit = messinaSurv(x, y, ...)

	detection_tp = sum(fit@fits@summary$passed & true_x)
	detection_fp = sum(fit@fits@summary$passed & !true_x)
	detection_fn = sum(!fit@fits@summary$passed & true_x)
	detection_tn = sum(!fit@fits@summary$passed & !true_x)
	perf = list(tp = detection_tp, fp = detection_fp, fn = detection_fn, tn = detection_tn)

	list(fit = fit, x = x, y = y, true_class = true_class == 1, true_x = true_x, perf = perf)
}




# Detection validation
detector_mediancut = function(x, y, p = 0.05)
{
	apply(x, 1, function(x1) {
		x1 = x1 > median(x1)
		test = survdiff(y ~ x1)
		pval = pchisq(test$chisq, df = 1, lower.tail = FALSE)
		pval < p
	})
}

detector_singlecox = function(x, y, p = 0.05)
{
	apply(x, 1, function(x1) pchisq(diff(coxph(y ~ x1)$loglik)*2, df = 1, lower.tail = FALSE) < p)
}

detector_multicut = function(x, y, p = 0.05, ncuts = 10, correct = c("none", "BY"))
{
	correct = match.arg(correct)
	apply(x, 1, function(x1) {
		cutpoints = quantile(x1, probs = (1:ncuts)/(ncuts + 1))
		pvals = sapply(cutpoints, function(c) {
			x1c = x1 > c
			test = survdiff(y ~ x1c)
			pval = pchisq(test$chisq, df = 1, lower.tail = FALSE)
			pval
		})
		if (correct == "BY")	{ pvals = p.adjust(pvals, "BY") }
		min(pvals) < p
	})
}

detector_messinaSurv = function(x, y, minmargin = 2, mincox = 1)
{
	fits = messinaSurv(x, y, objective = messinaSurvObj.CoxCoef(mincox))
	fits@fits@summary$passed & fits@fits@summary$margin >= minmargin
}


detection_experiments = expand.grid(
	delta.x = 5,
	noise.sd = 1,
	class.fraction = c(0.2, 0.5),
	cohort.size = c(25, 100, 400),
	hazard.ratio = c(1, 1.5, 3, 6),
	censoring.rate = c(0.2, 0.5, 0.8),
	experiment.run = 1:200)
nrow(detection_experiments)



apply(detection_experiments, 1, function(exp) {

}





library(messina)
library(antiProfilesData)
data(apColonData)
x = exprs(apColonData)
y = pData(apColonData)$SubType
sel = y %in% c("normal", "tumor")
x = x[,sel]
y = y[sel]
fit.apColon = messina(x, y == "tumor", min_sens = 0.95, min_spec = 0.85, seed = 1234, silent = TRUE)
plot(fit.apColon, plot_type = "point")
# fit.apColon = messina(x, y == "tumor", min_sens = 0.95, min_spec = 0.6, seed = 1234, silent = TRUE)
# plot(fit.apColon, plot_type = "point")
# fit.apColon = messina(x, y == "tumor", min_sens = 0.8, min_spec = 0.6, seed = 1234, silent = TRUE)
# plot(fit.apColon, plot_type = "point")

# x = x[rep("204719_at", 3),]
# #rownames(x) <- NULL
# #debug(messina)
# fit.apColon = messina(x, y == "tumor", min_sens = 0.95, min_spec = 0.85, seed = 1234, silent = TRUE)
# fit.apColon@fits@summary

# plot(fit.apColon, indices = 3503, sort_features = FALSE, plot_type = "point")

temp = fit.apColon@fits@summary
temp[order(temp$pvalpass, temp$margin, decreasing = TRUE),][1:20,]

# x_truth = rep(c(0, 1), c(40, 40))
# x_obs = matrix(rnorm(length(x_truth) * 100), nrow = 10, ncol = length(x_truth))
# x_obs[1,] = x_truth*10 + rnorm(length(x_truth))
# x_obs = (x_obs - rowMeans(x_obs)) / apply(x_obs, 1, sd) + rowMeans(x_obs)
# y = (x_truth + rnorm(length(x_truth), sd = 0.5) > 0.5)*1
# messina(x_obs, y==1, min_sens = 0.8, min_spec = 0.8)

# Two concepts:
#   * Detection rate -- at what rate does the method do what it says on the tin?
#   * Robustness -- are survival classifiers more robust than those from other methods?

# For detection rate:
#   * Stick to CoxCoef, because this is the one for which it's easiest to quantify 
#     the survival difference required, from a distribution sense.
#   * Detectability will be affected by cohort size, class balance, censoring rate, 
#     censoring pattern, group SNR, UGH...
#   * Maybe show the performance for a reasonable combination of values, then
#     demonstrate that it's relatively robust to edge cases on the others.

# Other methods:
#   * Cox coefficient filtering -- Tibs (?) suggested
#   * Cox P-value filtering -- in common use
#   * Logrank all-cuts testing -- in common use

simulator = function(delta_x, noise_sd, surv_dists, censor_dist, balance, n, p1, p2, ...)
{
	n1 = round(balance * n)
	n2 = n - n1
	true_class = rep(c(0, 1), c(n1, n2))
	p = p1 + p2

	x = matrix(NA, nrow = p, ncol = n)
	for (i in 1:p1)	{ x[i,] = true_class*delta_x + rnorm(n, sd = noise_sd) }
	x[(p1+1):nrow(x),] = rnorm(n*p2)
	x = (x - rowMeans(x)) / apply(x, 1, sd) * median(apply(x[1:p1,], 1, sd)) + rowMeans(x)
	true_x = rep(c(TRUE, FALSE), c(p1, p2))
	
	time_event = rep(NA, n)
	time_event[true_class == 0] = surv_dists[["0"]](n1)
	time_event[true_class == 1] = surv_dists[["1"]](n2)
	time_cens = censor_dist(n)

	time_observed = pmin(time_event, time_cens)
	event_observed = time_event <= time_cens

	y = Surv(time_observed, event_observed)

	fit = messinaSurv(x, y, ...)

	detection_tp = sum(fit@fits@summary$passed & true_x)
	detection_fp = sum(fit@fits@summary$passed & !true_x)
	detection_fn = sum(!fit@fits@summary$passed & true_x)
	detection_tn = sum(!fit@fits@summary$passed & !true_x)
	perf = list(tp = detection_tp, fp = detection_fp, fn = detection_fn, tn = detection_tn)

	list(fit = fit, x = x, y = y, true_class = true_class == 1, true_x = true_x, perf = perf)
}



set.seed(20150301)
temp = simulator(
	delta_x = 15, 
	noise_sd = 1, 
	surv_dists = list("0" = function(n) rweibull(n, shape = 1/0.7253, scale = exp(6.9)/sqrt(exp(1))), "1" = function(n) rweibull(n, shape = 1/0.7253, scale = exp(6.9)*sqrt(exp(1)))), 
	censor_dist = function(n) rweibull(n, shape = 1/0.3216, scale = exp(7.1)), 
	balance = 0.5, 
	n = 200, 
	p1 = 50, p2 = 50, 
	objective = messinaSurvObj.CoxCoef(1))


experiments = expand.grid(
	delta_x = c(1, 2, 4, 8),
	noise_sd = 1,
	surv_dist = c("A", "B", "C"),
	cens_dist = "A",
	balance = c(0.2, 0.5, 0.8),
	n = c(20, 50, 100),
	p1 = 50, p2 = 50)


surv_dist_library = list(
	"A" = list("0" = function(n) rweibull(n, shape = 1/0.7253, scale = exp(6.9)/sqrt(exp(0.5))), "1" = function(n) rweibull(n, shape = 1/0.7253, scale = exp(6.9)*sqrt(exp(0.5)))),
	"B" = list("0" = function(n) rweibull(n, shape = 1/0.7253, scale = exp(6.9)/sqrt(exp(1))), "1" = function(n) rweibull(n, shape = 1/0.7253, scale = exp(6.9)*sqrt(exp(1)))),
	"C" = list("0" = function(n) rweibull(n, shape = 1/0.7253, scale = exp(6.9)/sqrt(exp(1.5))), "1" = function(n) rweibull(n, shape = 1/0.7253, scale = exp(6.9)*sqrt(exp(1.5)))))
cens_dist_library = list("A" = function(n) rweibull(n, shape = 1/0.3216, scale = exp(7.1)))


set.seed(20150301)
temp = simulator(
	delta_x = 15, 
	noise_sd = 1, 
	surv_dists = list("0" = function(n) rweibull(n, shape = 1/0.7253, scale = exp(6.9)/sqrt(exp(1))), "1" = function(n) rweibull(n, shape = 1/0.7253, scale = exp(6.9)*sqrt(exp(1)))), 
	censor_dist = function(n) rweibull(n, shape = 1/0.3216, scale = exp(7.1)), 
	balance = 0.5, 
	n = 200, 
	p1 = 50, p2 = 50, 
	objective = messinaSurvObj.CoxCoef(1))

plot(temp$fit)
#plot(survfit(temp$y ~ temp$truex))
