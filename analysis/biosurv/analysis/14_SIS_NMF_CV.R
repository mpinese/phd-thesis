######################################################################
# LIBRARIES
######################################################################
options(java.parameters = "-Xmx4G", echo = TRUE)
library(survival)
library(NMF)
library(glmnet)
library(glmulti)
library(nnls)

nmf.options(cores = 32, pbackend = "par", gc = 1, shared.memory = FALSE)


######################################################################
# DATA
######################################################################
load("../data/07_data_for_SIS.rda")
source("../common/08_SIS_common_funcs.R")


######################################################################
# HIGH LEVEL PARAMETERS
######################################################################
tau = 0.72
theta = 0.05
seed = 1234567890
nmf.nrun.rank = 50
nmf.nrun.fit = 500
nmf.rankrange = 2:10
nmf.rankrandcount = 5


sis_nmf_fitpred = function(x, y, xtest, theta, tau, nmf.nrun.rank, nmf.nrun.fit, nmf.rankrange, nmf.rankrandcount, seed)
{
	######################################################################
	# DERIVED VARIABLES
	######################################################################
	xlin = 2^x
	xtestlin = 2^xtest
	xtestlin.scaled = (xtestlin - apply(xlin, 1, min)) / as.vector(diff(apply(xlin, 1, range)))
	xlin.scaled = (xlin - apply(xlin, 1, min)) / as.vector(diff(apply(xlin, 1, range)))

	######################################################################
	# PROBE SELECTION
	######################################################################
	set.seed(seed)
	cpss.sis = CPSS(x, y, SIS.FAST, tau, 50, nsis = floor(theta*nrow(x)))

	x.sel = x[cpss.sis$sel,]
	xlin.sel = xlin[cpss.sis$sel,]
	xlin.scaled.sel = xlin.scaled[cpss.sis$sel,]
	xtestlin.scaled.sel = xtestlin.scaled[cpss.sis$sel,]

	######################################################################
	# RANK ESTIMATION
	######################################################################
	message("Initial factorizations...")
	temp.nmf.rank = nmf(
		x = xlin.scaled.sel, 
		rank = nmf.rankrange, 
		method = "snmf/l", 
		seed = seed, nrun = nmf.nrun.rank, 
		.options = list(verbose = 1, track = FALSE, parallel = TRUE, keep.all = TRUE))
	message("Random factorizations...")
	temp.nmf.rank.random = lapply(1:nmf.rankrandcount, function(i) {
		message(i)
		nmf(x = randomize(xlin.scaled.sel), 
			rank = nmf.rankrange, 
			method = "snmf/l", 
			seed = seed, nrun = nmf.nrun.rank, 
			.options = list(verbose = 1, track = FALSE, parallel = TRUE, keep.all = TRUE))
		})
	temp.orig_resids = sapply(temp.nmf.rank$fit, residuals)
	temp.perm_resids = sapply(temp.nmf.rank.random, function(rep) sapply(rep$fit, residuals))
	temp.perm_resids_mean = rowMeans(temp.perm_resids)
	temp.orig_resids.spline = splinefun(nmf.rankrange, temp.orig_resids, method = "natural")
	temp.perm_resids.spline = apply(temp.perm_resids, 2, function(r) splinefun(nmf.rankrange, r, method = "natural"))
	temp.perm_resids_mean.spline = splinefun(nmf.rankrange, temp.perm_resids_mean)
	temp.orig_resids.d1 = temp.orig_resids.spline(nmf.rankrange, deriv = 1)
	temp.perm_resids_mean.d1 = temp.perm_resids_mean.spline(nmf.rankrange, deriv = 1)
	temp.perm_resids.d1 = sapply(temp.perm_resids.spline, function(x) x(nmf.rankrange, deriv = 1))
	temp.perm_resids.d1.mean = rowMeans(temp.perm_resids.d1)
	temp.perm_resids.d1.sd = apply(temp.perm_resids.d1, 1, sd)
	nmf.rank = nmf.rankrange[min(which(temp.orig_resids.d1 >= temp.perm_resids.d1.mean - 2*temp.perm_resids.d1.sd))]
	if (length(nmf.rank) == 0)	{ nmf.rank = max(nmf.rankrange) }

	######################################################################
	# NMF FACTORIZATION
	######################################################################
	xlin.scaled.sel.nmf = nmf(
		x = xlin.scaled.sel, 
		rank = nmf.rank, 
		method = "snmf/l", 
		seed = seed, nrun = nmf.nrun.fit, 
		.options = list(verbose = 2, track = FALSE, parallel = TRUE, keep.all = TRUE))

	######################################################################
	# TRAIN X SCORING
	######################################################################
	xlin.scores = t(apply(xlin.scaled.sel, 2, function(xcol) nnls(basis(xlin.scaled.sel.nmf), xcol)$x))
#	xlin.scores = t(coef(xlin.scaled.sel.nmf))
	colnames(xlin.scores) = paste("mg", 1:ncol(xlin.scores), sep = ".")

	######################################################################
	# TEST X SCORING
	######################################################################
	xtestlin.scores = t(apply(xtestlin.scaled.sel, 2, function(xcol) nnls(basis(xlin.scaled.sel.nmf), xcol)$x))
#	xtestlin.scores = t(xtestlin.scaled.sel) %*% basis(xlin.scaled.sel.nmf)
	colnames(xtestlin.scores) = paste("mg", 1:ncol(xtestlin.scores), sep = ".")

	######################################################################
	# BEST SUBSET SELECTION
	######################################################################
	asreg.data <<- cbind(data.frame(time = y[,1], event = y[,2]), xlin.scores)
	nobs.coxph <<- function(obj)	{ obj$nevent }
	asreg.result = glmulti(Surv(time, event) ~ ., data = asreg.data, fitfunction = "coxph", level = 2, marginality = TRUE, crit = bic, plotty = FALSE, report = FALSE)
#	rm(nobs.coxph)
	pred.bs.best = predict(asreg.result@objects[[1]], newdata = as.data.frame(xtestlin.scores))
	pred.bs.average = as.vector(predict(asreg.result, newdata = as.data.frame(xtestlin.scores))$averages)

	######################################################################
	# LASSO
	######################################################################
	glmnet.fit.cv = cv.glmnet(x = xlin.scores, y = cbind(time = y[,1], status = y[,2]*1), family = "cox", nfolds = 10)
	pred.lasso.1se = as.vector(predict(glmnet.fit.cv$glmnet.fit, newx = xtestlin.scores, s = glmnet.fit.cv$lambda.1se, type = "link"))
	pred.lasso.min = as.vector(predict(glmnet.fit.cv$glmnet.fit, newx = xtestlin.scores, s = glmnet.fit.cv$lambda.min, type = "link"))

	######################################################################
	# ADAPTIVE LASSO
	######################################################################
	adaglmnet.weights = 1/abs(coef(coxph(y ~ xlin.scores)))
	adaglmnet.x = t(t(xlin.scores) * adaglmnet.weights)
	adaglmnet.newx = t(t(xtestlin.scores) * adaglmnet.weights)
	adaglmnet.fit.cv = cv.glmnet(x = adaglmnet.x, y = cbind(time = y[,1], status = y[,2]*1), family = "cox", nfolds = 10, standardize = FALSE)
	pred.adalasso.1se = as.vector(predict(adaglmnet.fit.cv$glmnet.fit, newx = adaglmnet.newx, s = adaglmnet.fit.cv$lambda.1se))
	pred.adalasso.min = as.vector(predict(adaglmnet.fit.cv$glmnet.fit, newx = adaglmnet.newx, s = adaglmnet.fit.cv$lambda.min))

	as.matrix(cbind(
		bs.best = pred.bs.best, 
		bs.average = pred.bs.average, 
		lasso.1se = pred.lasso.1se, 
		lasso.min = pred.lasso.min, 
		adalasso.1se = pred.adalasso.1se, 
		adalasso.min = pred.adalasso.min), ncol = 6)
}


doCV = function(x, y, fitpred, K, seed)
{
	set.seed(seed)
	n = nrow(y)
	folds = sample(rep(1:K, ceiling(n/K))[1:n])

	yhat_runs = lapply(1:K, function(k) {
		message(sprintf("Fold %d of %d", k, K))
		ind.test = folds == k
		ind.train = !ind.test

		x.test = x[,ind.test,drop = FALSE]
		x.train = x[,ind.train,drop = FALSE]
		y.test = y[ind.test,,drop = FALSE]
		y.train = y[ind.train,,drop = FALSE]

		yhat.test = fitpred(x.train, y.train, x.test)
		yhat.test
	})

	yhats = matrix(NA, nrow = ncol(yhat_runs[[1]]), ncol = n)
	rownames(yhats) = colnames(yhat_runs[[1]])
	for (k in 1:K)
	{
		yhats[,folds == k] = yhat_runs[[k]]
		colnames(yhats)[folds == k] = colnames(x)[folds == k]
	}

	yhats
}


cv.diag_dsd = doCV(x.diag_dsd, y.diag_dsd, 
	function(x, y, x.test) sis_nmf_fitpred(x, y, x.test, theta, tau, nmf.nrun.rank, nmf.nrun.fit, nmf.rankrange, nmf.rankrandcount, seed), 
	10, seed)


sessioninfo = sessionInfo()
sessioninfo

saveRDS(cv.diag_dsd, file = "14_SIS_NMF_CV_results.rds")

