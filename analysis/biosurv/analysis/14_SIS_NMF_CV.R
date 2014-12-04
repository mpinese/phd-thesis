######################################################################
# LIBRARIES
######################################################################
options(java.parameters = "-Xmx4G", echo = TRUE, warn = 1)
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
tau = 0.9
gamma = 7
seed = 1234567890
nmf.nrun.rank = 50
nmf.nrun.fit = 500
nmf.rankrange = 2:9
nmf.rankrandcount = 5
nmf.beta = 0.01


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
	cpss.sis = CPSS(x, y, SIS.FAST, tau, 50, gamma = gamma, scale = TRUE)

	x.sel = x[cpss.sis$sel,]
	xlin.sel = xlin[cpss.sis$sel,]
	xlin.scaled.sel = xlin.scaled[cpss.sis$sel,]
	xtestlin.scaled.sel = xtestlin.scaled[cpss.sis$sel,]

	######################################################################
	# RANK ESTIMATION
	######################################################################
	message("Initial factorizations...")
	nmf.runs.rank = nmf(
	        x = xlin.sel,
        	rank = nmf.rankrange,
	        method = "snmf/l",
        	seed = seed, nrun = nmf.nrun.rank,
	        .options = list(verbose = 1, track = FALSE, parallel = TRUE, keep.all = FALSE),
        	beta = nmf.beta)
	message("Random factorizations...")
	nmf.runs.rank.random = lapply(1:nmf.rankrandcount, function(i) {
        	message(i)
	        nmf(x = randomize(xlin.sel),  
                	rank = nmf.rankrange, 
        	        method = "snmf/l",  
	                seed = seed, nrun = nmf.nrun.rank,
                	.options = list(verbose = 1, track = FALSE, parallel = TRUE, keep.all = FALSE),
        	        beta = nmf.beta)
	        })

	temp.orig_resids = sapply(nmf.runs.rank$fit, residuals)
	temp.perm_resids = sapply(nmf.runs.rank.random, function(rep) sapply(rep$fit, residuals))
	temp.orig_resids.delta = diff(temp.orig_resids)
	temp.perm_resids.delta = apply(temp.perm_resids, 2, diff)
	temp.perm_resids.delta.mean = rowMeans(temp.perm_resids.delta)
	temp.perm_resids.delta.sd = apply(temp.perm_resids.delta, 1, sd)
	temp.perm_resids.delta.threshold = temp.perm_resids.delta.mean - 2*temp.perm_resids.delta.sd
	temp.perm_resids.delta.above_threshold = temp.orig_resids.delta >= temp.perm_resids.delta.threshold
	if (all(temp.perm_resids.delta.above_threshold))                { nmf.rank.auto = min(nmf.rankrange)
	} else if (all(!(temp.perm_resids.delta.above_threshold)))	{ nmf.rank.auto = max(nmf.rankrange)
	} else                                                          { nmf.rank.auto = min(nmf.rankrange[temp.perm_resids.delta.above_threshold]) }
        nmf.rank = nmf.rank.auto

	######################################################################
	# NMF FACTORIZATION
	######################################################################
	xlin.scaled.sel.nmf = nmf(
		x = xlin.scaled.sel, 
		rank = nmf.rank, 
		method = "snmf/l", 
		seed = seed, nrun = nmf.nrun.fit, 
		.options = list(verbose = 2, track = FALSE, parallel = TRUE, keep.all = TRUE),
		beta = nmf.beta)

	######################################################################
	# TRAIN X SCORING
	######################################################################
	xlin.scores = t(apply(xlin.scaled.sel, 2, function(xcol) nnls(basis(xlin.scaled.sel.nmf), xcol)$x))
	colnames(xlin.scores) = paste("mg", 1:ncol(xlin.scores), sep = ".")

	######################################################################
	# TEST X SCORING
	######################################################################
	xtestlin.scores = t(apply(xtestlin.scaled.sel, 2, function(xcol) nnls(basis(xlin.scaled.sel.nmf), xcol)$x))
	colnames(xtestlin.scores) = paste("mg", 1:ncol(xtestlin.scores), sep = ".")

	######################################################################
	# LASSO
	######################################################################
	glmnet.fit.cv = cv.glmnet(x = xlin.scores, y = cbind(time = y[,1], status = y[,2]*1), family = "cox", nfolds = 10)
	pred.lasso.1se = as.vector(predict(glmnet.fit.cv$glmnet.fit, newx = xtestlin.scores, s = glmnet.fit.cv$lambda.1se, type = "link"))
	pred.lasso.min = as.vector(predict(glmnet.fit.cv$glmnet.fit, newx = xtestlin.scores, s = glmnet.fit.cv$lambda.min, type = "link"))

	as.matrix(cbind(
		lasso.1se = pred.lasso.1se, 
		lasso.min = pred.lasso.min), ncol = 2)
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

