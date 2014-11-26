######################################################################
# LIBRARIES
######################################################################
options(java.parameters = "-Xmx4G")
library(glmulti)
library(glmnet)
library(energy)
library(RColorBrewer)
library(gplots)
library(ahaz)
library(survival)
library(bigmemory)


######################################################################
# DATA
######################################################################
message("Loading data")
load("../../data/07_data_for_SIS.rda")
source("../../src/08_SIS_common_funcs.R")

######################################################################
# HIGH LEVEL PARAMETERS
######################################################################
rownames(y.diag_dsd) = colnames(x.diag_dsd)
rownames(y.diag_rec) = colnames(x.diag_rec)
rownames(y.recr_dsd) = colnames(x.recr_dsd)
samps = samples[colnames(x),]
y_list = list(diag_dsd = y.diag_dsd, diag_rec = y.diag_rec, recr_dsd = y.recr_dsd)
tau = 0.72
theta1 = 0.025
x0 = 6.335
all_sigs = x.msigdb
seed = 1234567890
nmf.nrun.rank = 50
nmf.nrun.fit = 500
nmf.rank = "auto"
nmf.rankrange = 2:15
nmf.rankrandcount = 10
sig.corr.threshold = 0.5


######################################################################
# DERIVED VARIABLES
######################################################################
message("Calculating derived variables")
xlin = 2^(x-x0)
xlin.scaled = (xlin - apply(xlin, 1, min)) / as.vector(diff(apply(xlin, 1, range)))
sigs = all_sigs[,colnames(x)]


######################################################################
# PROBE SELECTION
######################################################################
subset2.SIS.FAST = function(x, y, subset, ...)
{
	y_sub = y[intersect(rownames(y), subset),,drop = FALSE]
	x_sub = x[,rownames(y_sub)]
	valid = apply(!is.na(x_sub), 2, all) & (!is.na(y_sub[,1])) & (!is.na(y_sub[,2]))
	x_sub = x_sub[,valid,drop = FALSE]
	y_sub = y_sub[valid,,drop = FALSE]
	SIS.FAST(x_sub, y_sub, ...)
}

subset.SIS.FAST = function(x, ys, subset, ...)
{
	sels = sapply(ys, function(y) subset2.SIS.FAST(x, y, subset, ...))
	apply(sels, 1, any)
}

group.CPSS = function(x, ys, selfunc, tau, B = 50, ...)
{
	n = ncol(x)
	samp_ids = colnames(x)

	samps = lapply(1:B, function(i) sort(sample(samp_ids, floor(n/2))))

	Shats = sapply(1:(2*B), function(i) {
		A = samps[[ceiling(i/2)]]
		if (i %% 2 == 0)	{ A = setdiff(samp_ids, A) }

		selfunc(x, ys, subset = A, ...)
	})

	Pihat = rowSums(Shats) / (2*B)
	qs = colSums(Shats)

	result = list(sel = Pihat > tau, Pihat = Pihat, samples = samps, qs = qs, qhat = mean(qs))
	result
}

message("CPSS")
set.seed(seed)
cpss.sis = group.CPSS(x, y_list, subset.SIS.FAST, tau, 50, nsis = floor(theta1*nrow(x)))
x.sel = x[cpss.sis$sel,]
xlin.sel = xlin[cpss.sis$sel,]
xlin.scaled.sel = xlin.scaled[cpss.sis$sel,]
message(sprintf("%d/%d selected, observed theta = %.3f", sum(cpss.sis$sel), length(cpss.sis$sel), cpss.sis$qhat / length(cpss.sis$sel)))


######################################################################
# RANK ESTIMATION
######################################################################
source("snmfl.R")
message("NMF rank estimation")
set.seed(seed)
temp.nmf.rank = snmfl(
	A = xlin.scaled.sel, 
	ks = nmf.rankrange, 
	nrun = nmf.nrun.rank, 
	maxiter = 1e4, cores = 8)
message("  Randomized")
set.seed(seed)
temp.nmf.rank.random = lapply(1:nmf.rankrandcount, function(i) {
	message(i)
	Aperm = apply(xlin.scaled.sel, 2, sample)
	snmfl(
		A = Aperm, 
		ks = nmf.rankrange, 
		nrun = nmf.nrun.rank, 
		maxiter = 1e4, cores = 8)
	}
)

save.image("temp.rda")

message("Automatic rank calculation")
temp.resids = sapply(temp.nmf.rank, function(kf) kf$norms)
temp.resids_rel = t(t(temp.resids) / apply(temp.resids, 2, min))
temp.resids_scaled = t((t(temp.resids) - apply(temp.resids, 2, min)) / (apply(temp.resids, 2, max) - apply(temp.resids, 2, min)))
temp.orig_resids = sapply(temp.nmf.rank, function(kf) min(kf$norms))
temp.perm_resids = sapply(temp.nmf.rank.random, function(rep) sapply(rep, function(kf) min(kf$norms)))
temp.orig_resids.delta = diff(temp.orig_resids)
temp.perm_resids.delta = apply(temp.perm_resids, 2, diff)
temp.perm_resids.delta.mean = rowMeans(temp.perm_resids.delta)
temp.perm_resids.delta.sd = apply(temp.perm_resids.delta, 1, sd)
temp.perm_resids.delta.threshold = temp.perm_resids.delta.mean - 2*temp.perm_resids.delta.sd
temp.perm_resids.delta.above_threshold = temp.orig_resids.delta >= temp.perm_resids.delta.threshold
if (all(temp.perm_resids.delta.above_threshold))			{ nmf.rank.auto = min(nmf.rankrange) 
} else if (all(!(temp.perm_resids.delta.above_threshold)))	{ nmf.rank.auto = max(nmf.rankrange)
} else 														{ nmf.rank.auto = min(nmf.rankrange[temp.perm_resids.delta.above_threshold]) }
if (nmf.rank == "auto")
{
	nmf.rank = nmf.rank.auto
}

######################################################################
# FACTORIZATION
######################################################################
message("Full factorization")
xlin.scaled.sel.nmf = snmfl(
	A = xlin.scaled.sel, 
	k = nmf.rank, 
	nrun = nmf.nrun.fit, 
	maxiter = 1e4, cores = 8)

save.image("temp2.rda")

coefs = t(xlin.scaled.sel.nmf$best_fit$H)
colnames(coefs) = paste("mg", 1:ncol(coefs), sep = ".")
coefs.diag_dsd = coefs[rownames(y.diag_dsd),]
coefs.diag_rec = coefs[rownames(y.diag_rec),]
coefs.recr_dsd = coefs[rownames(y.recr_dsd),]


######################################################################
# EXPRESSION CORRELATION
######################################################################
message("Correlation")
x.sel.kcor = cor(t(x.sel), method = "kendall")
x.sel.dcor = sapply(1:(nrow(x.sel)-1), function(i) c(rep(NA, i), sapply((i+1):nrow(x.sel), function(j) dcor(x.sel[i,], x.sel[j,]))))
x.sel.dcor = cbind(x.sel.dcor, NA)
diag(x.sel.dcor) = 1
x.sel.dcor[upper.tri(x.sel.dcor)] = t(x.sel.dcor)[upper.tri(x.sel.dcor)]


######################################################################
# SIGNATURE-METAGENE CORRELATION
######################################################################
xlin.scaled.sel.nmf.msigdb.corr = cor(coefs, t(sigs), method = "kendall")


######################################################################
# ALL-SUBSETS REGRESSION
######################################################################
diag_dsd.asreg.data = as.data.frame(cbind(time = y.diag_dsd[,1], event = y.diag_dsd[,2], coefs.diag_dsd))
diag_rec.asreg.data = as.data.frame(cbind(time = y.diag_rec[,1], event = y.diag_rec[,2], coefs.diag_rec))
recr_dsd.asreg.data = as.data.frame(cbind(time = y.recr_dsd[,1], event = y.recr_dsd[,2], coefs.recr_dsd))
nobs.coxph = function(obj)	{ obj$nevent }
diag_dsd.asreg.result = glmulti(Surv(time, event) ~ ., data = diag_dsd.asreg.data, fitfunction = "coxph", level = 2, marginality = TRUE, crit = bic, plotty = FALSE, report = FALSE)
diag_rec.asreg.result = glmulti(Surv(time, event) ~ ., data = diag_rec.asreg.data, fitfunction = "coxph", level = 2, marginality = TRUE, crit = bic, plotty = FALSE, report = FALSE)
recr_dsd.asreg.result = glmulti(Surv(time, event) ~ ., data = recr_dsd.asreg.data, fitfunction = "coxph", level = 2, marginality = TRUE, crit = bic, plotty = FALSE, report = FALSE)
rm(nobs.coxph)


######################################################################
# LASSO
######################################################################
diag_dsd.glmnet.fit.cv = cv.glmnet(x = coefs.diag_dsd, y = cbind(time = y.diag_dsd[,1], status = y.diag_dsd[,2]*1), family = "cox", nfolds = 10)
diag_rec.glmnet.fit.cv = cv.glmnet(x = coefs.diag_rec, y = cbind(time = y.diag_rec[,1], status = y.diag_rec[,2]*1), family = "cox", nfolds = 10)
recr_dsd.glmnet.fit.cv = cv.glmnet(x = coefs.recr_dsd, y = cbind(time = y.recr_dsd[,1], status = y.recr_dsd[,2]*1), family = "cox", nfolds = 10)
diag_dsd.glmnet.coef.1se = coef(diag_dsd.glmnet.fit.cv$glmnet.fit, s = diag_dsd.glmnet.fit.cv$lambda.1se)
diag_dsd.glmnet.coef.min = coef(diag_dsd.glmnet.fit.cv$glmnet.fit, s = diag_dsd.glmnet.fit.cv$lambda.min)
diag_rec.glmnet.coef.1se = coef(diag_rec.glmnet.fit.cv$glmnet.fit, s = diag_rec.glmnet.fit.cv$lambda.1se)
diag_rec.glmnet.coef.min = coef(diag_rec.glmnet.fit.cv$glmnet.fit, s = diag_rec.glmnet.fit.cv$lambda.min)
recr_dsd.glmnet.coef.1se = coef(recr_dsd.glmnet.fit.cv$glmnet.fit, s = recr_dsd.glmnet.fit.cv$lambda.1se)
recr_dsd.glmnet.coef.min = coef(recr_dsd.glmnet.fit.cv$glmnet.fit, s = recr_dsd.glmnet.fit.cv$lambda.min)


######################################################################
# ADAPTIVE LASSO
######################################################################
diag_dsd.adaglmnet.weights = 1/abs(coef(coxph(y.diag_dsd ~ coefs.diag_dsd)))
diag_rec.adaglmnet.weights = 1/abs(coef(coxph(y.diag_rec ~ coefs.diag_rec)))
recr_dsd.adaglmnet.weights = 1/abs(coef(coxph(y.recr_dsd ~ coefs.recr_dsd)))
diag_dsd.adaglmnet.x = t(t(coefs.diag_dsd) * diag_dsd.adaglmnet.weights)
diag_rec.adaglmnet.x = t(t(coefs.diag_rec) * diag_rec.adaglmnet.weights)
recr_dsd.adaglmnet.x = t(t(coefs.recr_dsd) * recr_dsd.adaglmnet.weights)
diag_dsd.adaglmnet.fit.cv = cv.glmnet(x = diag_dsd.adaglmnet.x, y = cbind(time = y.diag_dsd[,1], status = y.diag_dsd[,2]*1), family = "cox", nfolds = 10, standardize = FALSE)
diag_rec.adaglmnet.fit.cv = cv.glmnet(x = diag_rec.adaglmnet.x, y = cbind(time = y.diag_rec[,1], status = y.diag_rec[,2]*1), family = "cox", nfolds = 10, standardize = FALSE)
recr_dsd.adaglmnet.fit.cv = cv.glmnet(x = recr_dsd.adaglmnet.x, y = cbind(time = y.recr_dsd[,1], status = y.recr_dsd[,2]*1), family = "cox", nfolds = 10, standardize = FALSE)
diag_dsd.adaglmnet.coef.1se = coef(diag_dsd.adaglmnet.fit.cv$glmnet.fit, s = diag_dsd.adaglmnet.fit.cv$lambda.1se)
diag_dsd.adaglmnet.coef.min = coef(diag_dsd.adaglmnet.fit.cv$glmnet.fit, s = diag_dsd.adaglmnet.fit.cv$lambda.min)
diag_rec.adaglmnet.coef.1se = coef(diag_rec.adaglmnet.fit.cv$glmnet.fit, s = diag_rec.adaglmnet.fit.cv$lambda.1se)
diag_rec.adaglmnet.coef.min = coef(diag_rec.adaglmnet.fit.cv$glmnet.fit, s = diag_rec.adaglmnet.fit.cv$lambda.min)
recr_dsd.adaglmnet.coef.1se = coef(recr_dsd.adaglmnet.fit.cv$glmnet.fit, s = recr_dsd.adaglmnet.fit.cv$lambda.1se)
recr_dsd.adaglmnet.coef.min = coef(recr_dsd.adaglmnet.fit.cv$glmnet.fit, s = recr_dsd.adaglmnet.fit.cv$lambda.min)

session_info = sessionInfo()
save.image("image.rda")
