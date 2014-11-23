######################################################################
# LIBRARIES
######################################################################
options(java.parameters = "-Xmx4G")
library(xtable)
library(glmulti)
library(glmnet)
library(energy)
library(RColorBrewer)
library(NMF)
library(gplots)
library(stargazer)

nmf.options(cores = 8, pbackend = "par", gc = 1, shared.memory = TRUE)


######################################################################
# DATA
######################################################################
load("../../data/07_data_for_SIS.rda")
source("../../src/08_SIS_common_funcs.R")


######################################################################
# HIGH LEVEL PARAMETERS
######################################################################
data_source = commandArgs(TRUE)[1]
x = get(sprintf("x.%s", data_source))
y = get(sprintf("y.%s", data_source))
samps = get(sprintf("samps.%s", data_source))
tau = 0.72
theta = 0.05
x0 = 6.335
all_sigs = x.msigdb
seed = 1234567890
nmf.nrun.rank = 50
nmf.nrun.fit = 500
nmf.rank = "auto"
nmf.rankrange = 2:10
nmf.rankrandcount = 10
sig.corr.threshold = 0.5


######################################################################
# DERIVED VARIABLES
######################################################################
xlin = 2^(x-x0)
xlin.scaled = (xlin - apply(xlin, 1, min)) / as.vector(diff(apply(xlin, 1, range)))
sigs = all_sigs[,colnames(x)]


######################################################################
# PROBE SELECTION
######################################################################
set.seed(seed)
cpss.sis = CPSS(x, y, SIS.FAST, tau, 50, nsis = floor(theta*nrow(x)))

x.sel = x[cpss.sis$sel,]
xlin.sel = xlin[cpss.sis$sel,]
xlin.scaled.sel = xlin.scaled[cpss.sis$sel,]


######################################################################
# EXPRESSION CORRELATION
######################################################################
x.sel.kcor = cor(t(x.sel), method = "kendall")
x.sel.dcor = sapply(1:(nrow(x.sel)-1), function(i) c(rep(NA, i), sapply((i+1):nrow(x.sel), function(j) dcor(x.sel[i,], x.sel[j,]))))
x.sel.dcor = cbind(x.sel.dcor, NA)
diag(x.sel.dcor) = 1
x.sel.dcor[upper.tri(x.sel.dcor)] = t(x.sel.dcor)[upper.tri(x.sel.dcor)]


######################################################################
# RANK ESTIMATION
######################################################################
temp.nmf.rank = nmf(
	x = xlin.scaled.sel, 
	rank = nmf.rankrange, 
	method = "snmf/l", 
	seed = seed, nrun = nmf.nrun.rank, 
	.options = list(verbose = 1, track = TRUE, parallel = TRUE, keep.all = TRUE))
temp.nmf.rank.random = lapply(1:nmf.rankrandcount, function(i) {
	message(i)
	nmf(x = randomize(xlin.scaled.sel), 
		rank = nmf.rankrange, 
		method = "snmf/l", 
		seed = seed, nrun = nmf.nrun.rank, 
		.options = list(verbose = 1, track = TRUE, parallel = TRUE, keep.all = TRUE))
	})

temp.resids = sapply(temp.nmf.rank$fit, function(f) sapply(f, residuals))
temp.resids_rel = t(t(temp.resids) / apply(temp.resids, 2, min))
temp.resids_scaled = t((t(temp.resids) - apply(temp.resids, 2, min)) / (apply(temp.resids, 2, max) - apply(temp.resids, 2, min)))
temp.orig_resids = sapply(temp.nmf.rank$fit, residuals)
temp.perm_resids = sapply(temp.nmf.rank.random, function(rep) sapply(rep$fit, residuals))
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
xlin.scaled.sel.nmf = nmf(
	x = xlin.scaled.sel, 
	rank = nmf.rank, 
	method = "snmf/l", 
	seed = seed, nrun = nmf.nrun.fit, 
	.options = list(verbose = 0, track = TRUE, parallel = TRUE, keep.all = TRUE))


######################################################################
# SIGNATURE-METAGENE CORRELATION
######################################################################
xlin.scaled.sel.nmf.msigdb.corr = cor(t(coef(xlin.scaled.sel.nmf)), t(sigs), method = "kendall")


######################################################################
# ALL-SUBSETS REGRESSION
######################################################################
asreg.data = data.frame(time = y[,1], event = y[,2], mg = t(coef(xlin.scaled.sel.nmf)))
nobs.coxph = function(obj)	{ obj$nevent }
asreg.result = glmulti(Surv(time, event) ~ ., data = asreg.data, fitfunction = "coxph", level = 2, marginality = TRUE, crit = bic, plotty = FALSE, report = FALSE)
rm(nobs.coxph)


######################################################################
# LASSO
######################################################################
glmnet.x = t(coef(xlin.scaled.sel.nmf))
colnames(glmnet.x) = paste("mg", 1:ncol(glmnet.x), sep = ".")
glmnet.fit.cv = cv.glmnet(x = glmnet.x, y = cbind(time = y[,1], status = y[,2]*1), family = "cox", nfolds = 10)
glmnet.coef.1se = coef(glmnet.fit.cv$glmnet.fit, s = glmnet.fit.cv$lambda.1se)
glmnet.coef.min = coef(glmnet.fit.cv$glmnet.fit, s = glmnet.fit.cv$lambda.min)


######################################################################
# ADAPTIVE LASSO
######################################################################
adaglmnet.weights = 1/abs(coef(coxph(y ~ glmnet.x)))
adaglmnet.x = t(t(glmnet.x) * adaglmnet.weights)
adaglmnet.fit.cv = cv.glmnet(x = adaglmnet.x, y = cbind(time = y[,1], status = y[,2]*1), family = "cox", nfolds = 10, standardize = FALSE)
adaglmnet.coef.1se = coef(adaglmnet.fit.cv$glmnet.fit, s = adaglmnet.fit.cv$lambda.1se)
adaglmnet.coef.min = coef(adaglmnet.fit.cv$glmnet.fit, s = adaglmnet.fit.cv$lambda.min)


session_info = sessionInfo()
save.image("image.rda")
