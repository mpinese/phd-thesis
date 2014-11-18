#!~/bin/Rscript

# 09_SIS_diag_dsd.R

options(echo = TRUE)
X11.options(type = "dbcairo")

load("../data/07_data_for_SIS.rda")
source("08_SIS_common_funcs.R")



######################################################################
# HIGH LEVEL PARAMETERS
######################################################################

x = x.diag_dsd
y = y.diag_dsd
tau = 0.72
theta = 0.05
x0 = 6.335
ncomp = 5
all_sigs = x.msigdb.c123467.merged
seed = 1234



######################################################################
# DERIVED VARIABLES
######################################################################

xlin = 2^(x-x0)
xlin = (xlin - apply(xlin, 1, min)) / as.vector(diff(apply(xlin, 1, range)))
sigs = all_sigs[,colnames(x)]




pdf("../results/08_SIS_diag_dsd.pdf", height = 10, width = 10)

library(energy)
library(RColorBrewer)
library(apcluster)
library(fastICA)
library(fdrtool)
library(NMF)
library(MASS)


# From the tables in the CPSS pub, if we have theta = 0.01
# (=> floor(nrow(x)*0.01) = 130 vars selected), and tau = 0.50,
# then we expect fewer than nrow(x)*1.31e-4 = 1.7 incorrect
# vars to be chosen.
set.seed(seed)
cpss.sis = CPSS(x, y, SIS.FAST, tau, 50, nsis = floor(theta*nrow(x)))


# We now have screened variables.  These are only from
# marginal screening, but it's not clear how to intelligently
# apply ISIS to this partiular dimension of problem.  SIS
# should still be able to pick out orthogonal signatures.

# Look at the expression of these variables, and try to separate
# into components.

x.sel = x[cpss.sis$sel,]
xlin.sel = xlin[cpss.sis$sel,]

x.sel.std = (x.sel - rowMeans(x.sel))
x.sel.kcor = cor(t(x.sel.std), method = "kendall")
x.sel.dcor = sapply(1:(nrow(x.sel.std)-1), function(i) c(rep(NA, i), sapply((i+1):nrow(x.sel.std), function(j) dcor(x.sel.std[i,], x.sel.std[j,]))))
x.sel.dcor = cbind(x.sel.dcor, NA)
diag(x.sel.dcor) = 1
x.sel.dcor[upper.tri(x.sel.dcor)] = t(x.sel.dcor)[upper.tri(x.sel.dcor)]

corPlot(x.sel.kcor, main = sprintf("%s: Correlation Clusters of CPSS-SIS-FAST Probes\nKendall log", title))
corPlot(abs(x.sel.kcor), zlim = c(0, 1), pal = "GnBu", main = sprintf("%s: Correlation Clusters of CPSS-SIS-FAST Probes\nAbsolute Kendall log", title))
corPlot(x.sel.dcor, zlim = c(0, 1), pal = "GnBu", main = sprintf("%s: Correlation Clusters of CPSS-SIS-FAST Probes\ndcor log", title))

x.sel.dcor.apc = apclusterK(s = x.sel.dcor, details = TRUE, K = ncomp)
heatmap(x.sel.dcor.apc, x.sel.dcor)


# Perform ICA to try to separate the genes into modules.  
# FastICA finds S, A so that X ~= SA, and columns of S are
# as independent and non-gaussian as possible.  Here S is
# the gene x module 'signature matrix', and A is the module x
# patient score matrix.

set.seed(1234)
xlin.sel.ica = fastICA(xlin.sel, n.comp = ncomp, row.norm = FALSE)
ica.S = xlin.sel.ica$S
ica.A = xlin.sel.ica$A
colnames(ica.A) = colnames(xlin.sel)
ica.k = apply(ica.S, 2, sd)
ica.S = t(t(ica.S) / ica.k)
ica.A = ica.A * ica.k

tests.linica_surv = apply(xlin.sel.ica$A, 1, function(xc) coxph(y ~ xc))
tests.linica_purity = apply(xlin.sel.ica$A, 1, cor.test, y = samples[colnames(x.sel),]$purity_qpure, method = "kendall")


# Try to correlate the loadings with MSigDB GSVA
# Combine all sigs and do the correlations
msigdb.ica.cors = cor(t(ica.A), t(x.msigdb), method = "kendall")
msigdb.ica.cors.fdr = fdrtool(as.vector(msigdb.ica.cors), statistic = "correlation")

# Ok, now map these back to the correlation matrix
msigdb.ica.cors.lfdr = matrix(msigdb.ica.cors.fdr$lfdr, nrow = nrow(msigdb.ica.cors), ncol = ncol(msigdb.ica.cors))
colnames(msigdb.ica.cors.lfdr) = colnames(msigdb.ica.cors)
rownames(msigdb.ica.cors.lfdr) = rownames(msigdb.ica.cors)

pairs(cbind(t(ica.A), t(x.msigdb[colnames(msigdb.ica.cors.lfdr)[apply(msigdb.ica.cors.lfdr, 2, min) < 0.2],])))

	list(
		cpss.sis = cpss.sis,
		kcor = x.sel.kcor,
		dcor = x.sel.dcor,
		ica = xlin.sel.ica,
		ica.S = ica.S,
		ica.A = ica.A,
		msigdb.ica.cor = list(cors = msigdb.ica.cors, lfdr = msigdb.ica.cors.lfdr),
		tests = list(ica_surv = tests.linica_surv, ica_purity = tests.linica_purity)
		)
}


sa.diag_dsd = doSignatureAnalysis(x.diag_dsd, y.diag_dsd, x.msigdb.c123467.merged[,colnames(x.diag_dsd)], "Diagnosis to DSD", 			tau = 0.72, theta = 0.05, x0 = 6.34, ncomp = 5)
sa.surg_dsd = doSignatureAnalysis(x.surg_dsd, y.surg_dsd, x.msigdb.c123467.merged[,colnames(x.surg_dsd)], "Surgery to DSD", 			tau = 0.72, theta = 0.05, x0 = 6.34, ncomp = 5)
sa.recr_dsd = doSignatureAnalysis(x.recr_dsd, y.recr_dsd, x.msigdb.c123467.merged[,colnames(x.recr_dsd)], "Recurrence to DSD", 			tau = 0.72, theta = 0.05, x0 = 6.34, ncomp = 5)
sa.diag_rec = doSignatureAnalysis(x.diag_rec, y.diag_rec, x.msigdb.c123467.merged[,colnames(x.diag_rec)], "Diagnosis to Recurrence", 	tau = 0.72, theta = 0.05, x0 = 6.34, ncomp = 5)
sa.surg_rec = doSignatureAnalysis(x.surg_rec, y.surg_rec, x.msigdb.c123467.merged[,colnames(x.surg_rec)], "Surgery to Recurrence", 		tau = 0.72, theta = 0.05, x0 = 6.34, ncomp = 5)





xlin.diag_dsd = 2^(x.diag_dsd-6.34)
xlin.diag_dsd.sel = xlin.diag_dsd[sa.diag_dsd$cpss.sis$sel,]


temp.nmf = nmf(
	x = xlin.diag_dsd.sel, 
	rank = 6, 
	method = "snmf/l", 
	seed = "random", nrun = 200, 
	.options = list(verbose = 2, track = TRUE, parallel = FALSE))
plot(temp.nmf)
consensusmap(temp.nmf)
basismap(temp.nmf)
coefmap(temp.nmf)

apply(coef(temp.nmf), 1, function(xc) coxph(y.diag_dsd ~ xc))
apply(coef(temp.nmf), 1, function(xc) cor.test(samps.diag_dsd$purity_qpure, xc, method = "kendall"))
temp.nullfit = coxph(y.diag_dsd ~ 1)
temp.nullresids = residuals(temp.nullfit, type = "martingale")
par(mfrow = c(2, 3))
for (i in 1:nrow(coef(temp.nmf)))
{
	scatter.smooth(temp.nullresids ~ coef(temp.nmf)[i,])
}
par(mfrow = c(1, 1))
par(mfrow = c(2, 3))
for (i in 1:nrow(coef(temp.nmf)))
{
	scatter.smooth(samps.diag_dsd$purity_qpure ~ coef(temp.nmf)[i,])
}
par(mfrow = c(1, 1))
temp.nmf.pvals = data.frame(
	p.surv = apply(coef(temp.nmf), 1, function(xc) pchisq(2*diff(coxph(y.diag_dsd ~ xc)$loglik), df = 1, lower.tail = FALSE)),
	p.pure = apply(coef(temp.nmf), 1, function(xc) cor.test(samps.diag_dsd$purity_qpure, xc, method = "kendall")$p.value)
)
temp.qvals = p.adjust(c(temp.nmf.pvals$p.surv, temp.nmf.pvals$p.pure), "holm")
temp.nmf.pvals$q.surv = temp.qvals[1:(length(temp.qvals)/2)]
temp.nmf.pvals$q.pure = temp.qvals[(length(temp.qvals)/2 + 1):length(temp.qvals)]


library(fdrci)
library(parallel)
set.seed(1234)
B = 2000
stats.perm = mclapply(1:B, function(i) {
	cat(i, " ", sep = "")
	coef.perm = coef(temp.nmf)[,sample.int(ncol(coef(temp.nmf)))]
	stats = 1-abs((cor(t(coef.perm), t(x.msigdb.c46.merged[,colnames(coef(temp.nmf))]), method = "kendall")))
	stats
	}, mc.cores = 2)

stats.obs = 1-abs((cor(t(coef(temp.nmf)), t(x.msigdb.c23467.merged[,colnames(coef(temp.nmf))]), method = "kendall")))

stats.threshold = 0.5
rowSums(stats.obs < 1 - stats.threshold)

stats.fdr = t(sapply(1:nrow(coef(temp.nmf)), function(i) {
	obs = stats.obs[i,]
	perm = lapply(stats.perm, function(p) data.frame(statistic = p[i,]))
	fdr = fdr_od(obs, perm, "statistic", length(obs), thres = 1-stats.threshold)
	if (all(is.na(fdr)))	{ return(rep(NA, 7)) }
	return(fdr)
}))
colnames(stats.fdr) = c("FDR", "ll", "ul", "pi0", "c1", "ro", "vp1")
stats.fdr


msigdb.nmf.cors.pvals = apply(coef(temp.nmf), 1, function(c1) apply(x.msigdb.c46.merged[,colnames(coef(temp.nmf))], 1, function(m1) cor.test(c1, m1, method = "kendall")$p.value))
apply(msigdb.nmf.cors.pvals, 2, p.adjust, method = "BY")


par(mfrow = c(3, 3))
for (i in seq(0, 1, length.out = 9))
{
	temp.x = mvrnorm(1e4, c(0, 0), matrix(c(1, i, i, 1), nrow = 2))
	plot(temp.x[1:100,], main = sprintf("%.3g  rho = %.3g", i, cor(temp.x[,1], temp.x[,2], method = "kendall")), xlim = c(-3, 3), ylim = c(-3, 3))
}


# msigdb.nmf.cors = cor(t(coef(temp.nmf)), t(x.msigdb.c46.merged[,colnames(coef(temp.nmf))]), method = "kendall")

# library(fdrtool)
# msigdb.nmf.cors.lfdr_pvals = apply(msigdb.nmf.cors.pvals, 2, function(p) fdrtool(p, "pvalue")$lfdr)
# msigdb.nmf.cors.lfdr_corrs = apply(msigdb.nmf.cors, 1, function(s) fdrtool(s, "correlation")$lfdr)

# apply(msigdb.nmf.cors.lfdr_pvals < 0.2, 2, sum)
# apply(msigdb.nmf.cors.lfdr_corrs < 0.2, 2, sum)


# P values for these correlation statistics are a bit silly.  Perhaps just set
# an effect size limit and take it from there?
#msigdb.nmf.cors.fdr = fdrtool(as.vector(msigdb.nmf.cors), statistic = "correlation")
# Ok, now map these back to the correlation matrix
# msigdb.nmf.cors.lfdr = matrix(msigdb.nmf.cors.fdr$lfdr, nrow = nrow(msigdb.nmf.cors), ncol = ncol(msigdb.nmf.cors))
# colnames(msigdb.nmf.cors.lfdr) = colnames(msigdb.nmf.cors)
# rownames(msigdb.nmf.cors.lfdr) = rownames(msigdb.nmf.cors)


# library(CCA)
# correl = matcor(t(coef(temp.nmf)), t(x.msigdb.merged[,colnames(coef(temp.nmf))]))
# img.matcor(correl, type = 2)
# res.regul <- estim.regul(t(coef(temp.nmf)), t(x.msigdb.merged[,colnames(coef(temp.nmf))]), plt = TRUE)


# library(fdrci)
# library(parallel)
# stats.obs = 1-abs(as.vector(cor(t(coef(temp.nmf)), t(x.msigdb.merged[,colnames(coef(temp.nmf))]), method = "kendall")))
# stats.perm = mclapply(1:1e2, function(i) {
# 	cat(i, " ", sep = "")
# 	coef.perm = coef(temp.nmf)[,sample.int(ncol(coef(temp.nmf)))]
# 	data.frame(statistic = 1-abs(as.vector(cor(t(coef(temp.nmf)), t(x.msigdb.merged[,colnames(coef(temp.nmf))]), method = "kendall"))))
# 	}, mc.cores = 16)


pairs(cbind(t(coef(temp.nmf)), t(x.msigdb.merged[,colnames(coef(temp.nmf))])[,apply(abs(msigdb.nmf.cors), 2, max) >= 0.5]))

pairs(cbind(t(coef(temp.nmf)), t(x.msigdb[,colnames(coef(temp.nmf))])[,rank(-apply(abs(msigdb.nmf.cors), 2, max)) <= 10]))
pairs(cbind(log2(t(coef(temp.nmf)) + 1e-3), t(x.msigdb[,colnames(coef(temp.nmf))])[,rank(-apply(abs(msigdb.nmf.cors), 2, max)) <= 10]))

#library(bnlearn)
#x.sel.bn = boot.strength(as.data.frame(t(x.sel)), algorithm = "hc", debug = TRUE)
#graphviz.plot(x.sel.bn, layout = "dot")


dev.off()

sessionInfo()
