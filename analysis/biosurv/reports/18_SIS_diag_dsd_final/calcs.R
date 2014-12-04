######################################################################
# LIBRARIES
######################################################################
options(java.parameters = "-Xmx4G", warn = 1)
library(glmulti)
library(glmnet)
library(energy)
library(NMF)

nmf.options(cores = 32, pbackend = "par", gc = 1, shared.memory = FALSE)


######################################################################
# DATA
######################################################################
load("../../data/07_data_for_SIS.rda")
source("../../common/08_SIS_common_funcs.R")

######################################################################
# HIGH LEVEL PARAMETERS
######################################################################
x = x.diag_dsd
y = y.diag_dsd
sigs = x.diag_dsd.msigdb.merged
samps = samps.diag_dsd
tau = 0.8
gamma = 7
seed = 1234567890
nmf.beta = 0.01
nmf.nrun.rank = 50
nmf.nrun.fit = 500
nmf.rank = "auto"
nmf.rankrange = 2:9
nmf.rankrandcount = 5
sig.corr.threshold = 0.5


######################################################################
# DERIVED VARIABLES
######################################################################
linearizeX = function(x)
{
	xlin = 2^x
	(xlin - apply(xlin, 1, min)) / as.vector(diff(apply(xlin, 1, range)))
}

xlin.diag_dsd = linearizeX(x.diag_dsd)
xlin.diag_rec = linearizeX(x.diag_rec)
xlin.recr_dsd = linearizeX(x.recr_dsd)
xlin.pdac_au = linearizeX(x.pdac_au)
xlin = linearizeX(x)


######################################################################
# PROBE SELECTION
######################################################################
x.std = (x - rowMeans(x)) / apply(x, 1, sd)
set.seed(seed)
cpss.sis = CPSS(x.std, y, SIS.FAST, tau, 50, gamma = gamma, scale = FALSE)
cpss.sis.permuted = simplify2array(mclapply(1:50, function(i) {
	cat(i)
	set.seed(seed + i)
	CPSS(x.std, y[sample.int(nrow(y)),], SIS.FAST, tau, 50, gamma = gamma, scale = FALSE)$sel
}, mc.cores = 32))
sum(cpss.sis$sel)
cpss.sis$qhat / nrow(x.std)
median(apply(cpss.sis.permuted, 2, sum))

x.sel = x[cpss.sis$sel,]
xlin.sel = xlin[cpss.sis$sel,]
xlin.diag_dsd.sel = xlin.diag_dsd[cpss.sis$sel,]
xlin.diag_rec.sel = xlin.diag_rec[cpss.sis$sel,]
xlin.recr_dsd.sel = xlin.recr_dsd[cpss.sis$sel,]
xlin.pdac_au.sel = xlin.pdac_au[cpss.sis$sel,]


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
nmf.runs.rank = nmf(
	x = xlin.sel, 
	rank = nmf.rankrange, 
	method = "snmf/l", 
	seed = seed, nrun = nmf.nrun.rank, 
	.options = list(verbose = 1, track = FALSE, parallel = TRUE, keep.all = FALSE),
	beta = nmf.beta)
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
#temp.perm_resids.delta.threshold = temp.perm_resids.delta.mean*1.1
temp.perm_resids.delta.above_threshold = temp.orig_resids.delta >= temp.perm_resids.delta.threshold
if (all(temp.perm_resids.delta.above_threshold))			{ nmf.rank.auto = min(nmf.rankrange) 
} else if (all(!(temp.perm_resids.delta.above_threshold)))	{ nmf.rank.auto = max(nmf.rankrange)
} else 														{ nmf.rank.auto = min(nmf.rankrange[temp.perm_resids.delta.above_threshold]) }
nmf.rank.wasauto = FALSE
if (nmf.rank == "auto")
{
	nmf.rank = nmf.rank.auto
	nmf.rank.wasauto = TRUE
}


######################################################################
# FACTORIZATION
######################################################################
nmf.final = nmf(
	x = xlin.sel, 
	rank = nmf.rank, 
	method = "snmf/l", 
	seed = seed, nrun = nmf.nrun.fit, 
	.options = list(verbose = 0, track = FALSE, parallel = TRUE, keep.all = FALSE),
	beta = nmf.beta)

######################################################################
# SIGNATURE-METAGENE CORRELATION
######################################################################
cat("Correlation")
nmf.final.msigdb.corr = cor(t(coef(nmf.final)), t(sigs), method = "kendall")


######################################################################
# LASSO
######################################################################
glmnet.x = t(coef(nmf.final))
colnames(glmnet.x) = paste("mg", 1:ncol(glmnet.x), sep = ".")
glmnet.fit.cv = cv.glmnet(x = glmnet.x, y = cbind(time = y[,1], status = y[,2]*1), family = "cox", nfolds = 10)
glmnet.coef.1se = coef(glmnet.fit.cv$glmnet.fit, s = glmnet.fit.cv$lambda.1se)
glmnet.coef.min = coef(glmnet.fit.cv$glmnet.fit, s = glmnet.fit.cv$lambda.min)


######################################################################
# META-PCNA SCORING
######################################################################
metapcna.sig = c("ADAMTS13", "ALAS2", "APOBEC3B", "ARID3A", "ASF1B", "AURKA", 
	"AURKB", "BIRC5", "BPGM", "BUB1B", "BZRPL1", "C21orf45", "CCNA2", "CCNB1", 
	"CCNB2", "CDC2", "CDC20", "CDC45L", "CDCA3", "CDCA4", "CDCA8", "CDKN3", 
	"CDT1", "CENPA", "CHAF1A", "CKLF", "CKS1B", "CKS2", "DDX39", 
	"DKFZp762E1312", "DTL", "EPB42", "ERAF", "ESPL1", "FBXO5", "FBXO7", 
	"FECH", "FEN1", "FOXM1", "GATA1", "GINS1", "GINS2", "GTPBP2", "GTSE1", 
	"GYPA", "GYPB", "H3F3A", "HMBS", "HMGB2", "HMGN2", "KEL", "KIAA0101", 
	"KIF20A", "KIF22", "KIF2C", "KIF4A", "KLF1", "KLF15", "LBR", "LIG1", 
	"LMNB1", "LOC146909", "LSM6", "LYL1", "MAD2L1", "MCM2", "MCM3", "MCM4", 
	"MCM5", "MCM6", "MCM7", "MELK", "MICB", "MKI67", "MLF1IP", "NCAPD2", 
	"NCAPD3", "NCAPG2", "NFE2", "NUDT1", "NUP210", "NUP37", "NUSAP1", "OIP5", 
	"ORC6L", "PCNA", "PF4", "PGD", "PLEK", "POLE2", "PPBP", "PPIH", "PRC1", 
	"PSMD9", "PTTG1", "RACGAP1", "RAD51AP1", "RFC3", "RFC4", "RFWD3", "RHAG", 
	"RHCE", "RHD", "RPA3", "RPIA", "RPP30", "RRM2", "SFRS2", "SHCBP1", "SMC4", 
	"SNF8", "SNRPB", "SNRPD1", "SPTA1", "TACC3", "TAL1", "TCF3", "TFDP1", 
	"TIMELESS", "TOP2A", "TPX2", "TRIM10", "TRIM58", "TRMT5", "TROAP", "TYMS", 
	"UBE2C", "VRK1", "WHSC1", "ZWINT")
metapcna.scores = apply(x[rownames(x) %in% metapcna.sig,], 2, median)


session_info = sessionInfo()
save.image("image.rda")



