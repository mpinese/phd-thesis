load("../../data/07_data_for_SIS.rda")
source("../../common/08_SIS_common_funcs.R")

x = x.diag_dsd
y = y.diag_dsd
seed = 1234567890

library(plyr)
library(doMC)
registerDoMC(32)

grid.tau = seq(0.6, 0.99, 0.01)
grid.gamma = seq(5, 15, 1)

grid.tasks = expand.grid(tau = grid.tau, gamma = grid.gamma)


x.std = (x - rowMeans(x)) / apply(x, 1, sd)

set.seed(seed)
results = aaply(grid.tasks, 1, function(params) {
	print(params)
	nonperm = CPSS(x.std, y, SIS.FAST, params$tau, 50, gamma = params$gamma, scale = FALSE)$sel
	# nonperm_boot = sapply(1:24, function(i) {
	# 	samp = sample.int(ncol(x), replace = TRUE)
	# 	CPSS(x.std[,samp], y[samp,], SIS.FAST, params$tau, 50, gamma = params$gamma, scale = FALSE)$sel
	# })
	perm = sapply(1:25, function(i) {
		CPSS(x.std, y[sample.int(nrow(y)),], SIS.FAST, params$tau, 50, gamma = params$gamma, scale = FALSE)$sel
	})
#	c(sum(nonperm), apply(nonperm_boot, 2, sum), apply(perm, 2, sum))
	c(sum(nonperm), apply(perm, 2, sum))
}, .parallel = TRUE)


nvar.total = aaply(results[,,1], 1:2, median)
nvar.false = aaply(results[,,-1], 1:2, median)
nvar.true = nvar.total - nvar.false
fdr = nvar.false / nvar.total

library(fields)
library(reshape2)
data = melt(nvar.total)
colnames(data)[3] = "nTotal"
data = cbind(data, "nFalse" = melt(nvar.false)[,3])
data$fdr = data$nFalse / data$nTotal
data$nTrue = data$nTotal - data$nFalse
data$nTotalSmth = predict(Tps(x = data[,1:2], Y = data[,3]))
data$nFalseSmth = predict(Tps(x = data[,1:2], Y = data[,4]))
data$fdrSmth = data$nFalseSmth / data$nTotalSmth
data$nTrueSmth = data$nTotalSmth - data$nFalseSmth
nvar.trueSmth = acast(data, tau ~ gamma, value.var = "nTrueSmth")
fdrSmth = acast(data, tau ~ gamma, value.var = "fdrSmth")

pdf("varsel.pdf", height = 10, width = 10)
contour(x = as.numeric(rownames(fdr)), y = as.numeric(colnames(fdr)), z = fdr, zlim = c(0, 0.4), nlevels = 10, col = "red")
contour(x = as.numeric(rownames(nvar.true)), y = as.numeric(colnames(nvar.true)), z = nvar.true, add = TRUE, col = "blue")
contour(x = as.numeric(rownames(fdrSmth)), y = as.numeric(colnames(fdr)), z = fdrSmth, zlim = c(0, 0.4), nlevels = 10, col = "red")
contour(x = as.numeric(rownames(nvar.trueSmth)), y = as.numeric(colnames(nvar.trueSmth)), z = nvar.true, add = TRUE, col = "blue")
plot(fdr ~ nTrue, data)
plot(fdrSmth ~ nTrueSmth, data, ylim = c(0, 1))
dev.off()

save.image("varsel.rda")

