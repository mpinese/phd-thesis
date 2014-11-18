#!~/bin/Rscript

# 03_combine_gex_cpv.R

options(echo = TRUE)

load("../data/03_gex_cpv.rda")

pdf("../results/04_norm_gex.pdf", height = 10, width = 10)

library(lumi)

qc.detrate = colMeans(detection(probe_raw) < 0.01)
qc.controldat = controlData(probe_raw)
qc.negmed = log2(apply(qc.controldat[qc.controldat[,1] == "negative", c(-1,-2)], 2, median))
qc.hkgmed = log2(apply(qc.controldat[qc.controldat[,1] == "housekeeping", c(-1,-2)], 2, median))
qc.beadnum.q10 = apply(beadNum(probe_raw), 2, quantile, prob = 0.1)
qc.cv = log(se.exprs(probe_raw)) / log(exprs(probe_raw))
qc.det = detection(probe_raw)
qc.detected.cvmed = sapply(1:ncol(probe_raw), function(i) {
	det = qc.det[,i]
	cv = qc.cv[,i]
	median(cv[det < 0.01])
})

pairs(cbind("Detection rate" = qc.detrate, "Negative signal" = qc.negmed, "HKG signal" = qc.hkgmed, "Median CV of detected" = qc.detected.cvmed), col = sample_data$scanner_group + 1)

# Looks mostly OK.  Put a threshold on detection rate of 0.30.

gex.qc = probe_raw[,which(qc.detrate >= 0.30)]
sample_data = sample_data[sampleNames(gex.qc),]

# Normalize these qc-ed arrays
gex.qc.b = lumiB(gex.qc, method = "bgAdjust.affy")
gex.qc.bt = lumiT(gex.qc.b, method = "vst")
gex.qc.btn = lumiN(gex.qc.bt, method = "quantile")

# Produce some basic diagnostic plots.
qc.raw.idr = apply(log2(exprs(gex.qc)), 1, quantile, prob = c(0.1, 0.9))
qc.raw.idr = qc.raw.idr[2,] - qc.raw.idr[1,]
qc.norm.idr = apply(exprs(gex.qc.btn), 1, quantile, prob = c(0.1, 0.9))
qc.norm.idr = qc.norm.idr[2,] - qc.norm.idr[1,]
qc.raw.med = apply(log2(exprs(gex.qc)), 1, median)
qc.norm.med = apply(log2(exprs(gex.qc.btn)), 1, median)
qc.detrate2 = rowMeans(detection(gex.qc) < 0.01)
plot(qc.raw.idr ~ qc.raw.med, pch = ".", col = c("lightblue", "black")[I(qc.detrate2 >= 0.1) + 1], main = "Pre-normalization", xlab = "Median expression (log2)", ylab = "Expression IDR (log2)")
legend("topright", legend = c("Detected < 10% of samples", "Detected in >= 10% of samples"), fill = c("lightblue", "black"), inset = 0.05)
plot(qc.norm.idr ~ qc.norm.med, pch = ".", col = c("lightblue", "black")[I(qc.detrate2 >= 0.1) + 1], main = "rma-vst-quantile", xlab = "Median expression (vst)", ylab = "Expression IDR (vst)")
legend("topright", legend = c("Detected < 10% of samples", "Detected in >= 10% of samples"), fill = c("lightblue", "black"), inset = 0.05)

# We have duplicated samples.  Look at these to get an idea of reproducibility
# and technical variation.
temp.array_counts = table(sample_data$patient_id)
temp.array_counts = temp.array_counts[temp.array_counts > 1]
sample_data[sample_data$patient_id %in% names(temp.array_counts),]

#              patient_id scanner_group
# 5445316014_A  APGI_1959  TRUE
# 5491021082_B  APGI_1959  TRUE
# 5624980044_B  APGI_2222  TRUE
# 5624980044_C  APGI_2222  TRUE
# 8986311029_L  APGI_2315 FALSE
# 9020374003_H  APGI_2315 FALSE
# 7976997017_C  APGI_2340 FALSE
# 9020374009_B  APGI_2340 FALSE

# APGI_2222 is a good test for intra-slide variation.  The others
# will display inter- and intra-slide variability.
plot(log2(exprs(gex.qc)[,c("5624980044_B", "5624980044_C")]), pch = 16, cex = 0.5, col = rgb(I(qc.detrate2 < 0.1)*0.7, I(qc.detrate2 < 0.1)*0.7, qc.detrate2 < 0.1, 0.1), main = "APGI_2222 (H) Pre-normalization")
abline(0, 1, col = rgb(1, 0, 0, 0.7), lwd = 2)
plot(log2(exprs(gex.qc.btn)[,c("5624980044_B", "5624980044_C")]), pch = 16, cex = 0.5, col = rgb(I(qc.detrate2 < 0.1)*0.7, I(qc.detrate2 < 0.1)*0.7, qc.detrate2 < 0.1, 0.1), main = "APGI_2222 (H) rma-vst-quantile")
abline(0, 1, col = rgb(1, 0, 0, 0.7), lwd = 2)

plot(log2(exprs(gex.qc)[,c("5445316014_A", "5491021082_B")]), pch = 16, cex = 0.5, col = rgb(I(qc.detrate2 < 0.1)*0.7, I(qc.detrate2 < 0.1)*0.7, qc.detrate2 < 0.1, 0.1), main = "APGI_1959 (H) Pre-normalization")
abline(0, 1, col = rgb(1, 0, 0, 0.7), lwd = 2)
plot(log2(exprs(gex.qc.btn)[,c("5445316014_A", "5491021082_B")]), pch = 16, cex = 0.5, col = rgb(I(qc.detrate2 < 0.1)*0.7, I(qc.detrate2 < 0.1)*0.7, qc.detrate2 < 0.1, 0.1), main = "APGI_1959 (H) rma-vst-quantile")
abline(0, 1, col = rgb(1, 0, 0, 0.7), lwd = 2)

plot(log2(exprs(gex.qc)[,c("8986311029_L", "9020374003_H")]), pch = 16, cex = 0.5, col = rgb(I(qc.detrate2 < 0.1)*0.7, I(qc.detrate2 < 0.1)*0.7, qc.detrate2 < 0.1, 0.1), main = "APGI_2315 (L) Pre-normalization")
abline(0, 1, col = rgb(1, 0, 0, 0.7), lwd = 2)
plot(log2(exprs(gex.qc.btn)[,c("8986311029_L", "9020374003_H")]), pch = 16, cex = 0.5, col = rgb(I(qc.detrate2 < 0.1)*0.7, I(qc.detrate2 < 0.1)*0.7, qc.detrate2 < 0.1, 0.1), main = "APGI_2315 (L) rma-vst-quantile")
abline(0, 1, col = rgb(1, 0, 0, 0.7), lwd = 2)

plot(log2(exprs(gex.qc)[,c("7976997017_C", "9020374009_B")]), pch = 16, cex = 0.5, col = rgb(I(qc.detrate2 < 0.1)*0.7, I(qc.detrate2 < 0.1)*0.7, qc.detrate2 < 0.1, 0.1), main = "APGI_2340 (L) Pre-normalization")
abline(0, 1, col = rgb(1, 0, 0, 0.7), lwd = 2)
plot(log2(exprs(gex.qc.btn)[,c("7976997017_C", "9020374009_B")]), pch = 16, cex = 0.5, col = rgb(I(qc.detrate2 < 0.1)*0.7, I(qc.detrate2 < 0.1)*0.7, qc.detrate2 < 0.1, 0.1), main = "APGI_2340 (L) rma-vst-quantile")
abline(0, 1, col = rgb(1, 0, 0, 0.7), lwd = 2)


# Try and use these plots to get an idea of noise levels.
noisePlot = function(xy, span, ...)
{
	muhat = rowMeans(xy)
	sdhat = apply(xy, 1, sd)
	plot(sdhat ~ muhat, xlab = "Mean expression value", ylab = "Noise SD", ...)
	sdfit = loess(sdhat ~ muhat, family = "symmetric", span = span)
	plot.x = seq(min(xy), max(xy), length.out = 100)
	lines(predict(sdfit, newdata = plot.x) ~ plot.x, lty = "solid", lwd = 2, col = "red")
	legend("topright", legend = "LOESS fit", lty = "solid", lwd = 2, col = "red", inset = 0.05)
}

noisePlot(log2(exprs(gex.qc)[,c("5624980044_B", "5624980044_C")]), span = 0.1, pch = 16, cex = 0.5, ylim = c(0, 1.5), main = "APGI_2222 (H) Pre-normalization",  col = rgb(I(qc.detrate2 < 0.1)*0.7, I(qc.detrate2 < 0.1)*0.7, qc.detrate2 < 0.1, 0.1))
noisePlot(log2(exprs(gex.qc.btn)[,c("5624980044_B", "5624980044_C")]), span = 0.1, pch = 16, cex = 0.5, ylim = c(0, 0.2), main = "APGI_2222 (H) rma-vst-quantile",  col = rgb(I(qc.detrate2 < 0.1)*0.7, I(qc.detrate2 < 0.1)*0.7, qc.detrate2 < 0.1, 0.1))

noisePlot(log2(exprs(gex.qc)[,c("5445316014_A", "5491021082_B")]), span = 0.1, pch = 16, cex = 0.5, ylim = c(0, 1.5), main = "APGI_1959 (H) Pre-normalization", col = rgb(I(qc.detrate2 < 0.1)*0.7, I(qc.detrate2 < 0.1)*0.7, qc.detrate2 < 0.1, 0.1))
noisePlot(log2(exprs(gex.qc.btn)[,c("5445316014_A", "5491021082_B")]), span = 0.1, pch = 16, cex = 0.5, ylim = c(0, 0.2), main = "APGI_1959 (H) rma-vst-quantile", col = rgb(I(qc.detrate2 < 0.1)*0.7, I(qc.detrate2 < 0.1)*0.7, qc.detrate2 < 0.1, 0.1))

noisePlot(log2(exprs(gex.qc)[,c("8986311029_L", "9020374003_H")]), span = 0.1, pch = 16, cex = 0.5, ylim = c(0, 1.5), main = "APGI_2315 (L) Pre-normalization", col = rgb(I(qc.detrate2 < 0.1)*0.7, I(qc.detrate2 < 0.1)*0.7, qc.detrate2 < 0.1, 0.1))
noisePlot(log2(exprs(gex.qc.btn)[,c("8986311029_L", "9020374003_H")]), span = 0.1, pch = 16, cex = 0.5, ylim = c(0, 0.2), main = "APGI_2315 (L) rma-vst-quantile", col = rgb(I(qc.detrate2 < 0.1)*0.7, I(qc.detrate2 < 0.1)*0.7, qc.detrate2 < 0.1, 0.1))

noisePlot(log2(exprs(gex.qc)[,c("7976997017_C", "9020374009_B")]), span = 0.1, pch = 16, cex = 0.5, ylim = c(0, 1.5), main = "APGI_2340 (L) Pre-normalization", col = rgb(I(qc.detrate2 < 0.1)*0.7, I(qc.detrate2 < 0.1)*0.7, qc.detrate2 < 0.1, 0.1))
noisePlot(log2(exprs(gex.qc.btn)[,c("7976997017_C", "9020374009_B")]), span = 0.1, pch = 16, cex = 0.5, ylim = c(0, 0.2), main = "APGI_2340 (L) rma-vst-quantile", col = rgb(I(qc.detrate2 < 0.1)*0.7, I(qc.detrate2 < 0.1)*0.7, qc.detrate2 < 0.1, 0.1))

# Looks like the norm is doing its job.  We can infer something useful about
# expected minimum SD -- it takes on a lowest inter-array value of around 0.015,
# so going by the rule-of-thumb of 2*sd for around 5% false positives, set
# min SD to 0.03.

library("illuminaHumanv4.db")

sel.probe.unsup.det = qc.detrate2 >= 0.1
sel.probe.unsup.sdabovenoise = apply(exprs(gex.qc.btn), 1, sd) >= 0.03
sel.probe.unsup.quality = grepl("Perfect|Good", unlist(mget(featureNames(gex.qc.btn), illuminaHumanv4PROBEQUALITY)))
sel.probe.unsup = sel.probe.unsup.det & sel.probe.unsup.sdabovenoise & sel.probe.unsup.quality

features = data.frame(
	probe_id = featureNames(gex.qc.btn),
	loc = unlist(mget(featureNames(gex.qc.btn), illuminaHumanv4GENOMICLOCATION)),
	entrezid = unlist(mget(featureNames(gex.qc.btn), illuminaHumanv4ENTREZID)),
	symbol = unlist(mget(featureNames(gex.qc.btn), illuminaHumanv4SYMBOL)),
	name = unlist(mget(featureNames(gex.qc.btn), illuminaHumanv4GENENAME)),
	quality = unlist(mget(featureNames(gex.qc.btn), illuminaHumanv4PROBEQUALITY))
	)
rownames(features) = features$probe_id

library(VennDiagram)
draw.triple.venn(area1 = sum(sel.probe.unsup.det), area2 = sum(sel.probe.unsup.sdabovenoise), area3 = sum(sel.probe.unsup.quality),
	n12 = sum(sel.probe.unsup.det & sel.probe.unsup.sdabovenoise), n23 = sum(sel.probe.unsup.sdabovenoise & sel.probe.unsup.quality), n13 = sum(sel.probe.unsup.det & sel.probe.unsup.quality),
	n123 = sum(sel.probe.unsup.det & sel.probe.unsup.sdabovenoise & sel.probe.unsup.quality),
	category = c("Detected in >= 10%", "SD >= 2*noise", "Quality perfect or good"), overrideTriple = TRUE)
library(venneuler)
plot(venneuler(cbind(sel.probe.unsup.det, sel.probe.unsup.sdabovenoise, sel.probe.unsup.quality)))

temp = ls()
rm(list = temp[grepl("^qc\\.", temp)])
rm(list = temp[grepl("^temp\\.", temp)])
rm(probe_raw, temp)
rm(noisePlot)

gex = exprs(gex.qc.btn)
rm(gex.qc.b, gex.qc.bt, gex.qc, gex.qc.btn)

save.image("../data/04_gexnormed_cpv.rda")

dev.off()

sessionInfo()
