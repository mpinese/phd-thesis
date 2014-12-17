#!~/bin/Rscript

# 05_match_gex_cpv.R

options(echo = TRUE)

load("04_gexnormed_cpv.rda")

pdf("../results/05_match_gex_cpv.pdf", height = 10, width = 10)

library(lumi)

# The main issue is to combine the duplicate arrays into single ones.
# Either can actually average them, or just choose one arbitrarily.

# First see how they stack up as replicates

temp.sd = apply(gex, 1, sd)
temp.selsd = temp.sd > quantile(temp.sd[sel.probe.unsup], 0.95) & sel.probe.unsup
temp.selgex = gex[temp.selsd,]
temp.selgex = (temp.selgex - rowMeans(temp.selgex)) / apply(temp.selgex, 1, sd)
temp.dists = dist(t(temp.selgex))
temp.mds = cmdscale(temp.dists)
temp.col = rep("grey", ncol(temp.selgex))
temp.reps = table(sample_data$patient_id)
temp.reps = temp.reps[temp.reps > 1]
temp.reps = sample_data[sample_data$patient_id %in% names(temp.reps),]
library(RColorBrewer)
temp.pal = brewer.pal(4, "Paired")
names(temp.pal) = unique(temp.reps$patient_id)
for (temp in 1:length(temp.pal)) { temp.col[sample_data$patient_id == names(temp.pal)[temp]] = temp.pal[temp] }
plot(temp.mds, col = temp.col, pch = 16, main = "Replicate consistency: MDS of top 5% by SD")
legend("bottomright", legend = names(temp.pal), fill = temp.pal, inset = 0.05)

# Replicates are very close.  May as well merge them.
temp.gex.merged = c()
temp.samples.merged = c()
sample_data2 = sample_data
sample_data2$purity_qpure = as.numeric(sample_data2$purity_qpure)
sample_data2$exp_array_id = as.character(sample_data2$exp_array_id)
sample_data2$exp_slide_id = as.character(sample_data2$exp_slide_id)
gex2 = gex
for (temp.patid in unique(sample_data$patient_id)) {
	if (temp.patid %in% temp.reps$patient_id)
	{
		cat(temp.patid, "\n")
		gex2 = cbind(gex2, apply(gex[,rownames(temp.reps)[temp.reps$patient_id == temp.patid]], 1, median))
		sample_data2 = rbind(sample_data2, c(temp.patid, "7:Primary tumour", temp.reps$purity_qpure[temp.reps$patient_id == temp.patid][1],
			paste(temp.reps$exp_array_id[temp.reps$patient_id == temp.patid], collapse = "-"),
			paste(unique(temp.reps$exp_slide_id[temp.reps$patient_id == temp.patid]), collapse = "-"),
			paste(temp.reps$path[temp.reps$patient_id == temp.patid], collapse = "-"),
			temp.reps$scanner_group[temp.reps$patient_id == temp.patid][1]))
		colnames(gex2)[ncol(gex2)] = paste(rownames(temp.reps)[temp.reps$patient_id == temp.patid], collapse = "-")
		rownames(sample_data2)[nrow(sample_data2)] = paste(rownames(temp.reps)[temp.reps$patient_id == temp.patid], collapse = "-")
	}
}

temp = !(colnames(gex2) %in% rownames(temp.reps))
gex2 = gex2[,temp]
sample_data2 = sample_data2[temp,]
sample_data2 = as.data.frame(sample_data2)
sample_data2$scanner_group = as.logical(sample_data2$scanner_group)
sample_data2$purity_qpure = as.numeric(sample_data2$purity_qpure)

sample_data2 = sample_data2[colnames(gex2),]

table(table(cpvs$Patient.ID))
# Two patients have two entries in CPV.  Exclude them.
temp = table(cpvs$Patient.ID)
temp = names(temp)[temp == 1]
cpvs = cpvs[cpvs$Patient.ID %in% temp,]

temp = intersect(cpvs$Patient.ID, sample_data2$patient_id)
rownames(cpvs) = cpvs$Patient.ID
cpvs = cpvs[temp,]
gex2 = gex2[,match(temp, sample_data2$patient_id)]
sample_data2 = sample_data2[match(temp, sample_data2$patient_id),]

# Sanity checks
stopifnot(cpvs$Patient.ID == sample_data2$patient_id)

temp = ls()
rm(list = temp[grepl("^temp", temp)])
gex = gex2
samples = sample_data2
rm(sample_data, sample_data2, gex2)

# Sanity checks
stopifnot(cpvs$Patient.ID == samples$patient_id)
stopifnot(colnames(gex) == rownames(samples))

save.image("05_gexnormed_cpv_matched.rda")

dev.off()

sessionInfo()
