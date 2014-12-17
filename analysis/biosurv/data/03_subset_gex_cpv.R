#!~/bin/Rscript

# 03_combine_gex_cpv.R

options(echo = TRUE)

load("02_gex.rda")
cpvs = readRDS("01_cpvs.rds")

pdf("../results/03_combine_gex_cpv.pdf", height = 10, width = 10)

library(VennDiagram)
draw.pairwise.venn(
	area1 = length(unique(sample_data$patient_id[sample_data$sample_code == "7:Primary tumour"])), 
	area2 = length(unique(cpvs$Patient.ID)), 
	cross.area = length(intersect(cpvs$Patient.ID, sample_data$patient_id[sample_data$sample_code == "7:Primary tumour"])),
	category = c("GEX", "CPV"),
	offset = 1)

library(lumi)

temp.sel = sample_data$sample_code == "7:Primary tumour"
sample_data = sample_data[temp.sel,]
probe_raw = probe_raw[,which(temp.sel)]

temp.sel = sample_data$patient_id %in% cpvs$Patient.ID
sample_data = sample_data[temp.sel,]
probe_raw = probe_raw[,which(temp.sel)]

cpvs = cpvs[cpvs$Patient.ID %in% sample_data$patient_id,]

plot(probe_raw)

temp.cd = controlData(probe_raw)
temp.negmed = apply(temp.cd[temp.cd$controlType == "negative", c(-1,-2)], 2, median)
plot(density(temp.negmed), main = "Negative control probe median signal")
rug(temp.negmed)
abline(v = 100, col = "red")
sample_data$scanner_group = temp.negmed > 100

temp.sample_multiplicity = table(sample_data$patient_id)
temp.sample_multiplicity[temp.sample_multiplicity > 1]
sample_data[sample_data$patient_id %in% names(temp.sample_multiplicity)[temp.sample_multiplicity > 1],]
# Ah FFS the replicates are confounded with scanner settings.  No way to verify 
# if the norm. is working or not.  FUCK.  At least they aren't all on the same
# fucking slide too.

temp = NA 
temp = ls()
rm(list = temp[grepl("^temp", temp)])

save.image("03_gex_cpv.rda")

dev.off()

sessionInfo()
