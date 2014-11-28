#!~/bin/Rscript

# 15_prep_validation.R -- load validation GEX and CPV data, and save
# as a processed rda.

options(stringsAsFactors = FALSE, echo = TRUE)


######################################################################
# GSE21501
######################################################################

GSE21501.samp = t(read.table(pipe("xz -dc ../data/validation/GSE21501_series_matrix.txt.xz | grep -E '^!Sample'"), sep = "\t"))
colnames(GSE21501.samp) = gsub("^!", "", GSE21501.samp[1,])
GSE21501.samp = data.frame(GSE21501.samp[-1,], check.names = TRUE)
GSE21501.samp = GSE21501.samp[,c("Sample_title", "Sample_geo_accession", "Sample_characteristics_ch2", "Sample_characteristics_ch2.1", "Sample_characteristics_ch2.2", "Sample_characteristics_ch2.3")]
GSE21501.samp = GSE21501.samp[grepl("^os", GSE21501.samp$Sample_characteristics_ch2),]
colnames(GSE21501.samp) = c("Sample_title" = "id", "Sample_geo_accession" = "acc", "Sample_characteristics_ch2" = "time", "Sample_characteristics_ch2.1" = "event", "Sample_characteristics_ch2.2" = "tstage", "Sample_characteristics_ch2.3" = "nstage")[colnames(GSE21501.samp)]
GSE21501.samp$time = as.numeric(gsub("os time: ", "", GSE21501.samp$time))
GSE21501.samp$event = as.numeric(gsub("os event: ", "", GSE21501.samp$event))
GSE21501.samp$tstage = as.numeric(gsub("t stage: ", "", GSE21501.samp$tstage))
GSE21501.samp$nstage = as.numeric(gsub("n stage: ", "", GSE21501.samp$nstage))
GSE21501.samp = GSE21501.samp[order(GSE21501.samp$acc),]
GSE21501.samp$event = GSE21501.samp$event == 1

GSE21501.gex = read.table("../data/validation/GSE21501_series_matrix.txt.xz", sep = "\t", comment.char = "!", header = TRUE)
GSE21501.gex = GSE21501.gex[,GSE21501.samp$acc]
GSE21501.gex = as.matrix(GSE21501.gex)


GSE21501.feat = read.table(pipe("xz -dc ../data/validation/GPL4133.annot.xz | grep -Ev '^[!#^]'"), sep = "\t", comment.char = "", quote = "", header = TRUE)

temp.sd = apply(GSE21501.gex, 1, sd)
temp.perm = order(-temp.sd)
GSE21501.gex = GSE21501.gex[temp.perm,]
GSE21501.feat = GSE21501.feat[temp.perm,]
temp.sel = !duplicated(GSE21501.feat$Gene.symbol) & (GSE21501.feat$Gene.symbol != "")
GSE21501.gex = GSE21501.gex[temp.sel,]
GSE21501.feat = GSE21501.feat[temp.sel,]

# Hard threshold at median SD
temp.sd = apply(GSE21501.gex, 1, sd)
temp.sel = temp.sd > quantile(temp.sd, 0.5)
GSE21501.gex = GSE21501.gex[temp.sel,]
GSE21501.feat = GSE21501.feat[temp.sel,]

# Transform to 0-1 linear
GSE21501.lingex = 2^GSE21501.gex
GSE21501.lingex = (GSE21501.lingex - apply(GSE21501.lingex, 1, min)) / (apply(GSE21501.lingex, 1, max) - apply(GSE21501.lingex, 1, min))





######################################################################
# GSE28735
######################################################################

GSE28735.samp = t(read.table(pipe("xz -dc ../data/validation/GSE28735_series_matrix.txt.xz | grep -E '^!Sample'"), sep = "\t"))
colnames(GSE28735.samp) = gsub("^!", "", GSE28735.samp[1,])
GSE28735.samp = data.frame(GSE28735.samp[-1,], check.names = TRUE)
GSE28735.samp = GSE28735.samp[,c("Sample_title", "Sample_geo_accession", "Sample_characteristics_ch1", "Sample_characteristics_ch1.1", "Sample_characteristics_ch1.2", "Sample_description")]
GSE28735.samp = GSE28735.samp[GSE28735.samp$Sample_characteristics_ch1 == "tissue: T",]
GSE28735.samp = GSE28735.samp[!grepl("na$", GSE28735.samp$Sample_characteristics_ch1.1),]
GSE28735.samp$Sample_title = gsub(".* ", "", GSE28735.samp$Sample_title)
GSE28735.samp$Sample_characteristics_ch1.1 = as.numeric(gsub("survival_month: ", "", GSE28735.samp$Sample_characteristics_ch1.1))
GSE28735.samp$Sample_characteristics_ch1.2 = GSE28735.samp$Sample_characteristics_ch1.2 == "cancer_death: 1"
colnames(GSE28735.samp) = c("Sample_title" = "id", "Sample_geo_accession" = "acc", "Sample_characteristics_ch1.1" = "time", "Sample_characteristics_ch1.2" = "event")[colnames(GSE28735.samp)]
GSE28735.samp = GSE28735.samp[order(GSE28735.samp$acc),]

GSE28735.gex = read.table("../data/validation/GSE28735_series_matrix.txt.xz", sep = "\t", comment.char = "!", header = TRUE)
GSE28735.gex = GSE28735.gex[,GSE28735.samp$acc]
GSE28735.gex = as.matrix(GSE28735.gex)


GSE28735.feat = read.table(pipe("xz -dc ../data/validation/GPL6244.annot.xz | grep -Ev '^[!#^]'"), sep = "\t", comment.char = "", quote = "", header = TRUE)

temp.sd = apply(GSE28735.gex, 1, sd)
temp.perm = order(-temp.sd)
GSE28735.gex = GSE28735.gex[temp.perm,]
GSE28735.feat = GSE28735.feat[temp.perm,]
temp.sel = !duplicated(GSE28735.feat$Gene.symbol) & (GSE28735.feat$Gene.symbol != "")
GSE28735.gex = GSE28735.gex[temp.sel,]
GSE28735.feat = GSE28735.feat[temp.sel,]

# Hard threshold at median SD
temp.sd = apply(GSE28735.gex, 1, sd)
temp.sel = temp.sd > quantile(temp.sd, 0.5)
GSE28735.gex = GSE28735.gex[temp.sel,]
GSE28735.feat = GSE28735.feat[temp.sel,]

# Transform to 0-1 linear
GSE28735.lingex = 2^GSE28735.gex
GSE28735.lingex = (GSE28735.lingex - apply(GSE28735.lingex, 1, min)) / (apply(GSE28735.lingex, 1, max) - apply(GSE28735.lingex, 1, min))



temp = NA
temp = ls()
rm(list = temp[grepl("^temp", temp)])


save.image("../data/15_validation.rda")

sessionInfo()

