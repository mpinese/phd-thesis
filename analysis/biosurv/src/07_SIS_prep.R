#!~/bin/Rscript

# 07_SIS_prep.R

options(echo = TRUE)

load("../data/05_gexnormed_cpv_matched.rda")


# Subset to as homogeneous patient set as possible.  For now, concentrate
# on well-annotated PDACs from the Australian cohort only, who did not
# suffer from a sugical death.
sel.patient.au = cpvs$Patient.Country == "Australia"
sel.patient.pdac = cpvs$Path.HistoType == "Pancreatic Ductal Adenocarcinoma"
sel.patient.surgok = (cpvs$History.Death.Cause != "Surgical Death") | is.na(cpvs$History.Death.Cause)
sel.patient.diag_dsd = sel.patient.au & sel.patient.pdac & sel.patient.surgok & (!is.na(cpvs$Surv.EventTimeFromDiag.Death)) & (!is.na(cpvs$Surv.Event.DSDeath))
sel.patient.surg_dsd = sel.patient.au & sel.patient.pdac & sel.patient.surgok & (!is.na(cpvs$Surv.EventTimeFromSurg.Death)) & (!is.na(cpvs$Surv.Event.DSDeath))
sel.patient.recr_dsd = sel.patient.au & sel.patient.pdac & sel.patient.surgok & (!is.na(cpvs$Surv.EventTimeFromRec.Death)) & (!is.na(cpvs$Surv.Event.DSDeath))
sel.patient.diag_rec = sel.patient.au & sel.patient.pdac & sel.patient.surgok & (!is.na(cpvs$Surv.EventTimeFromDiag.Rec)) & (!is.na(cpvs$Surv.Event.Rec))
sel.patient.surg_rec = sel.patient.au & sel.patient.pdac & sel.patient.surgok & (!is.na(cpvs$Surv.EventTimeFromSurg.Rec)) & (!is.na(cpvs$Surv.Event.Rec))

sum(sel.patient.diag_dsd)
sum(sel.patient.surg_dsd)
sum(sel.patient.recr_dsd)
sum(sel.patient.diag_rec)
sum(sel.patient.surg_rec)


# Collapse gene expression measurements to one per gene, by 
# selecting the probe with the highest SD for each gene.
x = gex[sel.probe.unsup & !is.na(features$symbol),]
feats = features[sel.probe.unsup & !is.na(features$symbol),]
temp.xperm = order(-apply(x, 1, sd))
x = x[temp.xperm,]
feats = feats[temp.xperm,]
x = x[!duplicated(feats$symbol),]
feats = feats[!duplicated(feats$symbol),]
stopifnot(rownames(x) == rownames(feats))
rownames(x) = feats$symbol


# Construct x and y objects for subsequent work
library(survival)
x.diag_dsd = x[, sel.patient.diag_dsd]
x.surg_dsd = x[, sel.patient.surg_dsd]
x.recr_dsd = x[, sel.patient.recr_dsd]
x.diag_rec = x[, sel.patient.diag_rec]
x.surg_rec = x[, sel.patient.surg_rec]
feats.diag_dsd = features[rownames(x.diag_dsd),]
feats.surg_dsd = features[rownames(x.surg_dsd),]
feats.recr_dsd = features[rownames(x.recr_dsd),]
feats.diag_rec = features[rownames(x.diag_rec),]
feats.surg_rec = features[rownames(x.surg_rec),]
samps.diag_dsd = samples[colnames(x.diag_dsd),]
samps.surg_dsd = samples[colnames(x.surg_dsd),]
samps.recr_dsd = samples[colnames(x.recr_dsd),]
samps.diag_rec = samples[colnames(x.diag_rec),]
samps.surg_rec = samples[colnames(x.surg_rec),]
cpvs.diag_dsd = cpvs[sel.patient.diag_dsd,]
cpvs.surg_dsd = cpvs[sel.patient.surg_dsd,]
cpvs.recr_dsd = cpvs[sel.patient.recr_dsd,]
cpvs.diag_rec = cpvs[sel.patient.diag_rec,]
cpvs.surg_rec = cpvs[sel.patient.surg_rec,]
y.diag_dsd = Surv(cpvs.diag_dsd$Surv.EventTimeFromDiag.DSDeath, 	cpvs.diag_dsd$Surv.Event.DSDeath)
y.surg_dsd = Surv(cpvs.surg_dsd$Surv.EventTimeFromSurg.DSDeath, 	cpvs.surg_dsd$Surv.Event.DSDeath)
y.recr_dsd = Surv(cpvs.recr_dsd$Surv.EventTimeFromRec.DSDeath, 		cpvs.recr_dsd$Surv.Event.DSDeath)
y.diag_rec = Surv(cpvs.diag_rec$Surv.EventTimeFromDiag.Recurrence, 	cpvs.diag_rec$Surv.Event.Recurrence)
y.surg_rec = Surv(cpvs.surg_rec$Surv.EventTimeFromSurg.Recurrence, 	cpvs.surg_rec$Surv.Event.Recurrence)
stopifnot(samps.diag_dsd$patient_id == cpvs.diag_dsd$Patient.ID)
stopifnot(samps.surg_dsd$patient_id == cpvs.surg_dsd$Patient.ID)
stopifnot(samps.recr_dsd$patient_id == cpvs.recr_dsd$Patient.ID)
stopifnot(samps.diag_rec$patient_id == cpvs.diag_rec$Patient.ID)
stopifnot(samps.surg_rec$patient_id == cpvs.surg_rec$Patient.ID)


calcGSVAScores = function(x, rnaseq = FALSE)
{
	require(GSVA)
	require(GSEABase)
	require("illuminaHumanv4.db")

	msigdb.c1 = getGmt("../data/MSigDB/c1.all.v4.0.symbols.gmt", geneIdType = SymbolIdentifier("illuminaHumanv4"), collectionType = BroadCollection(category = "c1"))
	msigdb.c2 = getGmt("../data/MSigDB/c2.all.v4.0.symbols.gmt", geneIdType = SymbolIdentifier("illuminaHumanv4"), collectionType = BroadCollection(category = "c2"))
	msigdb.c3 = getGmt("../data/MSigDB/c3.all.v4.0.symbols.gmt", geneIdType = SymbolIdentifier("illuminaHumanv4"), collectionType = BroadCollection(category = "c3"))
	msigdb.c4 = getGmt("../data/MSigDB/c4.all.v4.0.symbols.gmt", geneIdType = SymbolIdentifier("illuminaHumanv4"), collectionType = BroadCollection(category = "c4"))
	msigdb.c6 = getGmt("../data/MSigDB/c6.all.v4.0.symbols.gmt", geneIdType = SymbolIdentifier("illuminaHumanv4"), collectionType = BroadCollection(category = "c6"))
	msigdb.c7 = getGmt("../data/MSigDB/c7.all.v4.0.symbols.gmt", geneIdType = SymbolIdentifier("illuminaHumanv4"), collectionType = BroadCollection(category = "c7"))

	# Focus on c4 and c6 for now.  These are the most likely
	# to yield useful hits.  c2 will probably give heaps of stuff, but
	# interpretation and FDR control will be nasty.  Get into c7 only
	# if there's some indication of immunological involvement from
	# the previous tests.  c3 always seems promising, but I've never
	# really got it to work well, and c4/6 capture much of it anyway.

	gsva.c1 = gsva(x, msigdb.c1, abs.ranking = TRUE)
	gsva.c2 = gsva(x, msigdb.c2, abs.ranking = TRUE)
	gsva.c3 = gsva(x, msigdb.c3, abs.ranking = TRUE)
	gsva.c4 = gsva(x, msigdb.c4, abs.ranking = TRUE, rnaseq = rnaseq)
	gsva.c6 = gsva(x, msigdb.c6, abs.ranking = FALSE, rnaseq = rnaseq)
	gsva.c7 = gsva(x, msigdb.c7, abs.ranking = FALSE)

	# Convert the signed scores (UP/DN) into single UP-DN combined scores.
	combineSignedScores = function(scores)
	{
		score_stems = gsub("_(UP|DN)$", "", rownames(scores))
		combined = tapply(1:nrow(scores), score_stems, function(is) {
			if (length(is) == 1)	{ return(scores[is,,drop = FALSE]) }
			this_scores = scores[is,]
			up_scores = grepl("_UP$", rownames(this_scores))
			dn_scores = grepl("_DN$", rownames(this_scores))
			other_scores = !(up_scores | dn_scores)
			if (sum(up_scores) == 0 || sum(dn_scores) == 0)		{ return(scores[is,,drop == FALSE]) }
			combined_up_score = apply(scores[is,][up_scores,,drop = FALSE], 2, median)
			combined_dn_score = apply(scores[is,][dn_scores,,drop = FALSE], 2, median)
			signed_score = combined_up_score - combined_dn_score
			result = matrix(signed_score, nrow = 1)
			rownames(result) = gsub("_UP$", "_SIGNED", rownames(scores)[is][which(up_scores)[1]])
			if (sum(other_scores) != 0)
			{
				result = rbind(result, scores[is,][other_scores,,drop = FALSE])
			}
			return(result)
		})
		combined = do.call(rbind, combined)
		colnames(combined) = colnames(scores)
		combined
	}

	gsva.c6.signed = combineSignedScores(gsva.c6$es.obs)
	gsva.c7.signed = combineSignedScores(gsva.c7$es.obs)

	gsva.c123467 = rbind(gsva.c1$es.obs, gsva.c2$es.obs, gsva.c3$es.obs, gsva.c4$es.obs, gsva.c6.signed, gsva.c7.signed)
	gsva.c46 = rbind(gsva.c4$es.obs, gsva.c6.signed)

	list(c46 = gsva.c46, c123467 = gsva.c123467)
}


mergeSimilarGSVAScores = function(scores, threshold)
{
	scores = scale(t(scores))
	dists = as.dist((1 - cor(scores))/2)
	clust = hclust(dists)
	cuts = cutree(clust, h = 1 - threshold)
	merged = tapply(1:ncol(scores), cuts, function(is) {
		newval = matrix(apply(scores[,is,drop = FALSE], 1, median), nrow = 1)
		rownames(newval) = paste(colnames(scores)[is], collapse = "/")
		newval
		})
	merged = do.call(rbind, merged)
	colnames(merged) = rownames(scores)
	merged
}


x.msigdb.sets = calcGSVAScores(x, rnaseq = FALSE)
x.msigdb.c46 = x.msigdb.sets$c46
x.msigdb.c123467 = x.msigdb.sets$c123467
x.msigdb.c46.merged = mergeSimilarGSVAScores(x.msigdb.c46, threshold = 0.98)
x.msigdb.c123467.merged = mergeSimilarGSVAScores(x.msigdb.c123467, threshold = 0.98)


save.image("../data/07_data_for_SIS.rda")


sessionInfo()
