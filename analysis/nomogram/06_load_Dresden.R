options(stringsAsFactors = FALSE)

# Read original data
cpvs = read.csv("data/clin_data_rev_David15012015_DC.csv")
scores = read.csv("data/Dresden_S100A2_S100A4_Scores.csv")


# Tidy up CPVs
temp.cpvs_translate = c(
	"NumericID" = 								"Dresden.NID",
	"Unique_ID" = 								"Dresden.UID",
	"Single_sample_ID" = 						"Dresden.SSID",
	"Number" = 									"Dresden.NUM",
	"Type" = 									"Dresden.Type",
	"Type2" = 									"Dresden.Type2",
	"Typ2_n" = 									"Dresden.Typ2_n",

	"Sexf2m1" = 								"Patient.Gender",
	"Origin" = 									"Patient.Origin",
	
	"Age_at_OP" = 								"History.Surgery.AgeAtYears",
	"post_OP_Survdays" = 						"History.Death.EventTimeDays",
	"Dead0no1yes" = 							"History.Death.Event",
	"Death.from" = 								"History.Death.Cause",

	"R..Margin." = 								"Treat.Surgery.ExcisionStatus",
	"Adjuvant_Therapy1yes0no" = 				"Treat.Chemo.Adjuvant",

	"G..Grade." = 								"Path.Grade",

	"pT" = 										"Staging.pT",
	"pN" = 										"Staging.pN",
	"pM" = 										"Staging.pM",
	"Vascular.space.inavasion.0.no_1.yes" = 	"Path.Invasion.VS",
	"Perineural.invasion.0.no_1.yes" = 			"Path.Invasion.PN",
	"Tumour.location.0.head_1.tail" = 			"Path.TumourLocation",
	"Tumour.size.maximal.diameter.in.mm" = 		"Path.TumourSizeMm"
)
cpvs = cpvs[,names(temp.cpvs_translate)]
colnames(cpvs) = temp.cpvs_translate[colnames(cpvs)]

cpvs$History.Death.Event = cpvs$History.Death.Event == 1
cpvs$History.Death.Cause[cpvs$History.Death.Cause == ""] = NA
cpvs$History.Death.Cause = factor(cpvs$History.Death.Cause)
cpvs$History.DSDeath.Event = cpvs$History.Death.Event & cpvs$History.Death.Cause == "PaCa"
cpvs$Patient.Gender = factor(c("1" = "M", "2" = "F")[as.character(cpvs$Patient.Gender)])
cpvs$Patient.Origin = factor(cpvs$Patient.Origin)
cpvs$Treat.Chemo.Adjuvant = cpvs$Treat.Chemo.Adjuvant == 1
cpvs$Treat.Surgery.ExcisionStatus = ordered(paste("R", cpvs$Treat.Surgery.ExcisionStatus, sep = ""), levels = c("R0", "R1", "R2"))

# Hard to exactly figure out the Type fields.
# Type seems to contain a lot of information, but the classes are mostly obscure to me.
# Type2 is simplified; major classes are PaCa, N (normal? neuroendocrine?)
# Typ2_n is simplified further, where class 1 is almost always Type2 == PaCa, class 2 is always Type2 == N, and class 3 is mostly other.
# Type2 is almost always the last part of Single_sample_ID.
cpvs$Dresden.Type = factor(cpvs$Dresden.Type)
cpvs$Dresden.Type2 = factor(cpvs$Dresden.Type2)
cpvs$Dresden.Typ2_n = factor(cpvs$Dresden.Typ2_n)

cpvs$Staging.pT = ordered(paste("T", cpvs$Staging.pT, sep = ""), levels = c("Tis", "T1", "T2", "T3", "T4"))
cpvs$Staging.pN = ordered(paste("N", cpvs$Staging.pN, sep = ""), levels = c("N0", "N1", "N2"))			# Looks like an old system is in use
cpvs$Staging.pM = ordered(paste("M", cpvs$Staging.pM, sep = ""), levels = c("M0", "M1"))
cpvs$Path.Invasion.VS = cpvs$Path.Invasion.VS == 1
cpvs$Path.Invasion.PN = cpvs$Path.Invasion.PN == 1
cpvs$Path.TumourLocation = factor(c("0" = "Head", "1" = "Tail")[as.character(cpvs$Path.TumourLocation)])


# Tidy up staining scores
temp.scores_translate = c(
	"Eindeutige.Koordinate" = 					"Dresden.Unique.coordinate",
	"alte.Koordinate" = 						"Dresden.Old.coordinate",
	"ID.Nummer" = 								"Dresden.Numeric.ID",
	"TMA" = 									"Dresden.TMA.ID",
	"Herkunft" = 								"Origin",
	"S100A2.cytoplasmic.I...0.3" =				"Molec.S100A2.CytoInt",
	"S100A2.cytoplasmic.." = 					"Molec.S100A2.CytoPercent",
	"S100A2.score" = 							"Molec.S100A2.DCThresh",
	"S100A4.cytoplasmic......below3.." = 		"Molec.S100A4.CytoAnySignal",
	"S100A4.nuclear......below3.." = 			"Molec.S100A4.NucAnySignal",
	"S100A4.score" = 							"Molec.S100A4.DCThresh"
)
scores = scores[,names(temp.scores_translate)]
colnames(scores) = temp.scores_translate[colnames(scores)]

scores$Molec.S100A2.CytoInt[scores$Molec.S100A2.CytoInt %in% c("", "e")] = NA
scores$Molec.S100A2.CytoInt = as.numeric(scores$Molec.S100A2.CytoInt)
scores$Molec.S100A2.CytoPercent[scores$Molec.S100A2.CytoPercent %in% c("", "e")] = NA
scores$Molec.S100A2.CytoPercent = as.numeric(scores$Molec.S100A2.CytoPercent)
scores$Molec.S100A2.DCThresh = (scores$Molec.S100A2.CytoInt >= 3) & (scores$Molec.S100A2.CytoPercent >= 30)
scores$Molec.S100A4.CytoAnySignal[scores$Molec.S100A4.CytoAnySignal %in% c("", "e")] = NA
scores$Molec.S100A4.CytoAnySignal = scores$Molec.S100A4.CytoAnySignal == "p"
scores$Molec.S100A4.NucAnySignal[scores$Molec.S100A4.NucAnySignal %in% c("", "e")] = NA
scores$Molec.S100A4.NucAnySignal = scores$Molec.S100A4.NucAnySignal == "p"
scores$Molec.S100A4.DCThresh = (is.na(scores$Molec.S100A4.CytoAnySignal) | scores$Molec.S100A4.CytoAnySignal) | (is.na(scores$Molec.S100A4.NucAnySignal) | scores$Molec.S100A4.NucAnySignal)


# Combine CPVs and staining scores.  If multiple core scores are available
# for each patient, take the worst score (most positive) for each of the two markers.

# Links between data:
# scores$Dresden.Unique.coordinate <-> cpvs$Dresden.UID
# scores$Origin <-> cpvs$Patient.Origin    (Dresden == 3, Regensburg == 2, Jena == 1)
# scores$Dresden.Numeric.ID <-> cpvs$Dresden.NUM   (Not perfect)

# From DC:
#  Hi Mark,
# 
# Took me a bit of time to understand what the Dresden cohort is all about. Unfortunately 
# can you please do some clean ups and a bit of upfront processing before you can plug into 
# your scripts.
#
# Please find attached 2 excel files.:
# 
#     S100A2/S100A4 score: the ID number is the column that you will need to match to the 
#       file below with clinico-path outcome data. It is marked in red. As for S100A2 and 
#       S100A4 IHC score column marked in yellow, you will need to generate from the raw 
#       scores that were given by Dresden. For S100A2, if it is greater than 3+ at 30%, 
#       then its positive. For S100A4, its either nuclear or cyotplasmic positive, then 
#       its classified as positive. There are usually 5 cores per patients/ID number, so 
#       please take the highest score as the final score (as if the patient has bad prognosis, 
#       we want to know about it). The code for the score are as the second sheet of the file.
#     Clin data: to match up the the S100A2/S10014 scores in the S100A2/S100A4 score, you 
#       will need to match the numeric in the singlesampleID to the ID number of the file 
#       above for the patients with unique ID starting at DD. For patients starting with R 
#       and J, you will need to map back to "Eindeutige Koordinate” column in the 
#       S100A2/S100A4 score file.
# 
# Please not include any M1 patients (metastatic) patients in your nomogram analysis.
# I think this should work, its a bit hard without going through this myself, but when I 
# looked at this again today, I just could not get this done by tonight to get it to you.
# Let me know if this works.
# APGI to follow.
# Cheers
# D

# Mark's notes:
#   The cpv Unique_ID (recoded to Dresden.UID) has the following form:
#     cc_x_y_z
#   cc is the cohort: DD == Dresden, J == Jena, R == Regensburg
#   x is the TMA ID (recoded to Dresden.TMA.ID)
#   y and z seem to form a TMA coordinate, though it's different (transposed?) from that
#     in Dresden.Old.coordinate
#
# 

# I think a simpler strategy than that proposed by DC in the email will work -- merge
# by unique coordinate, *then* collapse down to patients.

# FOLLOWING CODE IS BUGGED, USE UN-COMMENTED DRESDEN-ONLY FOR NOW.
# merged = merge(cpvs, scores, by.x = "Dresden.UID", by.y = "Dresden.Unique.coordinate", all.x = TRUE, all.y = TRUE)
# table(merged$Patient.Origin, merged$Origin)
# merged$Patient.Origin = merged$Origin
# merged$Origin <- NULL

# merged$Dresden.UID = factor(merged$Dresden.UID)
# merged$Dresden.NID = factor(merged$Dresden.NID)
# merged$Dresden.SSID = factor(merged$Dresden.SSID)
# merged$Dresden.NUM = factor(merged$Dresden.NUM)
# merged$Patient.Origin = factor(merged$Patient.Origin)
# merged$Path.Grade = ordered(merged$Path.Grade, levels = 1:4)
# merged$Dresden.Old.coordinate[merged$Dresden.Old.coordinate == ""] = NA
# merged$Dresden.Old.coordinate = factor(merged$Dresden.Old.coordinate)
# merged$Dresden.Numeric.ID[merged$Dresden.Numeric.ID %in% c("", "leer", "unk")] = NA
# merged$Dresden.Numeric.ID = factor(merged$Dresden.Numeric.ID)
# merged$Dresden.TMA.ID = factor(merged$Dresden.TMA.ID)
# merged$Molec.S100A2.CytoInt = ordered(merged$Molec.S100A2.CytoInt, levels = 0:3)


# Remove all non-cancer samples, as it doesn't make sense to merge across disease states.
# Do this quite stringently, as I'm not sure exactly how the Dresden samples are coded, 
# and what data to really use.
# merged = merged[merged$Dresden.Typ2_n == "1" & !is.na(merged$Dresden.Typ2_n) & merged$Dresden.Type2 == "PaCa" & !is.na(merged$Dresden.Type2) & gsub(".*_", "", merged$Dresden.SSID) == "PaCa",]


# # From my reading of DC's email, patient ID is in the scores "Single_sample_ID" = "Dresden.SSID"
# # field.  Collapse down based on this.
# vector_identical = function(x)
# {
# 	length(unique(x)) == 1
# }

# merge_function = function(this)
# {
# 	if (nrow(this) == 1)	{ return(this) }

# 	this_cpv = this[,!grepl("(Dresden\\.UID)|(Molec\\.)|(Dresden\\.Old\\.coordinate)|(Dresden\\.NID)|(Dresden\\.Numeric\\.ID)", colnames(this))]
# 	this_scores = this[,c("Molec.S100A2.DCThresh", "Molec.S100A4.DCThresh")]
# 	cpvs_concordant = all(apply(this_cpv, 2, function(col) vector_identical(col)))
# 	if (!cpvs_concordant) { print(this_cpv) }

# 	scores_collapsed = data.frame(apply(this_scores, 2, function(x) {
# 		if (all(is.na(x))) { return(NA) }
# 		return(max(x, na.rm = TRUE) == 1)
# 	}))

# 	na_count = apply(is.na(this_cpv), 1, sum)
# 	this_single = this[which.min(na_count),]

# 	this_single$Molec.S100A2.DCThresh = scores_collapsed$Molec.S100A2.DCThresh
# 	this_single$Molec.S100A4.DCThresh = scores_collapsed$Molec.S100A4.DCThresh
# 	return(this_single)
# }

# library(plyr)

# merged_collapsed = ddply(merged, ~ Dresden.SSID, merge_function)


# #write.csv(merged_collapsed, "data/Dresden_20150215.csv")
# saveRDS(merged_collapsed, "06_Dresden.rds")




library(plyr)

cpvs = cpvs[cpvs$Patient.Origin == "3",]
scores = scores[scores$Origin == "Dresden",]
merged = merge(cpvs, scores, by.x = "Dresden.UID", by.y = "Dresden.Unique.coordinate", all.x = TRUE, all.y = FALSE)
temp.a2 = tapply(1:nrow(merged), merged$Dresden.SSID, function(is) max(merged$Molec.S100A2.DCThresh[is], na.rm = TRUE))
temp.a2[temp.a2 == -Inf] = NA
temp.a2 = temp.a2 == 1
temp.a4 = tapply(1:nrow(merged), merged$Dresden.SSID, function(is) max(merged$Molec.S100A4.DCThresh[is], na.rm = TRUE))
temp.a4[temp.a4 == -Inf] = NA
temp.a4 = temp.a4 == 1
merged$Molec.S100A2.DCThresh = temp.a2[merged$Dresden.SSID]
merged$Molec.S100A4.DCThresh = temp.a4[merged$Dresden.SSID]
merged = merged[,c("Dresden.SSID", "Dresden.Type", "Dresden.Type2", "Patient.Gender", "History.Surgery.AgeAtYears", "History.Death.EventTimeDays", "History.Death.Event", "History.DSDeath.Event", "History.Death.Cause", "Treat.Surgery.ExcisionStatus", "Path.Grade", "Staging.pT", "Staging.pN", "Staging.pM", "Path.Invasion.VS", "Path.Invasion.PN", "Path.TumourLocation", "Path.TumourSizeMm", "Molec.S100A2.DCThresh", "Molec.S100A4.DCThresh")]
merged = merged[merged$Dresden.Type2 == "PaCa",]
merged = merged[!is.na(merged$History.DSDeath.Event) & !is.na(merged$History.Death.EventTimeDays),]
merged = merged[,!(colnames(merged) %in% c("Dresden.Type", "Dresden.Type2"))]
merged$Dresden.SSID = factor(merged$Dresden.SSID)
merged = merged[!is.na(merged$Path.TumourLocation) & !is.na(merged$Path.TumourSizeMm) & !is.na(merged$Molec.S100A2.DCThresh) & !is.na(merged$Molec.S100A4.DCThresh) & merged$Path.TumourSizeMm > 0,]
merged$Path.Grade = ordered(merged$Path.Grade, levels = 1:4)
modal_value <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}
merged2 = ddply(merged, ~ Dresden.SSID, function(x) apply(x, 2, modal_value))
for (i in 1:ncol(merged))
{
	if ("ordered" %in% class(merged[,i]))		{ merged2[,i] = ordered(merged2[,i]) }
	else if (class(merged[,i]) == "factor")		{ merged2[,i] = factor(merged2[,i]) }
	else if (class(merged[,i]) == "integer")	{ merged2[,i] = as.integer(merged2[,i]) }
	else if (class(merged[,i]) == "logical")	{ merged2[,i] = as.logical(merged2[,i]) }
	else if (class(merged[,i]) == "array")		{ merged2[,i] = as.logical(merged2[,i]) }
}
merged2 = merged2[merged2$Staging.pM != "M1" & !is.na(merged2$Staging.pM),]

saveRDS(merged2, "06_Dresden.rds")
write.csv(merged2, "data/Dresden_20150215.csv")
