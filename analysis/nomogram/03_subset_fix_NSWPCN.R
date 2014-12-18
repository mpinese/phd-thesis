load("01_NSWPCN.rda")

pdf("03_subset_NSWPCN.pdf")

data$Molec.S100A4.DCThresh = (data$Molec.S100A4.NucInt >= "1") | (data$Molec.S100A4.CytoInt >= "1")
data$Molec.S100A2.DCThresh = (data$Molec.S100A2.Int >= "3") & (data$Molec.S100A2.Percent > 30)

temp.sel = data$Treat.Resected & 
	!is.na(data$History.Death.EventTimeDays) & 
	!is.na(data$History.DSDeath.Event) & 
	data$History.Death.EventTimeDays > 14 & 
	data$Cohort.ICGC == FALSE & 
	data$History.Diagnosis.Date >= as.Date("1994-01-01") &		# Data become sparse pre-1994, making correction for date noisy
	!is.na(data$History.Diagnosis.Date)
data = data[temp.sel,]
rownames(data) = paste("NSWPCN_", data$Patient.ID, sep = "")

data$Stage.pT.Simplified = ordered(c("T1" = "T1", "T2" = "T2", "T3" = "T34", "T4" = "T34")[as.character(data$Stage.pT)], levels = c("T1", "T2", "T34"))

median(data$History.Diagnosis.AgeAt, na.rm = TRUE)		# 68
median(data$History.Smoking.PackYears, na.rm = TRUE)	# 30
median(data$Path.Size, na.rm = TRUE)					# 30
median(data$Path.Bilirubin.Preop, na.rm = TRUE)			# 3.45
data$Path.Ca199.Preop.Cent = log(data$Path.Ca199.Preop)
median(data$Path.Ca199.Preop.Cent, na.rm = TRUE)		# 5.38
median(data$Path.Bilirubin.Postop, na.rm = TRUE)		# 0.643
data$Path.Ca199.Postop.Cent = log(data$Path.Ca199.Postop)
median(data$Path.Ca199.Postop.Cent, na.rm = TRUE)		# 3.97
median(data$History.Diagnosis.Date, na.rm = TRUE)		# "2002-01-13"

data$History.Diagnosis.AgeAt.Cent = data$History.Diagnosis.AgeAt - 68
data$History.Smoking.PackYears.Cent = data$History.Smoking.PackYears - 30
data$Path.Size.Cent = data$Path.Size - 30
data$Path.Bilirubin.Preop.Cent = data$Path.Bilirubin.Preop - 3.45
data$Path.Ca199.Preop.Cent = data$Path.Ca199.Preop.Cent - 5.38
data$Path.Bilirubin.Postop.Cent = data$Path.Bilirubin.Postop - 0.643
data$Path.Ca199.Postop.Cent = data$Path.Ca199.Postop.Cent - 3.97
data$History.Diagnosis.Date.Cent = as.numeric(data$History.Diagnosis.Date - as.Date("2002-01-13"))
data$Path.LN.InvolvedFraction = data$Path.LN.Involved / data$Path.LN.Inspected
data$Path.LN.Negative = data$Path.LN.Inspected - data$Path.LN.Involved


data.x.preop.conventional = data[,c(
	"Patient.Sex", 
#	"History.PreviousMalignancy", 			# Not available in Glasgow cohort
#	"History.FdrWithPancCancer", 			# Not available in Glasgow cohort
#	"History.FdrWithAnyCancer", 			# Not available in Glasgow cohort
	"History.Diagnosis.AgeAt.Cent", 
#	"History.AlcoholLevel", 				# Not available in Glasgow cohort.  Surprisingly.
#	"History.Smoking.Status", 				# Not available in Glasgow cohort
#	"History.Smoking.PackYears.Cent", 		# Almost always missing
#	"History.Comorbid.Diabetes", 			# Not available in Glasgow cohort
#	"History.Comorbid.ChronicPancreatitis", # Not available in Glasgow cohort
#	"Treat.ProcedureWhipple", 				# Close to confounded with Path.LocationBody
	"Path.LocationBody", 
	"Path.Size.Cent", 
#	"Path.Bilirubin.Preop.Cent", 			# Tempting but often missing
#	"Path.Ca199.Preop.Cent", 				# Tempting but often missing
	"Stage.pT.Simplified"
)]

data.x.preop.molecular = data[,c(
	"Molec.S100A4.DCThresh", 
	"Molec.S100A2.DCThresh")]

data.x.postop.conventional = data[,c(
	"Treat.MarginPositive", 
#	"Path.Bilirubin.Postop.Cent", 			# Tempting but often missing
#	"Path.Ca199.Postop.Cent", 				# Tempting but often missing
	"Path.Subtype", 
	"Path.Differentiation", 
#	"Path.LN.Involved", 
	"Path.LN.InvolvedFraction",
	"Path.Invasion.Vascular", 
	"Path.Invasion.Perineural", 
#	"Stage.pM", 							# Missing for a lot of pts, difficult to measure confidently 
#	"Stage.Overall",						# Without Stage.M, can't be determined -- M1 trumps all others and => Stage IV
	"Stage.pN")]

data.x.management = data[,c(
	"Treat.Chemo.Adjuvant"
# 	"Treat.Chemo.Adjuvant.GE3Cycles", 		# Not available in Glasgow cohort
#	"Treat.Chemo.Palliative", 				# Often missing, unlikely to be survival-associated
#	"Treat.Chemo.PalliativeDC", 			# Interpretation unknown, unlikely to be survival-associated
#	"Treat.Chemo.Any", 						# Often missing
#	"Treat.Chemo.GEM", 						# Often missing
#	"Treat.Radio"							# Not available in Glasgow cohort
), drop = FALSE]

data.x.confounders = data[,c(
	"History.Diagnosis.Date.Cent"
#	"Cohort.ICGC"							# Enforced by subsetting
	), drop = FALSE]

data.x.extra = data[,c(
	"Path.LN.Inspected",
	"Path.LN.Involved",
	"Path.LN.Negative",
	"History.Diagnosis.AgeAt",
	"Path.Size")]


library(survival)
data.y = Surv(data$History.Death.EventTimeDays, data$History.DSDeath.Event)


library(gplots)
data.x.all = cbind(data.x.preop.conventional, data.x.preop.molecular, data.x.postop.conventional, data.x.management, data.x.confounders, data.x.extra)
heatmap.2(t(1-is.na(data.x.all)), trace = "none", scale = "none", useRaster = TRUE, col = grey(c(0.1, 0.9)))


# Subset to complete data
temp.sel = apply(!is.na(data.x.all), 1, all)
data.x.preop.conventional = data.x.preop.conventional[temp.sel,]
data.x.postop.conventional = data.x.postop.conventional[temp.sel,]
data.x.preop.molecular = data.x.preop.molecular[temp.sel,]
data.x.management = data.x.management[temp.sel,]
data.x.confounders = data.x.confounders[temp.sel,,drop=FALSE]
data.x.all = data.x.all[temp.sel,]
data.y = data.y[temp.sel,]

# Ok so here's how they combine:
# 					preop.conventional	preop.molecular		postop.conventional		management		counfounders
# Conv-postop		YES					NO 					YES 					YES				YES
# Conv-preop 		YES 				NO 					NO  					YES				YES
# Molec-postop		YES					YES					YES						YES				YES
# Molec-preop		YES					YES					NO 						YES				YES
#
# Always model the confounders, but don't include them in the nomograms.  Take
# predictions at some suitable point.
#
# Conv-postop is MSKCC-mimic.  Compare to MSKCC actual.
# Conv-preop is MSKCC-like but preop only.  Compare to MSKCC with marginalized unknowns.
# Molec-postop is We-Know-All.
# Molec-preop is The DC Special
#

data.x.conv_postop = cbind( data.x.preop.conventional,                         data.x.postop.conventional, data.x.management, data.x.confounders)
data.x.conv_preop = cbind(  data.x.preop.conventional,                                                     data.x.management, data.x.confounders)
data.x.molec_postop = cbind(data.x.preop.conventional, data.x.preop.molecular, data.x.postop.conventional, data.x.management, data.x.confounders)
data.x.molec_preop = cbind( data.x.preop.conventional, data.x.preop.molecular,                             data.x.management, data.x.confounders)

rm(temp.sel)


save.image("03_NSWPCN_subset.rda")

sessionInfo()
