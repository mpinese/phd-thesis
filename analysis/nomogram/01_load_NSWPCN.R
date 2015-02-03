options(stringsAsFactors = FALSE)

data = read.table("./data/Total_cohort_MP_05022013.txt.xz", sep = "\t", quote = "", header = TRUE)
data[data == "."] = NA
data[data == ""] = NA
origcols = scan("./data/Total_cohort_MP_05022013.txt.xz", sep = "\t", nlines = 1, what = character())

# summary(data)
# ifelse(sapply(apply(data, 2, table), length) < 20, apply(data, 2, table, useNA = "ifany"), length(apply(data, 2, table)))

temp.fields = list(
	"ID" = 																	"Patient.ID",
	"Sex" = 																"Patient.Sex",
	"ICGC...1" = 															"Cohort.ICGC",
	"Pre.malignancy...1" = 													"History.PreviousMalignancy",
	"FHx.PC.0.no..1.FDR" = 													"History.FdrWithPancCancer",
	"FHx.FDR.Ca...1" = 														"History.FdrWithAnyCancer",
#	"Path.Date.Date.of.diagnosis.surgery...D0." = 							"History.Diagnosis.Date1",
	"Path.Date.Date.of.diagnosis.surgery...D0..2" = 						"History.Diagnosis.Date",
	"Age.at.Diagnosis" = 													"History.Diagnosis.AgeAt",
	"Alcohol.0.never..1.0.2.2.3.4.3.gr.5" = 								"History.AlcoholLevel",
	"Cig.smoking.0.never..1.ceased..2.current" = 							"History.Smoking.Status",
	"Cig.smoking.pack.years" = 												"History.Smoking.PackYears",
	"Comorbid.DM.1" = 														"History.Comorbid.Diabetes",
	"Comorbid.CP.1" = 														"History.Comorbid.ChronicPancreatitis",
	"Recurence.1" = 														"History.Recurrence.Event",
	"Date.of.recurrence" = 													"History.Recurrence.Date",
	"Censor...Cancer.Death.1" = 											"History.DSDeath.Event",
#	"Cause.Death" = 														"History.Death.Cause",
	"Date.of.Death" = 														"History.Death.Date",
	"Date.Last.Followup" = 													"History.Followup.Date",
	"LOF..days." = 															"History.Death.EventTimeDays",
	"Status" = 																"History.Status",
	"CPa.resect.0.Biopsy.2" = 												"Treat.Resected",
	"X..Operation.Whipple.0.." = 											"Treat.ProcedureWhipple",
#	"Margin..mm." = 														"Treat.MarginMm",					# Variable unreliable
	"X..Margins.Pos.1.." = 													"Treat.MarginPositive",
	"X..Any.Chemo.1.." = 													"Treat.Chemo.Any",
	"X..Adjuvant.Chemo.1.." = 												"Treat.Chemo.Adjuvant",
	"X..Adjuvant.Chem....3.cycles.." = 										"Treat.Chemo.Adjuvant.GE3Cycles",
	"X..Palliative.Chemo..1.." = 											"Treat.Chemo.Palliative",
	"Pall.Chemo.DC" = 														"Treat.Chemo.PalliativeDC",
#	"Adj.Chemo.Agent" = 													"Treat.Chemo.Adjuvant.Agent",
	"Received.GEM" = 														"Treat.Chemo.GEM",
	"X..New.Any.RTX.." = 													"Treat.Radio",
	"Location.non.resected.head.0..body.1" = 								"Path.LocationBody",
	"Tumour.Size" = 														"Path.Size",
	"Pre.bili.mg.dL" = 														"Path.Bilirubin.Preop",
	"Pre.Ca19.9" = 															"Path.Ca199.Preop",
	"Post.bili.mg.dL" = 													"Path.Bilirubin.Postop",
	"Corr.Post.CA19.9" = 													"Path.Ca199.Postop",
	"Cancer.type" = 														"Path.Subtype", 
	"Differentiation" = 													"Path.Differentiation",
	"LN.Involve" = 															"Path.LN.Involved",
	"LN.main.resect.total" = 												"Path.LN.Inspected",
	"X..Vascular.Invasion.Pos.1.." = 										"Path.Invasion.Vascular",
	"X..PN.Invasion.Pos.1.." = 												"Path.Invasion.Perineural",
	"pT.Stage" = 															"Stage.pT",
	"pN.Stage" = 															"Stage.pN",
	"pM.Stage" = 															"Stage.pM",
#	"pOverallStage" = 														"Stage.Overall",
#	"Clinical.Stage" = 														"Stage.Clinical",
	"BNIP3.Nuc" = 															"Molec.BNIP3.NucInt",
	"BNIP3.Cyto" = 															"Molec.BNIP3.CytoInt",
	"D1.Cyto.Lo" = 															"Molec.CCND1.CytoLo",
	"D1.Cyto.Hi" = 															"Molec.CCND1.CytoHi",
	"D1.Memb.Lo" = 															"Molec.CCND1.MembLo",
	"D1.Memb.Hi" = 															"Molec.CCND1.MembHi",
	"Grb7.Intensity" = 														"Molec.Grb7.Int",
	"Grb7.." = 																"Molec.Grb7.Percent",
	"HCNT3.HENT1" = 														"Molec.HCNT3PlusHENT1",
	"X.hENT1..." = 															"Molec.HENT1.Percent",
	"X.hENT1.Int." = 														"Molec.HENT1.Int",
	"HER2.ISH.AC" = 														"Molec.HER2",
	"HOXB2.." = 															"Molec.HOXB2.Percent",
	"HOXB2.intensity" = 													"Molec.HOXB2.Int",
	"RON.VS.AG" = 															"Molec.RON.Int",
	"X..S100A2.Intensity.." = 												"Molec.S100A2.Int",
	"X..S100A2...." = 														"Molec.S100A2.Percent",
	"S100A2.stroma.score" = 												"Molec.S100A2.StromaScore",
	"S100A4.Cyto.Int" = 													"Molec.S100A4.CytoInt",
	"S100A4.Cyto.." = 														"Molec.S100A4.CytoPercent",
	"X..S100A4.Nu.Int.C.." = 												"Molec.S100A4.NucInt",
	"X..S100A4.Nu...C.." = 													"Molec.S100A4.NucPercent")
#	"X.hENT1.H." = 															"Molec.HENT1.H",
#	"RON.H" = 																"Molec.RON.H",
#	"cyto.H" = 																"Molec.S100A4.CytoH",
#	"X..S100A4.Nu.H.." = 													"Molec.S100A4.NucH",
#	"M.0.2..1.1.2..2.0.1..3.0" = 											"",		# Stratified margin
#	"Pre.bili" = 															"",
#	"Normal.bilirubin" = 													"",
#	"Post.bili" = 															"",
#	"Change.bili" = 														"",
#	"Comorbid.1.1.DM..2.CP..3.PUD..4.obes..5.COPD..6.CVD.HT.lip.7.oth" = 	"",
#	"Current.smoker" = 														"",
#	"Former.smoker" = 														"",
#	"Heavy.ETOH" = 															"",

temp.sel = match(names(temp.fields), colnames(data))
data = data[,temp.sel]
origcols = origcols[temp.sel]
colnames(data) = temp.fields[colnames(data)]
names(origcols) = colnames(data)

calculateStage = function(T, N, M)
{
	# Based on Edge SB, Byrd DR, Compton CC, et al., eds.: AJCC Cancer Staging Manual. 7th ed. New York, NY: Springer, 2010, pp 241-9.
	stage = rep(NA, length(T))
	stage[                         M == "M1"] = "IV"
	stage[T == "T4"  & M == "M0"            ] = "III"
	stage[T != "T4"  & N == "N1" & M == "M0"] = "IIB"
	stage[T == "T3"  & N == "N0" & M == "M0"] = "IIA"
	stage[T == "T2"  & N == "N0" & M == "M0"] = "IB"
	stage[T == "T1"  & N == "N0" & M == "M0"] = "IA"
	stage[T == "Tis" & N == "N0" & M == "M0"] = "0"

	ordered(stage, levels = c("0", "IA", "IB", "IIA", "IIB", "III", "IV"))
}

data2 = data
data2$Patient.Sex = as.factor(substr(data2$Patient.Sex, 1, 1))
data2$Cohort.ICGC = as.numeric(data2$Cohort.ICGC) == 1 & !is.na(data2$Cohort.ICGC)
data2$History.PreviousMalignancy = as.numeric(data2$History.PreviousMalignancy) == 1
data2$History.FdrWithPancCancer = as.numeric(data2$History.FdrWithPancCancer) == 1
data2$History.FdrWithAnyCancer = as.numeric(data2$History.FdrWithAnyCancer) == 1
data2$History.Diagnosis.AgeAt = as.numeric(data2$History.Diagnosis.AgeAt)
data2$History.AlcoholLevel = ordered(as.numeric(data2$History.AlcoholLevel))
data2$History.Smoking.Status = factor(c("0" = "Never", "1" = "Ceased", "2" = "Current")[data2$History.Smoking.Status], levels = c("Never", "Ceased", "Current"))
data2$History.Smoking.PackYears = as.numeric(data2$History.Smoking.PackYears)
data2$History.Comorbid.Diabetes = as.numeric(data2$History.Comorbid.Diabetes) == 1
data2$History.Comorbid.ChronicPancreatitis = as.numeric(data2$History.Comorbid.ChronicPancreatitis) == 1
data2$History.Diagnosis.Date = as.Date("1900-1-1") + as.numeric(data2$History.Diagnosis.Date)
data2$History.Recurrence.Event = as.numeric(data2$History.Recurrence.Event) == 1
data2$History.Recurrence.Date = as.Date(gsub(" ", "", data2$History.Recurrence.Date), format = "%m/%d/%y")
data2$History.DSDeath.Event = as.numeric(data2$History.DSDeath.Event) == 1
data2$History.Death.Date = as.Date("1900-1-1") + as.numeric(data2$History.Death.Date)
data2$History.Death.EventTimeDays = as.numeric(data2$History.Death.EventTimeDays)
data2$History.Followup.Date = as.Date("1900-1-1") + as.numeric(data2$History.Followup.Date)
data2$Treat.Resected = data2$Treat.Resected == "0"												# TRUE => Resection, FALSE => Biopsy
data2$Treat.ProcedureWhipple = data2$Treat.ProcedureWhipple == "0"
data2$Treat.MarginPositive = data2$Treat.MarginPositive == "1"
data2$Treat.Chemo.Any = data2$Treat.Chemo.Any == "1"
data2$Treat.Chemo.Adjuvant = data2$Treat.Chemo.Adjuvant == "1"
data2$Treat.Chemo.Adjuvant.GE3Cycles = data2$Treat.Chemo.Adjuvant.GE3Cycles == "1"
data2$Treat.Chemo.Palliative = data2$Treat.Chemo.Palliative == "1"
data2$Treat.Chemo.PalliativeDC = data2$Treat.Chemo.PalliativeDC == "1"
data2$Treat.Chemo.GEM = data2$Treat.Chemo.GEM == "1"
data2$Treat.Radio = data2$Treat.Radio == "1"
data2$Path.LocationBody = data2$Path.LocationBody == "1"
data2$Path.Size = as.numeric(data2$Path.Size)
data2$Path.Bilirubin.Preop = as.numeric(data2$Path.Bilirubin.Preop)
data2$Path.Ca199.Preop = as.numeric(data2$Path.Ca199.Preop)
data2$Path.Bilirubin.Postop = as.numeric(data2$Path.Bilirubin.Postop)
data2$Path.Ca199.Postop = as.numeric(data2$Path.Ca199.Postop)
data2$Path.Subtype[is.na(data2$Path.Subtype)] = "NotSpecified"
data2$Path.Subtype = as.factor(data2$Path.Subtype)
data2$Path.Differentiation = ordered(as.numeric(substr(data2$Path.Differentiation, 1, 1)), levels = 1:4)
data2$Path.LN.Involved = as.numeric(data2$Path.LN.Involved)
data2$Path.LN.Inspected = as.numeric(data2$Path.LN.Inspected)
data2$Path.Invasion.Vascular = data2$Path.Invasion.Vascular == "1"
data2$Path.Invasion.Perineural = data2$Path.Invasion.Perineural == "1"
data2$Stage.pT[toupper(data2$Stage.pT) == "TX"] = NA
data2$Stage.pT = ordered(data2$Stage.pT, levels = c("Tis", "T1", "T2", "T3", "T4"))
data2$Stage.pN[toupper(data2$Stage.pN) == "NX"] = NA
data2$Stage.pN[toupper(data2$Stage.pN) == "N1A"] = "N1"
data2$Stage.pN[toupper(data2$Stage.pN) == "N1B"] = "N1"
data2$Stage.pN = ordered(data2$Stage.pN, levels = c("N0", "N1"))
data2$Stage.pM[toupper(data2$Stage.pM) == "MX"] = NA
data2$Stage.pM = ordered(data2$Stage.pM, levels = c("M0", "M1"))
data2$Stage.Overall = calculateStage(data2$Stage.pT, data2$Stage.pN, data2$Stage.pM)
data2$Molec.BNIP3.NucInt = ordered(as.numeric(data2$Molec.BNIP3.NucInt), levels = 0:3)
data2$Molec.BNIP3.CytoInt = ordered(as.numeric(data2$Molec.BNIP3.CytoInt), levels = 0:3)
data2$Molec.CCND1.CytoLo = ordered(as.numeric(data2$Molec.CCND1.CytoLo), levels = 0:3)
data2$Molec.CCND1.CytoHi = ordered(as.numeric(data2$Molec.CCND1.CytoHi), levels = 0:3)
data2$Molec.CCND1.MembLo = ordered(as.numeric(data2$Molec.CCND1.MembLo), levels = 0:3)
data2$Molec.CCND1.MembHi = ordered(as.numeric(data2$Molec.CCND1.MembHi), levels = 0:3)
data2$Molec.Grb7.Int = ordered(as.numeric(data2$Molec.Grb7.Int), levels = 0:3)
data2$Molec.Grb7.Percent = as.numeric(data2$Molec.Grb7.Percent)
data2$Molec.HCNT3PlusHENT1 = data2$Molec.HCNT3PlusHENT1 == "1"
data2$Molec.HENT1.Int = ordered(as.numeric(data2$Molec.HENT1.Int), levels = 0:3)
data2$Molec.HENT1.Percent = as.numeric(data2$Molec.HENT1.Percent)
data2$Molec.HER2 = data2$Molec.HER2 == "1"
data2$Molec.HOXB2.Int = ordered(as.numeric(data2$Molec.HOXB2.Int), levels = 0:3)
data2$Molec.HOXB2.Percent = as.numeric(data2$Molec.HOXB2.Percent)
data2$Molec.RON.Int = ordered(as.numeric(data2$Molec.RON.Int), levels = 0:3)
data2$Molec.S100A2.Int = ordered(as.numeric(data2$Molec.S100A2.Int), levels = 0:3)
data2$Molec.S100A2.Percent = as.numeric(data2$Molec.S100A2.Percent)
data2$Molec.S100A2.StromaScore = data2$Molec.S100A2.StromaScore == 1
data2$Molec.S100A4.CytoInt = ordered(as.numeric(data2$Molec.S100A4.CytoInt), levels = 0:3)
data2$Molec.S100A4.CytoPercent = as.numeric(data2$Molec.S100A4.CytoPercent)
data2$Molec.S100A4.NucInt = ordered(as.numeric(data2$Molec.S100A4.NucInt), levels = 0:3)
data2$Molec.S100A4.NucPercent = as.numeric(data2$Molec.S100A4.NucPercent)
data2$History.Status = gsub("[ .-]+", " ", data2$History.Status)
data2$History.Status = gsub("(^|[[:space:]])([[:alpha:]])", "\\1\\U\\2", data2$History.Status, perl = TRUE)
data2$History.Status = gsub("(^ +)|( +)$", "", data2$History.Status)
data2$History.Status[data2$History.Status == "Decaesed Of Disease"] = "Deceased Of Disease"
data2$History.Status[data2$History.Status %in% c("Deceased Of Other Disease", "Deceased Of Unknown Cause", "Deceased Other Cause", "Deceased Surgical Death")] = "Deceased Of Other Cause"
data2$History.Death.Event = grepl("Deceased", data2$History.Status)
data2 = data2[,colnames(data2) != "History.Status"]

# Remove duplicated patients.
temp = table(data2$Patient.ID)
temp = data2$Patient.ID %in% names(temp)[temp > 1]
data2 = data2[!temp,]

# Integrate in DC's data fixes
fixes = read.csv("./data/NSWPCN_DC_fixes_20150130.csv")
fixes = fixes[fixes$Patient.ID %in% data2$Patient.ID,]
fixes_match = match(fixes$Patient.ID, data2$Patient.ID)
for (i in setdiff(colnames(fixes), "Patient.ID"))
{
	if (grepl("Date", i))
		fixes[,i] = as.Date(fixes[,i], format = "%d/%m/%Y")
	if (i %in% colnames(data2))
		data2[fixes_match, i] = fixes[,i]
}
if ("Exclude" %in% colnames(fixes))
{
	data2$Exclude = FALSE
	data2$Exclude[fixes_match] = as.logical(fixes$Exclude)
	data2 = data2[!data2$Exclude,colnames(data2) != "Exclude"]
}

fixes = read.csv("./data/NSWPCN_DC_fixes_20150131.csv")
fixes = fixes[fixes$Patient.ID %in% data2$Patient.ID,]
fixes_match = match(fixes$Patient.ID, data2$Patient.ID)
for (i in setdiff(colnames(fixes), "Patient.ID"))
{
	if (grepl("Date", i))
		fixes[,i] = as.Date(fixes[,i], format = "%d/%m/%Y")
	if (i %in% colnames(data2))
		data2[fixes_match, i] = fixes[,i]
}
if ("Exclude" %in% colnames(fixes))
{
	data2$Exclude = FALSE
	data2$Exclude[fixes_match] = as.logical(fixes$Exclude)
	data2 = data2[!data2$Exclude,colnames(data2) != "Exclude"]
}

temp = NA
temp = ls()
rm(list = temp[grep("^temp", temp)])

data = data2
rm(data2)
rm(calculateStage)
rm(i, fixes, fixes_match)

#data = data[,order(colnames(data))]
origcols = origcols[match(colnames(data), names(origcols))]

save.image("01_NSWPCN.rda")

sessionInfo()

