options(stringsAsFactors = FALSE)

data = read.csv("./originals/Glasgow_S100A2_minimal_20140414.csv.xz", quote = "", header = TRUE)

data$Stage.Overall.Coarse = ordered(c("2" = "II", "3" = "III", "4" = "IV")[as.character(data$AJCC)], levels = c("I", "II", "III", "IV"))

temp.fields = list(
	"ID" = 																	"Patient.ID",
	"Gender" = 																"Patient.Sex",
	"Age" = 																"History.Diagnosis.AgeAt",
	"Operation" = 															"Treat.Procedure",
	"ClinTumourLocationHead" = 												"Path.Location",
	"TumourType" = 															"Path.Type", 
	"Differentiation" = 													"Path.Differentiation",
	"Grade" = 																"Path.Grade", 
	"Tstage" = 																"Stage.pT",
	"LymphNodestatus" = 													"Stage.pN",
#	"LymphaticInvasioncode" = ???
	"perineuralinvasioncode" = 												"Path.Invasion.Perineural",
	"vascularinvasioncode" = 												"Path.Invasion.Vascular",
	"Lymphnodescheckedtotal" = 												"Path.LN.Inspected", 
	"Lymphnoderatio" = 														"Path.LN.InvolvedFraction", 
	"resectionmargincode.positivelessthan1mm." = 							"Treat.MarginPositive",
	"Veinresection" = 														"Treat.VeinResection", 
	"TmsizeMax" = 															"Path.Size",
	"Survivalmonths" = 														"History.Death.EventTimeDays", 
	"Alive0deadcancer1deadother2" = 										"History.Death.Cause", 
	"AdjChemotherapy" = 													"Treat.Chemo.Adjuvant",
	"neoadjuvantchemo" = 													"Treat.Chemo.Neoadjuvant",
	"S100A2combinedTMAscode" = 												"Molec.S100A2.DCThresh", 
	"S100A4anyscoreallTMA" = 												"Molec.S100A4.DCThresh")

data = data[,match(names(temp.fields), colnames(data))]
colnames(data) = temp.fields[colnames(data)]

data$Patient.ID = paste("Glasgow_", data$Patient.ID, sep = "")
rownames(data) = data$Patient.ID
data$Patient.Sex = as.factor(data$Patient.Sex)
data$Treat.ProcedureWhipple = data$Treat.Procedure == "Whipple"
data$Path.LocationBody = data$Path.Location != "HOP"
data$History.Death.EventTimeDays = data$History.Death.EventTimeDays * 365.25 / 12
data$History.DSDeath.Event = data$History.Death.Cause == "1"
data$History.ACDeath.Event = data$History.Death.Cause == "1" | data$History.Death.Cause == "2"
data$Path.Type = as.factor(data$Path.Type)
data$Path.Differentiation = c("MOD" = "Moderate", "MOD.POOR" = "Moderate-Poor", "MOD.WELL" = "Well-Moderate", 
	"MODERATE" = "Moderate", "MODERATE.POOR" = "Moderate-Poor", "MODERATE.WELL" = "Well-Moderate",
	"MODERATELY" = "Moderate", "POOR" = "Poor", "POORLY" = "Poor", "WELL" = "Well", "WELL.MOD" = "Well-Moderate",
	"WELL.MODERATE" = "Well-Moderate", "WELL.POOR" = "Well-Poor")[gsub(" +", ".", gsub("-", " ", gsub("/", " ", toupper(data$Path.Differentiation))))]
data$Path.Differentiation = c("Well" = 1, "Well-Moderate" = 2, "Moderate" = 2, "Moderate-Poor" = 3, "Poor" = 3, "Well-Poor" = 2)[data$Path.Differentiation]
data$Path.Differentiation = ordered(data$Path.Differentiation, levels = 1:4)
data$Path.Grade = ordered(data$Path.Grade, levels = c("Low", "High"))
data$Stage.pT = ordered(c("T1", "T2", "T3", "T4")[data$Stage.pT], levels = c("Tis", "T1", "T2", "T3", "T4"))
data$Stage.pN = ordered(c("N0", "N1")[data$Stage.pN + 1], levels = c("N0", "N1"))
data$Path.Invasion.Perineural = data$Path.Invasion.Perineural == 1
data$Path.Invasion.Vascular = data$Path.Invasion.Vascular == 1
data$Path.LN.InvolvedFraction = pmin(1, data$Path.LN.InvolvedFraction)
data$Path.LN.Involved = round(as.numeric(data$Path.LN.Inspected) * as.numeric(data$Path.LN.InvolvedFraction))
data$Treat.MarginPositive = data$Treat.MarginPositive == 1
data$Treat.VeinResection = data$Treat.VeinResection == 1
data$Treat.Chemo.Adjuvant = data$Treat.Chemo.Adjuvant == 1
data$Treat.Chemo.Neoadjuvant = data$Treat.Chemo.Neoadjuvant == 1
data$Molec.S100A2.DCThresh = data$Molec.S100A2.DCThresh == 1
data$Molec.S100A4.DCThresh = data$Molec.S100A4.DCThresh == 1
data$History.Death.Cause = as.factor(data$History.Death.Cause)
data$History.Diagnosis.AgeAt.Cent = data$History.Diagnosis.AgeAt - 68
data$Path.Size.Cent = data$Path.Size - 30

rm(temp.fields)
saveRDS(data, "05_Glasgow.rds")

sessionInfo()

