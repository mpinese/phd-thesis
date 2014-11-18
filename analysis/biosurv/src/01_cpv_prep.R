#!~/bin/Rscript

# 01_cpv_prep.R -- load CPVs from Clare Watson's database dump, and 
# perform basic recoding.

# Read the CPVs
options(stringsAsFactors = FALSE, echo = TRUE)
cpvs = read.csv("../data/originals/Data request_MPinese 1.2.csv")

levelToNA = function(x, nalevels)
{
	x = as.character(x)
	x[x %in% nalevels] = NA
	factor(x)
}

# Recode simple cases
cpvs$PatientID = paste("APGI_", cpvs$PatientID, sep = "")
cpvs$Gender = factor(cpvs$Gender)
cpvs$Ethnicity = levelToNA(cpvs$Ethnicity, "")
cpvs$Country = factor(cpvs$Country)
cpvs$Smoker = levelToNA(cpvs$Smoker, "")
cpvs$Operation.Date = as.Date(cpvs$Operation.Date, format = "%d/%m/%Y")
cpvs$Type.of.Procedure = factor(cpvs$Type.of.Procedure)
cpvs$Date.of.Diagnosis = as.Date(cpvs$Date.of.Diagnosis, format = "%d/%m/%Y")
cpvs$Histologic.Type = factor(cpvs$Histologic.Type)
cpvs$Histologic.Subtype = levelToNA(cpvs$Histologic.Subtype, "")
cpvs$Grade = ordered(as.numeric(substr(cpvs$Grade, 1, 1)), levels = 1:4)
cpvs$Tumour.Location = factor(cpvs$Tumour.Location)
cpvs$Tumor.Size..mm. = as.numeric(cpvs$Tumor.Size..mm.)
cpvs$Perineural.Inv = levelToNA(cpvs$Perineural.Inv, c("Not reported", ""))
cpvs$Vascular.Inv = levelToNA(cpvs$Vascular.Inv, c("Not reported", ""))
cpvs$Regional.LN.Total = as.numeric(cpvs$Regional.LN.Total)
cpvs$Regional.LN.Inv = as.numeric(cpvs$Regional.LN.Inv)
cpvs$SepRec.LN.Total = as.numeric(cpvs$Regional.LN.Total)
cpvs$SepRec.LN.Inv = as.numeric(cpvs$Regional.LN.Inv)
cpvs$Excision.Status = levelToNA(substr(cpvs$Excision.Status, 1, 2), "RX")
cpvs$Pancreatic.Margin = levelToNA(cpvs$Pancreatic.Margin, c("Not Applicable", "Not Reported"))
cpvs$Periunc.Margin = levelToNA(cpvs$Periunc.Margin, c("Not Applicable", "", "Not Reported"))
cpvs$PV.Groove.Margin = levelToNA(cpvs$PV.Groove.Margin, c("Not Applicable", "", "Not Reported"))
cpvs$Retrop.Margin = levelToNA(cpvs$Retrop.Margin, c("Not Applicable", "", "Not Reported"))
cpvs$CBD.Margin = levelToNA(cpvs$CBD.Margin, c("Not Applicable", "", "Not Reported"))
cpvs$Duodenal.Margin = levelToNA(cpvs$Duodenal.Margin, c("Not Applicable", "", "Not Reported"))
cpvs$Gastric.Margin = levelToNA(cpvs$Gastric.Margin, c("Not Applicable", "", "Not Reported"))
cpvs$Stage.Version = levelToNA(cpvs$Stage.Version, "")
cpvs$pM = ordered(levelToNA(cpvs$pM, "MX"), levels = c("M0", "M1"))
cpvs$pN = ordered(levelToNA(substr(cpvs$pN, 1, 2), "NX"), levels = c("N0", "N1"))
cpvs$pT = ordered(levelToNA(cpvs$pT, "TX"), levels = c("Tis", "T1", "T2", "T3", "T4"))
cpvs$Status = factor(cpvs$Status)
cpvs$Date.of.Death = as.Date(cpvs$Date.of.Death, format = "%d/%m/%Y")
cpvs$Cause.of.Death = levelToNA(cpvs$Cause.of.Death, c("", "Unknown"))
cpvs$Date.of.Recurrence = as.Date(cpvs$Date.of.Recurrence, format = "%d/%m/%Y")
cpvs$DOLFU = as.Date(cpvs$DOLFU, format = "%d-%b-%y")

# Perform more complex recodings
cpvs$No..Years.Smoked.Approx = c(
	"< 10" =   5, 
	"11-20" = 15, 
	"21-30" = 25, 
	"31-40" = 35, 
	"41-50" = 45,
	"> 50" =  55)[gsub("Smoked ", "", gsub(" years", "", cpvs$No..Years.Smoked))]
cpvs$No..Years.Smoked.Approx[cpvs$No..Years.Smoked == ""] = NA

cpvs$Overall.Stage[cpvs$Overall.Stage %in% c("0", "NA (Biopsy)")] = NA
cpvs$Overall.Stage = gsub(" .*", "", cpvs$Overall.Stage)
cpvs$Overall.Stage = ordered(cpvs$Overall.Stage, levels = c("IA", "IB", "IIA", "IIB", "III", "IV"))

cpvs$Recurrence[cpvs$Recurrence == ""] = "Not observed"
cpvs$Recurrence[cpvs$Recurrence %in% c("Not applicable", "Unknown")] = NA
cpvs$Recurrence = ordered(cpvs$Recurrence, levels = c("Not observed", "Suspected", "Confirmed"))

temp.sites = sort(unique(unlist(strsplit(cpvs$Recurrence.Site, ", "))), decreasing = TRUE)
for (temp.site in temp.sites)
{
	cpvs = cbind(cpvs, sapply(strsplit(cpvs$Recurrence.Site, ", "), function(x) temp.site %in% x))
	if (temp.site == "Other (please specify)") { temp.site = "Other" }
	colnames(cpvs)[ncol(cpvs)] = make.names(sprintf("Recurrence.Site.%s", gsub(" ", "", temp.site)))
}


# Calculate derived variables
cpvs$Pack.Years = cpvs$No..Years.Smoked.Approx * cpvs$No..of.Smokes.Day / 20


# Rename variabes for some consistency
temp.renamer = c(
"PatientID"							= "Patient.ID",
"Gender"							= "Patient.Gender",
"Ethnicity"							= "Patient.Ethnicity",
"Country"							= "Patient.Country",
"DOLFU"								= "History.LastFollowup.Date",
"Pack.Years"						= "History.Smoking.PackYears",
"Date.of.Diagnosis"					= "History.Diagnosis.Date",
"Age.at.Diagnosis..Years."			= "History.Diagnosis.AgeAtYears",
"Operation.Date"					= "History.Surgery.Date",
"Type.of.Procedure"					= "Treat.Surgery.Procedure",
"Excision.Status"					= "Treat.Surgery.ExcisionStatus",
"Pancreatic.Margin"					= "Treat.Surgery.Margin.Pancreatic",
"Pancreatic.Margin..mm."			= "Treat.Surgery.MarginSizeMm.Pancreatic",
"Periunc.Margin"					= "Treat.Surgery.Margin.Periunc",
"Periunc.Margin..mm."				= "Treat.Surgery.MarginSizeMm.Periunc", 
"PV.Groove.Margin"					= "Treat.Surgery.Margin.PVGroove",
"PV.Groove.Margin..mm."				= "Treat.Surgery.MarginSizeMm.PVGroove", 
"Retrop.Margin"						= "Treat.Surgery.Margin.Retrop",
"Retrop.Margin..mm."				= "Treat.Surgery.MarginSizeMm.Retrop", 
"CBD.Margin"						= "Treat.Surgery.Margin.CBD",
"CBD.Margin..mm."					= "Treat.Surgery.MarginSizeMm.CBD", 
"Duodenal.Margin"					= "Treat.Surgery.Margin.Duodenal",
"Duodenal.Margin..mm."				= "Treat.Surgery.MarginSizeMm.Duodenal", 
"Gastric.Margin"					= "Treat.Surgery.Margin.Gastric",
"Gastric.Margin..mm."				= "Treat.Surgery.MarginSizeMm.Gastric", 
"Margin.Comments"					= "Treat.Surgery.Margin.Comments",
"Histologic.Type"					= "Path.HistoType",
"Histologic.Subtype"				= "Path.HistoType.Subtype",
"Grade"								= "Path.Grade",
"Tumour.Location"					= "Path.TumourLocation",
"Tumor.Size..mm."					= "Path.TumourSizeMm",
"Perineural.Inv"					= "Path.Invasion.PN",
"Vascular.Inv"						= "Path.Invasion.VS",
"Regional.LN.Total"					= "Path.Nodes.Regional.Total",
"Regional.LN.Inv"					= "Path.Nodes.Regional.Involved",
"SepRec.LN.Total"					= "Path.Nodes.SepRec.Total",
"SepRec.LN.Inv"						= "Path.Nodes.SepRec.Involved",
"Stage.Version"						= "Staging.Version",
"pM"								= "Staging.pM",
"pN"								= "Staging.pN",
"pT"								= "Staging.pT",
"Overall.Stage"						= "Staging.Stage",
"Recurrence"						= "History.Recurrence",
"Date.of.Recurrence"				= "History.Recurrence.Date",
"Recurrence.Site.Stomach"			= "History.Recurrence.Site.Stomach",
"Recurrence.Site.Peritoneum"		= "History.Recurrence.Site.Peritoneum",
"Recurrence.Site.PancreaticRemnant" = "History.Recurrence.Site.PancRemnant",
"Recurrence.Site.PancreaticBed"		= "History.Recurrence.Site.PancBed",
"Recurrence.Site.Other"				= "History.Recurrence.Site.Other",
"Recurrence.Site.Omentum"			= "History.Recurrence.Site.Omentum",
"Recurrence.Site.Mesentery"			= "History.Recurrence.Site.Mesentery",
"Recurrence.Site.LymphNodes"		= "History.Recurrence.Site.LymphNodes",
"Recurrence.Site.Lung"				= "History.Recurrence.Site.Lung",
"Recurrence.Site.Liver"				= "History.Recurrence.Site.Liver",
"Recurrence.Site.Brain"				= "History.Recurrence.Site.Brain",
"Recurrence.Site.Bone"				= "History.Recurrence.Site.Bone",
"Status"							= "History.Status",
"Date.of.Death"						= "History.Death.Date",
"Cause.of.Death"					= "History.Death.Cause",
"Smoker"							= NA,
"No..Years.Smoked"					= NA,
"No..of.Smokes.Day"					= NA,
"Recurrence.Site"					= NA,
"Recurrence.Site..Other."			= NA,
"Length.FU..Days."					= NA,
"No..Years.Smoked.Approx"			= NA)

cpvs = cpvs[,names(temp.renamer)]
cpvs = cpvs[,colnames(cpvs) %in% names(temp.renamer)[!is.na(temp.renamer)]]
colnames(cpvs) = temp.renamer[colnames(cpvs)]


# Further derived variables for survival modelling
cpvs$Surv.Event.Death = grepl("Deceased", cpvs$History.Status)*1
cpvs$Surv.EventTimeFromDiag.Death = ifelse(cpvs$Surv.Event.Death == 1, cpvs$History.Death.Date - cpvs$History.Diagnosis.Date, cpvs$History.LastFollowup.Date - cpvs$History.Diagnosis.Date)
cpvs$Surv.EventTimeFromSurg.Death = ifelse(cpvs$Surv.Event.Death == 1, cpvs$History.Death.Date - cpvs$History.Surgery.Date, cpvs$History.LastFollowup.Date - cpvs$History.Surgery.Date)
cpvs$Surv.EventTimeFromRec.Death = ifelse(cpvs$Surv.Event.Death == 1, cpvs$History.Death.Date - cpvs$History.Recurrence.Date, cpvs$History.LastFollowup.Date - cpvs$History.Recurrence.Date)
cpvs$Surv.Event.DSDeath = (cpvs$History.Status == "Deceased - Of Disease")*1
cpvs$Surv.EventTimeFromDiag.DSDeath = ifelse(cpvs$Surv.Event.DSDeath == 1, cpvs$History.Death.Date - cpvs$History.Diagnosis.Date, cpvs$History.LastFollowup.Date - cpvs$History.Diagnosis.Date)
cpvs$Surv.EventTimeFromSurg.DSDeath = ifelse(cpvs$Surv.Event.DSDeath == 1, cpvs$History.Death.Date - cpvs$History.Surgery.Date, cpvs$History.LastFollowup.Date - cpvs$History.Surgery.Date)
cpvs$Surv.EventTimeFromRec.DSDeath = ifelse(cpvs$Surv.Event.DSDeath == 1, cpvs$History.Death.Date - cpvs$History.Recurrence.Date, cpvs$History.LastFollowup.Date - cpvs$History.Recurrence.Date)
cpvs$Surv.Event.Recurrence = (cpvs$History.Recurrence == "Confirmed")*1
cpvs$Surv.EventTimeFromDiag.Recurrence = ifelse(cpvs$Surv.Event.Recurrence == 1, cpvs$History.Recurrence.Date - cpvs$History.Diagnosis.Date, cpvs$History.LastFollowup.Date - cpvs$History.Diagnosis.Date)
cpvs$Surv.EventTimeFromSurg.Recurrence = ifelse(cpvs$Surv.Event.Recurrence == 1, cpvs$History.Recurrence.Date - cpvs$History.Surgery.Date, cpvs$History.LastFollowup.Date - cpvs$History.Surgery.Date)


rm(levelToNA, temp.renamer, temp.site, temp.sites)

saveRDS(cpvs, "../data/01_cpvs.rds")

sessionInfo()
