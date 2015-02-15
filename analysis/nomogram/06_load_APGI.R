scores.apgi = read.csv("./data/APGI_20150214.csv")
data.apgi = readRDS("../biosurv/data/01_cpvs.rds")
data.apgi$A2 = scores.apgi$A2[match(data.apgi$Patient.ID, scores.apgi$PatientID)]
data.apgi$A4 = scores.apgi$A4[match(data.apgi$Patient.ID, scores.apgi$PatientID)]
rm(scores.apgi)

data.apgi$Path.LN.Inspected = data.apgi$Path.Nodes.Regional.Total
data.apgi$Path.LN.Involved = data.apgi$Path.Nodes.Regional.Involved
data.apgi$Path.LN.Negative = data.apgi$Path.LN.Inspected - data.apgi$Path.LN.Involved
data.apgi$History.Diagnosis.AgeAt = data.apgi$History.Diagnosis.AgeAtYears
data.apgi$History.Diagnosis.AgeAt.Cent = data.apgi$History.Diagnosis.AgeAt - 68
data.apgi$Path.Size = data.apgi$Path.TumourSizeMm
data.apgi$Path.Size.Cent = data.apgi$Path.Size - 30
data.apgi$Patient.Sex = data.apgi$Patient.Gender
data.apgi$SexM = data.apgi$Patient.Sex == "M"
data.apgi$Treat.MarginPositive = data.apgi$Treat.Surgery.ExcisionStatus != "R0"
data.apgi$AgeCent = data.apgi$History.Diagnosis.AgeAt.Cent
data.apgi$SizeCent = data.apgi$Path.Size.Cent
data.apgi$A2 = data.apgi$A2 == 1
data.apgi$A4 = data.apgi$A4 == 1
data.apgi$Stage.pT = data.apgi$Staging.pT
data.apgi$Stage.pT.Simplified = c("T1" = "T1", "T2" = "T2", "T3" = "T34", "T4" = "T34")[as.character(data.apgi$Stage.pT)]
data.apgi$Path.LocationBody = !grepl("head", data.apgi$Path.TumourLocation, ignore.case = TRUE)
data.apgi$Path.LocationBody[data.apgi$Path.TumourLocation == ""] = NA
data.apgi$Path.Differentiation = data.apgi$Path.Grade
data.apgi$LocBody = data.apgi$Path.LocationBody
data.apgi$Time = data.apgi$Surv.EventTimeFromSurg.DSDeath
data.apgi$DSD = data.apgi$Surv.Event.DSDeath

temp.sel = apply(!is.na(data.apgi[,c("Path.LN.Inspected", "Path.LN.Involved", "Path.LN.Negative", "SexM", "AgeCent", "SizeCent", "A2", "A4", "LocBody", "Time", "DSD")]), 1, all) & data.apgi$Path.HistoType == "Pancreatic Ductal Adenocarcinoma"
data.apgi = data.apgi[temp.sel,]

saveRDS(data.apgi, file = "06_APGI.rds")
