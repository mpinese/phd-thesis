#!~/bin/Rscript

# 02_gex_process.R -- load 

options(stringsAsFactors = FALSE, echo = TRUE)

library(lumidat)

# Download all array files to a scratch location
temp.arrayids = system("smbclient -A ~/garvan.creds '\\\\gagri.garvan.unsw.edu.au\\GRIW' -D ICGCPancreas/icgc_data_mirror/array_raw_data/GEX_array -c 'ls' | grep ' D[A ]' | grep -E '[0-9]{10}' | sed -E 's/ +/ /g' | cut -f 2 -d ' '",
	ignore.stderr = TRUE, intern = TRUE)
temp.oldwd = getwd()
dir.create("/share/Temp/marpin/APGILumiGex", recursive = TRUE)
for (temp.arrayid in temp.arrayids)
{
	temp.dir = sprintf("/share/Temp/marpin/APGILumiGex/%s", temp.arrayid)
	dir.create(temp.dir)
	setwd(temp.dir)
	system(sprintf("smbclient -A ~/garvan.creds '\\\\gagri.garvan.unsw.edu.au\\GRIW' -D ICGCPancreas/icgc_data_mirror/array_raw_data/GEX_array/%s -c 'prompt;mget *.idat'", temp.arrayid))
}
setwd(temp.oldwd)


# Get a list of all array files we have
input_files = list.files(path = '/share/Temp/marpin/APGILumiGex',
	pattern = '^[0-9]{10}_[A-L]_Grn\\.idat$',
	full.names = TRUE,
	recursive = TRUE,
	include.dirs = TRUE)
length(input_files)

# Get the project annotation data
sample_data = as.data.frame(read.table("originals/garvan_sample_overview.tsv.xz", sep = "\t", quote = "", comment = "", header = TRUE))
colnames(sample_data)[colnames(sample_data) == "submitted_patient_id"] = "patient_id"

# Get purity estimates.  This is complicated by the fact that the purity estimates are
# calculated from DNA extractions, whereas the microarray data are from RNA extractions.
# Although the extractions are performed in tandem from the same piece of tissue, it's
# impossible with the data to hand to definitively conclude that a given DNA extraction
# (with purity estimate), and a given RNA extraction (with expression data), originated
# from the same tissue fragment.  For now address this problem in an approximate manner:
# if only one DNA and one RNA entry is present for a patient tumour, assume that they
# are from the same extraction (and therefore tissue fragment), and associate the DNA 
# entry's purity with the RNA.  Otherwise (in the case of many-to-one, one-to-many, or 
# many-to-many associations), let the purity be unknown.
temp.purity = sample_data[!is.na(sample_data$sample_qpure_score) & sample_data$sample_code == "7:Primary tumour",c("biospecimen_id", "patient_id", "sample_qpure_score")]
# Remove patients with more than one purity entry:
temp.purity = temp.purity[!(temp.purity$patient_id %in% names(which(table(temp.purity$patient_id) > 1))),]
temp.purity$purity_qpure = as.numeric(as.vector(temp.purity$sample_qpure_score)) / 100

sample_data = sample_data[sample_data$arrays_completed != "" & sample_data$material == "2:RNA", c("patient_id", "sample_code", "material", "arrays_completed")]

# Remove patients with more than one microarray:
temp.purity = temp.purity[!(temp.purity$patient_id %in% names(which(table(sample_data[sample_data$sample_code == "7:Primary tumour",]$patient_id) > 1))),]

# Add purities to remaining samples
sample_data$purity_qpure = temp.purity$purity_qpure[match(sample_data$patient_id, temp.purity$patient_id)]
# Purity only applies to tumour samples
sample_data$purity_qpure[sample_data$sample_code != "7:Primary tumour"] = NA
#sample_data$purity_qpure = as.numeric(as.vector(sample_data$purity_qpure))

# Parse the arrays_completed field, removing duplicate entries
temp = sapply(strsplit(sample_data$arrays_completed, ", "), unique)
temp = sapply(temp, function(x) sapply(x, function(y) strsplit(y, " ")[[1]][4]))

# Reconstruct the sample_data variable, with multiple entries for replicated samples.
sample_data_2 = c()
for (i in 1:nrow(sample_data))
{
	for (j in 1:length(temp[[i]]))
	{
		sample_data_2 = rbind(sample_data_2, unlist(c(sample_data[i,], temp[[i]][[j]])))
	}
}

sample_data = sample_data_2
rm(sample_data_2, i, j)
sample_data = sample_data[,colnames(sample_data) != "arrays_completed"]
colnames(sample_data)[ncol(sample_data)] = "exp_array_id"
sample_data = data.frame(sample_data)
sample_data$exp_slide_id = substr(sample_data$exp_array_id, 1, 10)
for (i in 1:ncol(sample_data))	{ if (colnames(sample_data)[i] != "purity_qpure") { sample_data[,i] = as.factor(sample_data[,i]) } }
sample_data = sample_data[,colnames(sample_data) != "material"]
rm(i)
rownames(sample_data) = sample_data$exp_array_id

temp = table(sample_data$patient_id, sample_data$sample_code)
temp[apply(temp, 1, max) > 1, apply(temp, 2, max) > 1]

input_files_id = gsub("_Grn\\.idat", "", gsub(".*/", "", input_files))
names(input_files) = input_files_id
rm(input_files_id)

complete_ids = intersect(names(input_files), rownames(sample_data))
input_files = input_files[complete_ids]
sample_data = sample_data[complete_ids,]
sample_data = cbind(sample_data, path = input_files)
rm(input_files, complete_ids)
nrow(sample_data)

system.time(probe_raw <- lumiR.idat(
	files = sample_data$path, 
	path = "",
	manifestfile = "HumanHT-12_V4_0_R2_15002873_B.txt",
	probeID = "ProbeID",
	verbose = TRUE,
	memory = "-Xmx8192m"))

save(sample_data, probe_raw, file = "02_gex.rda", compress = "gzip", compression_level = 4)

sessionInfo()
