.PHONY: all clean ultraclean
.SECONDARY: 

#RSCRIPT = ~/bin/Rscript
RSCRIPT = Rscript

all: 01_NSWPCN.rda 03_NSWPCN_subset.rda 05_train_NSWPCN_2.pdf 05_final_model.rds 06_Glasgow.rds 06_Dresden.rds 06_Dresden2.rds 06_APGI.rds 07_validate.pdf

clean:
	rm -f 01_NSWPCN.rda 03_NSWPCN_subset.rda 05_final_model.rds 06_Glasgow.rds 06_Dresden.rds 06_APGI.rds
	rm -f *.Rout
	rm -f *.log *.aux *.pdf *.rda glmulti.* *.tex

ultraclean: clean
	rm -rf figure/ cache/ tikzMetrics

01_NSWPCN.rda: data/Total_cohort_MP_05022013.txt.xz 01_load_NSWPCN.R
	$(RSCRIPT) 01_load_NSWPCN.R | tee 01_load_NSWPCN.Rout

03_NSWPCN_subset.rda: 01_NSWPCN.rda 03_subset_fix_NSWPCN.R
	$(RSCRIPT) 03_subset_fix_NSWPCN.R | tee 03_subset_fix_NSWPCN.Rout

05_final_model.rds 05_train_NSWPCN_2.pdf 05_train_NSWPCN_2.rda: 03_NSWPCN_subset.rda 05_train_NSWPCN_2.Rnw
	$(RSCRIPT) -e "knitr::knit2pdf('05_train_NSWPCN_2.Rnw')"

06_Glasgow.rds: data/Glasgow_S100A2_minimal_20140414.csv 06_load_Glasgow.R
	$(RSCRIPT) 06_load_Glasgow.R | tee 06_load_Glasgow.Rout

06_Dresden.rds 06_Dresden2.rds data/06_Dresden.csv: data/clin_data_rev_David15012015_DC.csv data/Dresden_S100A2_S100A4_Scores.csv 06_load_Dresden.R
	$(RSCRIPT) 06_load_Dresden.R | tee 06_load_Dresden.Rout

06_APGI.rds: data/APGI_20150214.csv ../biosurv/data/01_cpvs.rds 06_load_APGI.R
	$(RSCRIPT) 06_load_APGI.R | tee 06_load_APGI.Rout

07_validate.pdf: 05_final_model.rds 06_Glasgow.rds data/APGI_20150214.csv 06_Dresden2.rds 07_validate.Rnw
	$(RSCRIPT) -e "knitr::knit2pdf('07_validate.Rnw')"
