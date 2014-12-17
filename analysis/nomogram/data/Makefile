.PHONY: all
.SECONDARY: 

RSCRIPT = ~/bin/Rscript

all: 01_NSWPCN.rda 03_NSWPCN_subset.rda 04_NSWPCN_fits.rda 05_Glasgow.rds

01_NSWPCN.rda: originals/Total_cohort_MP_05022013.txt.xz 01_load_NSWPCN.R
	$(RSCRIPT) 01_load_NSWPCN.R | tee 01_load_NSWPCN.Rout

03_NSWPCN_subset.rda: 01_NSWPCN.rda 03_subset_fix_NSWPCN.R
	$(RSCRIPT) 03_subset_fix_NSWPCN.R | tee 03_subset_fix_NSWPCN.Rout

04_NSWPCN_fits.rda: 03_NSWPCN_subset.rda 04_train_NSWPCN.R
	$(RSCRIPT) 04_train_NSWPCN.R | tee 04_train_NSWPCN.Rout

05_Glasgow.rds: originals/Glasgow_S100A2_minimal_20140414.csv.xz 05_load_Glasgow.R
	$(RSCRIPT) 05_load_Glasgow.R | tee 05_load_Glasgow.Rout
