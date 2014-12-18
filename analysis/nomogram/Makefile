.PHONY: all
.SECONDARY: 

#RSCRIPT = ~/bin/Rscript
RSCRIPT = Rscript

all: 01_NSWPCN.rda 03_NSWPCN_subset.rda 04_NSWPCN_fits.rda 05_Glasgow.rds

01_NSWPCN.rda: data/Total_cohort_MP_05022013.txt.xz 01_load_NSWPCN.R
	$(RSCRIPT) 01_load_NSWPCN.R | tee 01_load_NSWPCN.Rout

03_NSWPCN_subset.rda: 01_NSWPCN.rda 03_subset_fix_NSWPCN.R
	$(RSCRIPT) 03_subset_fix_NSWPCN.R | tee 03_subset_fix_NSWPCN.Rout

04_NSWPCN_fits.rda: 04_molec_preop_bic.rds 04_molec_preop_aicc.rds 04_molec_postop_bic.rds 04_molec_postop_aicc.rds 04_conv_preop_bic.rds 04_conv_preop_aicc.rds 04_conv_postop_bic.rds 04_conv_postop_aicc.rds 04_rsf_molec_preop.rds 04_train_NSWPCN.R
	$(RSCRIPT) 04_train_NSWPCN.R 9 | tee 04_train_NSWPCN.Rout

04_molec_preop_bic.rds: 03_NSWPCN_subset.rda 04_train_NSWPCN.R
	$(RSCRIPT) 04_train_NSWPCN.R 1

04_molec_preop_aicc.rds: 03_NSWPCN_subset.rda 04_train_NSWPCN.R
	$(RSCRIPT) 04_train_NSWPCN.R 2

04_conv_preop_bic.rds: 03_NSWPCN_subset.rda 04_train_NSWPCN.R
	$(RSCRIPT) 04_train_NSWPCN.R 3

04_conv_preop_aicc.rds: 03_NSWPCN_subset.rda 04_train_NSWPCN.R
	$(RSCRIPT) 04_train_NSWPCN.R 4

04_molec_postop_bic.rds: 03_NSWPCN_subset.rda 04_train_NSWPCN.R
	$(RSCRIPT) 04_train_NSWPCN.R 5

04_molec_postop_aicc.rds: 03_NSWPCN_subset.rda 04_train_NSWPCN.R
	$(RSCRIPT) 04_train_NSWPCN.R 6

04_conv_postop_bic.rds: 03_NSWPCN_subset.rda 04_train_NSWPCN.R
	$(RSCRIPT) 04_train_NSWPCN.R 7

04_conv_postop_aicc.rds: 03_NSWPCN_subset.rda 04_train_NSWPCN.R
	$(RSCRIPT) 04_train_NSWPCN.R 8

04_rsf_molec_preop.rds: 03_NSWPCN_subset.rda 04_train_NSWPCN.R
	$(RSCRIPT) 04_train_NSWPCN.R 10

05_Glasgow.rds: data/Glasgow_S100A2_minimal_20140414.csv.xz 05_load_Glasgow.R
	$(RSCRIPT) 05_load_Glasgow.R | tee 05_load_Glasgow.Rout
