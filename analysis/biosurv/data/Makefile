.SECONDARY: 

.PHONY: all

all: 07_data_for_SIS.rda 15_validation.rda

07_data_for_SIS.rda: 05_gexnormed_cpv_matched.rda 07_SIS_prep.R
	Rscript 07_SIS_prep.R > runlogs/07_SIS_prep.Rout

05_gexnormed_cpv_matched.rda: 04_gexnormed_cpv.rda 05_match_gex_cpv.R
	Rscript 05_match_gex_cpv.R > runlogs/05_match_gex_cpv.Rout

04_gexnormed_cpv.rda: 03_gex_cpv.rda 04_norm_gex.R
	Rscript 04_norm_gex.R > runlogs/04_norm_gex.Rout

gex_cpv: 01_cpvs.rds 02_gex.rda 03_subset_gex_cpv.R
	Rscript 03_subset_gex_cpv.R > runlogs/03_subset_gex_cpv.Rout

02_gex.rda: 02_gex_process.R
	Rscript 02_gex_process.R > runlogs/02_gex_process.Rout

01_cpvs.rds: originals/Data_request_MPinese_1.2.csv 01_cpv_prep.R
	Rscript 01_cpv_prep.R > runlogs/01_cpv_prep.Rout

15_validation.rda: 15_prep_validation.R
	Rscript 15_prep_validation.R > runlogs/15_prep_validation.Rout
