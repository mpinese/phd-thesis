.SECONDARY: 

.PHONY: all

all: 06_EDA.done 14_SIS_NMF_CV_results.rds

06_EDA.done: ../data/05_gexnormed_cpv_matched.rda 06_EDA.R
	Rscript 06_EDA.R > runlogs/06_EDA.Rout && touch 06_EDA.done

14_SIS_NMF_CV_results.rds: ../data/07_data_for_SIS.rda 14_SIS_NMF_CV.R
	nice Rscript 14_SIS_NMF_CV.R > runlogs/14_SIS_NMF_CV.Rout
