.PHONY: all clean ultraclean figclean cacheclean

RSCRIPT = Rscript

all: 02_surv_exp1.rda 02_surv_exp1.pdf 03_surv_exp2.rda 03_surv_exp2.pdf

clean: 
	rm -f 01_surv_apgi_calcs.Rout 01_surv_apgi_calcs.rda 02_surv_exp1.rda 02_surv_exp1.pdf tikzMetrics .RData

figclean:
	rm -rf ./figure/

cacheclean:
	rm -rf ./cache/

ultraclean: clean figclean cacheclean

01_surv_apgi_calcs.rda 01_surv_apgi_calcs.Rout: 01_surv_apgi_calcs.R ../biosurv/data/07_data_for_SIS.rda
	$(RSCRIPT) 01_surv_apgi_calcs.R | tee 01_surv_apgi_calcs.Rout

02_surv_exp1.rda 02_surv_exp1.pdf: 02_surv_exp1.Rnw
	$(RSCRIPT) -e "knitr::knit2pdf('02_surv_exp1.Rnw')"

03_surv_exp2.rda 03_surv_exp2.pdf: 03_surv_exp2.Rnw ../biosurv/data/07_data_for_SIS.rda
	$(RSCRIPT) -e "knitr::knit2pdf('03_surv_exp2.Rnw')"
