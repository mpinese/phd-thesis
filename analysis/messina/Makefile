.PHONY: all clean ultraclean figclean cacheclean texclean

RSCRIPT = Rscript

all: 05_E1.pdf 06_E2A.pdf 06_E2B.pdf

clean: texclean
	rm -f 01_surv_apgi_calcs.Rout 01_surv_apgi_calcs.rda 02_surv_exp1.rda .RData

texclean:
	rm -f 02_surv_exp1.pdf 03_surv_exp2.pdf 04_surv_exp3.pdf *.aux *.log *.tex tikzMetrics 05_E1.pdf 06_E2A.pdf 06_E2B.pdf

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

04_surv_exp3.pdf: 04_surv_exp3.Rnw
	$(RSCRIPT) -e "knitr::knit2pdf('04_surv_exp3.Rnw')"
	
05_E1.pdf: 05_E1.Rnw
	$(RSCRIPT) -e "knitr::knit2pdf('05_E1.Rnw')"
	
06_E2A.pdf: 06_E2A.Rnw
	$(RSCRIPT) -e "knitr::knit2pdf('06_E2A.Rnw')"

06_E2B.pdf: 06_E2B.Rnw
	$(RSCRIPT) -e "knitr::knit2pdf('06_E2B.Rnw')"
