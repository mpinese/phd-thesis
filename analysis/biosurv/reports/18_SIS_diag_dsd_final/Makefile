.PRECIOUS: image.rda

.PHONY: all clean ultraclean

all: report.pdf

report.pdf: image.rda report.Rnw
	Rscript -e "knitr::knit2pdf('report.Rnw')"

image.rda: ../../data/07_data_for_SIS.rda ../../common/08_SIS_common_funcs.R calcs.R
	Rscript calcs.R

clean:
	rm -f *.aux
	rm -f *.log
	rm -f *.tex
	rm -f *.pdf
	rm -rf NMF_*

ultraclean: clean
	rm -f *.pdf
	rm -rf cache/
	rm -rf figure/
