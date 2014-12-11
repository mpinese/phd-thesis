.PRECIOUS: image.rda

.PHONY: all clean ultraclean

all: report.pdf

report.pdf: report.Rnw
	Rscript -e "knitr::knit2pdf('report.Rnw')"

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
