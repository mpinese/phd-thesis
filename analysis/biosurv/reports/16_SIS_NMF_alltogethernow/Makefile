.PRECIOUS: image.rda

all: report.pdf

report.pdf: image.rda
	/home/marpin/bin/Rscript -e "knitr::knit2pdf('report.Rnw')"

image.rda: 
	/home/marpin/bin/Rscript calcs.R

clean:
	rm -f *.aux
	rm -f *.log
	rm -f *.tex
	rm -f *.pdf
	rm -rf NMF_*

ultraclean: clean
	rm -f *.pdf
	rm -f *.rda
	rm -rf cache/
	rm -rf figure/
