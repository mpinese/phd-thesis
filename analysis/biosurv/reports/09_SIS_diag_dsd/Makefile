all: 09_SIS_diag_dsd.pdf

%.pdf : %.Rnw
	/home/marpin/bin/Rscript -e "knitr::knit2pdf('$<')"

clean:
	rm -f *.aux
	rm -f *.log
	rm -f *.tex
	rm -f *.pdf

ultraclean: clean
	rm -f *.pdf
	rm -rf cache/
	rm -rf figure/
