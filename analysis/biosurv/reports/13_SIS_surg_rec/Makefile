all: 13_SIS_surg_rec.pdf

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
