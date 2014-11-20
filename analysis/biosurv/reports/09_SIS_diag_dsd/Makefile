all: pdf

pdf: 09_SIS_diag_dsd.Rnw
	Rscript -e "knitr::knit2pdf('09_SIS_diag_dsd.Rnw')"

clean:
	rm -f *.aux
	rm -f *.log
	rm -f *.tex

ultraclean: clean
	rm -f 09_SIS_diag_dsd.pdf
	rm -rf cache/
	rm -rf figure/
