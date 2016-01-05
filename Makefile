.PHONY: all clean ultraclean dissertation rebuttal

all: dissertation rebuttal

dissertation: dissertation.pdf

rebuttal: rebuttal.pdf

rebuttal.pdf: rebuttal.tex
	-pdflatex -interaction=nonstopmode rebuttal.tex

%.pdf: %.tex dissertation.aux dissertation.glo
	-pdflatex -interaction=nonstopmode $*.tex

dissertation.aux dissertation.glo: *.tex dissertation.bib
	-pdflatex --enable-write18 -interaction=nonstopmode dissertation.tex
	bibtex dissertation.aux
	makeglossaries dissertation
	-pdflatex -interaction=nonstopmode dissertation.tex
	bibtex dissertation.aux
	makeglossaries dissertation

clean:
	rm -f *.aux *.bbl *.blg *.glg *.glo *.gls *.ist *.log *.out *.loa *.lof *.lot *.lox *.toc *.synctex.gz *.gnuplot *.table

ultraclean: clean
	rm dissertation.pdf rebuttal.pdf
