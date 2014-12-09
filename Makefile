.PHONY: all clean ultraclean dissertation

all: dissertation

dissertation: dissertation.pdf

%.pdf: %.tex dissertation.aux dissertation.glo
	-pdflatex -interaction=nonstopmode $*.tex

dissertation.aux dissertation.glo: *.tex
	-pdflatex -interaction=nonstopmode dissertation.tex
	bibtex dissertation.aux
	makeglossaries dissertation
	-pdflatex -interaction=nonstopmode dissertation.tex
	bibtex dissertation.aux
	makeglossaries dissertation

clean:
	rm -f *.aux *.bbl *.blg *.glg *.glo *.gls *.ist *.log *.out *.loa *.lof *.lot *.lox *.toc *.synctex.gz

ultraclean: clean
	rm *.pdf
