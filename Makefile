.PHONY: all clean ultraclean dissertation

all: dissertation

dissertation: dissertation.pdf

dissertation.pdf: *.tex
	mkdir -p build
	cp *.tex build
	cp *.bib build
	cd build
	-pdflatex -interaction=nonstopmode dissertation.tex
	bibtex dissertation.aux
	makeglossaries dissertation
	-pdflatex -interaction=nonstopmode dissertation.tex
	bibtex dissertation.aux
	makeglossaries dissertation
	-pdflatex -interaction=nonstopmode dissertation.tex
	mv dissertation.pdf ..

clean:
	rm -rf build
	rm -f *.aux *.bbl *.blg *.glg *.glo *.gls *.ist *.log *.out *.loa *.lof *.lot *.lox *.toc *.synctex.gz

ultraclean: clean
	rm *.pdf
