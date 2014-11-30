.PHONY: all clean ultraclean dissertation

all: dissertation

dissertation: dissertation.pdf

dissertation.pdf: *.tex
	mkdir -p build && \
	cp dissertation.bib build && \
	pdflatex --output-directory build -interaction=nonstopmode dissertation.tex && \
	cd build && \
	bibtex dissertation.aux && \
	makeglossaries dissertation && \
	cd .. && \
	pdflatex --output-directory build -interaction=nonstopmode dissertation.tex && \
	cd build && \
	bibtex dissertation.aux && \
	makeglossaries dissertation && \
	cd .. && \
	pdflatex --output-directory build -interaction=nonstopmode dissertation.tex && \
	mv build/dissertation.pdf .

clean:
	rm -rf build
	rm -f *.aux *.bbl *.blg *.glg *.glo *.gls *.ist *.log *.out *.loa *.lof *.lot *.lox *.toc *.synctex.gz

ultraclean: clean
	rm *.pdf
