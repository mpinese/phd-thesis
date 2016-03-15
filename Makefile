.PHONY: all clean ultraclean dissertation rebuttal

all: dissertation rebuttal final

final: dissertation_final_master.pdf dissertation_final_public.pdf

dissertation_final_master.pdf: final_cover.pdf final_statements.pdf dissertation_temp_front.pdf dissertation_temp_back.pdf
	pdftk final_cover.pdf dissertation_temp_front.pdf final_statements.pdf dissertation_temp_back.pdf cat output dissertation_master.pdf

dissertation_final_public.pdf: dissertation_final_master.pdf
	ln -s dissertation_master.pdf dissertation_public.pdf

dissertation_temp_front.pdf: dissertation.pdf
	pdftk dissertation.pdf 2-4 output dissertation_temp_front.pdf

dissertation_temp_back.pdf: dissertation.pdf
	pdftk dissertation.pdf 7-end output dissertation_temp_back.pdf

dissertation: dissertation.pdf

rebuttal: rebuttal.pdf

rebuttal.pdf: rebuttal.tex
	-pdflatex -interaction=nonstopmode rebuttal.tex

%.pdf: %.tex dissertation.aux dissertation.glo
	-pdflatex -interaction=nonstopmode $*.tex

dissertation.aux dissertation.glo: *.tex dissertation.bib
	-latex --enable-write18 -interaction=nonstopmode dissertation.tex
	bibtex dissertation.aux
	makeglossaries dissertation
	-latex -interaction=nonstopmode dissertation.tex
	bibtex dissertation.aux
	makeglossaries dissertation
	-latex -interaction=nonstopmode dissertation.tex
	-latex -interaction=nonstopmode dissertation.tex

clean:
	rm -f *.aux *.bbl *.blg *.glg *.glo *.gls *.ist *.log *.out *.loa *.lof *.lot *.lox *.toc *.synctex.gz *.gnuplot *.table

ultraclean: clean
	rm dissertation.pdf rebuttal.pdf dissertation_final_master.pdf dissertation_final_public.pdf
