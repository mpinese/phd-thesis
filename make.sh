#!/bin/bash
mkdir -p build
pdflatex --output-directory build -interaction=nonstopmode thesis.tex
bibtex thesis.aux
makeglossaries thesis
pdflatex --output-directory build -interaction=nonstopmode thesis.tex
bibtex thesis.aux
makeglossaries thesis
pdflatex --output-directory build -synctex=1 -interaction=nonstopmode thesis.tex
mv build/thesis.pdf .
