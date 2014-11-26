#!/bin/bash
mkdir -p build
cp thesis.bib build
pdflatex --output-directory build -interaction=nonstopmode thesis.tex
cd build
bibtex thesis.aux
makeglossaries thesis
cd ..
pdflatex --output-directory build -interaction=nonstopmode thesis.tex
cd build
bibtex build/thesis.aux
makeglossaries thesis
cd ..
pdflatex --output-directory build -interaction=nonstopmode thesis.tex
mv build/thesis.pdf .
