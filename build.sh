#!/bin/bash
pdflatex --output-directory out thesis.tex
mv out/thesis.pdf .
