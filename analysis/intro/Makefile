.PHONY: all clean ultraclean figclean
.SECONDARY: 

#RSCRIPT = ~/bin/Rscript
RSCRIPT = Rscript

all: intro.pdf

clean:
	rm -f *.log *.aux *.pdf *.tex

ultraclean: clean figclean
	rm -rf cache/

figclean: clean
	rm -rf figure/ tikzMetrics

intro.pdf: cpvs.20150119.RData intro.Rnw
	$(RSCRIPT) -e "knitr::knit2pdf('intro.Rnw')"
