TEXI2DVI = texi2dvi --pdf
SWEAVE = R CMD Sweave
RM = rm -rf
DISTR = probdistr-ref.bib 

all: probdistr-main.pdf clean

probdistr-main.pdf : probdistr-main.Rnw $(DISTR)
	$(SWEAVE) probdistr-main.Rnw
	$(TEXI2DVI) probdistr-main.tex
	$(TEXI2DVI) probdistr-main.tex

clean:
	$(RM) *.aux *.blg *.bbl *.brf *.out *.rel *.toc *~ auto

cleanall:
	$(RM) *.aux *.log *.blg *.bbl *.brf *.out *.rel *.toc *~ auto probdistr-main.pdf
