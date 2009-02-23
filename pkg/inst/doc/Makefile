TEXI2DVI = texi2dvi --pdf
SWEAVE = R CMD Sweave
RM = rm -rf

all: probdistr-main.pdf clean

probdistr-main.pdf : probdistr-main.Rnw
	$(SWEAVE) probdistr-main.Rnw
	$(TEXI2DVI) probdistr-main.tex

clean:
	$(RM) *.aux *.log *.blg *.bbl *.brf *.out *.rel *.toc *~ auto