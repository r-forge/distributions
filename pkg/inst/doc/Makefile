TEXI2DVI = texi2dvi --pdf
RM = rm -rf

all: probdistr-main.pdf clean

probdistr-main.pdf : probdistr-main.tex
	$(TEXI2DVI) probdistr-main.tex

clean:
	$(RM) *.aux *.log *.blg *.bbl *.brf *.out *.rel *.toc *~ auto