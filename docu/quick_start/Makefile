all: clean SICOPOLIS_V52_Quick_Start.pdf

pdf: SICOPOLIS_V52_Quick_Start.pdf

SICOPOLIS_V52_Quick_Start.pdf: SICOPOLIS_V52_Quick_Start.tex
	pdflatex SICOPOLIS_V52_Quick_Start.tex
	pdflatex SICOPOLIS_V52_Quick_Start.tex
	pdflatex SICOPOLIS_V52_Quick_Start.tex
	bibtex SICOPOLIS_V52_Quick_Start
	bibtex SICOPOLIS_V52_Quick_Start
	bibtex SICOPOLIS_V52_Quick_Start
	pdflatex SICOPOLIS_V52_Quick_Start.tex
	pdflatex SICOPOLIS_V52_Quick_Start.tex
	pdflatex SICOPOLIS_V52_Quick_Start.tex

clean:
	rm -f *.ps *.dvi *.aux *.toc *.idx *.ind *.ilg *.bbl *.blg \
              *.log *.out SICOPOLIS_V52_Quick_Start.pdf
