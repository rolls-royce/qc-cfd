#
# Makefile for Latex document
#

ROOTNAME = nozzle-util

default:
	pdflatex $(ROOTNAME).tex
	mv $(ROOTNAME).pdf ../.

index:
	makeindex nozzle-util.idx
	cp -f nozzle-util.ind index.tex

nomen:
	makeindex $(ROOTNAME).nlo -s nomencl.ist -o $(ROOTNAME).nls

bib:
	bibtex $(ROOTNAME)

pdf:
	#dvipdfm -p a4 review
	pdflatex $(ROOTNAME)

clean:
	rm  *.aux */*.aux *.log *.out *.ilg *.ind *.toc *.bbl *.blg index.tex

tar:
	cd ..; tar zcvf $(ROOTNAME)-latex.tar.gz $(ROOTNAME)

#
# Spellchecking from https://stackoverflow.com/questions/33847866/makefile-for-latex-aspell-with-multiple-files
#
# Type make spell
#
TEX_FILES = $(wildcard *.tex) $(wildcard */*.tex)
.PHONY: spell

spell: $(addsuffix .spchk,$(basename $(TEX_FILES)))

%.spchk: %.tex
	aspell --lang=en -t -c $<

