ALL: all-redirect
SHELL = /bin/sh

srcdir = .
.SUFFIXES: .pdf .dvi .tex

VPATH=.:${srcdir}

LATEX      = TEXINPUTS=".:$(srcdir):" ; export TEXINPUTS ; latex
BIBTEX     = BIBINPUTS=".:$(srcdir):" ; export BIBINPUTS ; bibtex
DVIPDFM    = TEXINPUTS=".:$(srcdir):" ; export TEXINPUTS ; dvipdf
LATEX2HTML = latex2html

.tex.dvi:
	-$(LATEX) $< </dev/null >/dev/null 2>&1
	-$(LATEX) $< </dev/null >/dev/null 2>&1
	-$(BIBTEX) $* </dev/null
	$(LATEX) $<

.dvi.pdf:
	$(DVIPDFM) $*

all-redirect: design.pdf

clean:
	-rm -f *.toc *.out *.dvi user.ps user.pdf *.aux *.blg *.bbl *.log

mandoc:

htmldoc:

latexdoc: ALL
