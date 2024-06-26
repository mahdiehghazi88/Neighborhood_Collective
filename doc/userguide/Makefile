ALL: all-redirect
SHELL = /bin/sh
srcdir = .
MPICH_VERSION = 3.1.4
.SUFFIXES: .pdf .ps .dvi .tex

# Set these to the location of the web and ftp locations for the documents
WEBDIR = /mcs/www/mpi/mpich/docs
FTPDIR = /home/ftp/pub/mpi/mpich

VPATH=.:${srcdir}

LATEX      = TEXINPUTS=".:$(srcdir):" ; export TEXINPUTS ; latex
BIBTEX     = BIBINPUTS=".:$(srcdir):" ; export BIBINPUTS ; bibtex
DVIPDFM    = TEXINPUTS=".:$(srcdir):" ; export TEXINPUTS ; dvipdf
LATEX2HTML = latex2html

# Update the %MPICH_VERSION% in user.tex with current version string.
user.tex: user.tex.vin
	sed -e "s/%MPICH_VERSION%/${MPICH_VERSION}/g" $? > $@

user.dvi: user.tex user.bib
	-$(LATEX) user.tex
	-$(BIBTEX) user 
	-$(LATEX) user.tex
	$(LATEX)  user.tex

user.pdf: user.dvi
	$(DVIPDFM) user.dvi

# Use latex2html to create the HTML version.  There are several 
# steps needed to workaround features of latex2html
# 1. latex2html needs the user.aux file, so we run a latex step if
#    necessary
# 2. latex2html is not compatible with VPATH builds (it runs in the
#    directory that contains the source file), so we copy the source
#    file.  If more source files are added, this step will need to 
#    be changed
user.aux: user.tex
	$(LATEX) user.tex

user/user.html: user.aux user.tex user.bib
	$(LATEX2HTML) user.tex

# Here are some thoughts on using tohtml
#	tohtml -debugfile -default -dosnl -gaudy \
#              -basedef $(srcdir)/../mpiman.def $(srcdir)/user.tex
#	tohtml -default -dosnl -gaudy -basedef $(srcdir)/../mpiman.def \
#	$(srcdir)/user.tex

all-redirect: user.pdf

mandoc:

htmldoc:

latexdoc: ALL

clean:
	-rm -f *.toc *.out *.dvi user.pdf *.aux *.blg *.bbl *.log

distclean: clean

install: user.pdf user/user.html
	@if [ -z "$(WEBDIR)" ] ; then echo "No WEBDIR set" ; exit 1 ; fi
	@if [ ! -d $(WEBDIR) ] ; then echo "Create $(WEBDIR) first" ; exit 1 ;fi
	@if [ -z "$(FTPDIR)" ] ; then echo "No FTPDIR set" ; exit 1 ; fi
	@if [ ! -d $(FTPDIR) ] ; then echo "Create $(FTPDIR) first" ; exit 1 ;fi
	cp user.pdf $(WEBDIR)/user.pdf
	cp -rp user $(WEBDIR)/user-new
	rm -rf $(WEBDIR)/user
	mv -f $(WEBDIR)/user-new $(WEBDIR)/user
	cp user.pdf $(FTPDIR)/user.pdf
	cp -rp user $(FTPDIR)/user-new
	rm -rf $(WEBDIR)/user
	mv -f $(FTPDIR)/user-new $(FTPDIR)/user
