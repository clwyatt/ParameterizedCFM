SOURCE=paper

all: pdf html

prepare:
	dexy
	cp style.sty output/

pdf: prepare
	${MAKE} -C output -f Makefile.latex

html: prepare
	cd output && htlatex $(SOURCE).tex "paper"

clean:
	dexy reset

realclean: clean
	rm -rf artifacts logs *.pyc


