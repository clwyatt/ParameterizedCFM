SOURCE=paper

all: pdf html

prepare:
	dexy

pdf: prepare
	${MAKE} -C output -f Makefile.latex

html: prepare
	cd output && htlatex $(SOURCE).tex

clean:
	dexy reset

realclean: clean
	rm -rf artifacts logs


