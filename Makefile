SOURCE=paper

paper:
	dexy
	${MAKE} -C output -f Makefile.latex $(SOURCE)

html: 
	cd output && htlatex $(SOURCE).tex

clean:
	dexy reset

realclean: clean
	rm -rf artifacts logs


