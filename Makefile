SOURCE=paper

paper:
	dexy
	cp Makefile.latex output/
	${MAKE} -C output -f Maefile.latex

html: 
	cd output && htlatex $(SOURCE).tex

clean:
	dexy reset

realclean: clean
	rm -rf artifacts log


