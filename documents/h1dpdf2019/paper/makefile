projectname = h1dpdf2019

all:
	pdflatex $(projectname)
	bibtex $(projectname)
	pdflatex $(projectname)

clean:
	rm -f *.aux *.out *.toc *.log $(projectname).blg $(projectname).bbl $(projectname).pdf

