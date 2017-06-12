#MAKEFILE

#Name of tex file
FILE_NAME  = paper
FILE_TEX = $(FILE_NAME).tex
FILE_PDF = $(FILE_NAME).pdf

#Latex compiler
LATEX = pdflatex


all: $(FILE_PDF)

$(FILE_PDF) : $(FILE_TEX) references.bib
	$(LATEX) $(FILE_TEX)
	bibtex $(FILE_NAME)
	$(LATEX) $(FILE_TEX)
	$(LATEX) $(FILE_TEX)

diff: 
	latexdiff old_letter.tex letter.tex > diff_letter.tex
	$(LATEX) diff_$(FILE_TEX)
	bibtex $(FILE_NAME)
	$(LATEX) diff_$(FILE_TEX)
	$(LATEX) diff_$(FILE_TEX)

submit:
	tar -czvf submit.tgz paper.tex  CLARA-TOL-main.pdf P_cl.pdf emcee_results.pdf sigma_cl.pdf vinf_cl.pdf

clean:
		rm -f $(FILE_NAME).aux
		rm -f $(FILE_NAME).out
		rm -f $(FILE_NAME).bbl
		rm -f $(FILE_NAME).log
		rm -f $(FILE_NAME).synctex.gz
		rm -f $(FILE_NAME).blg
		rm -f $(FILE_NAME)Notes.bib
		rm -f $(FILE_NAME).pdf

