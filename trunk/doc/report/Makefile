TARGET  = report
SRCS    = $(wildcard *.tex *.bib)
# FIGS    = $(wildcard figures/*.pdf)
FLAGS   = -m bibtex -Wrefs -Wmisc
# BIBS    = $(wildcard *.bib)

all: clean mixture.tex ps

ps:
	rubber --ps $(FLAGS) $(TARGET)
# need to manually activate the bibtex module because of the premature \end{document}

pdf:
	GS_OPTIONS=-dPDFSETTINGS=/prepress rubber --pdf $(FLAGS) $(SRCS) # embeds fonts

clean:
	rubber --ps --clean $(FLAGS) $(SRCS)
# obsolete: rm -f *.aux *.log *.out *.bbl *.blg *~ *.bak $(TARGET).ps $(TARGET).pdf
