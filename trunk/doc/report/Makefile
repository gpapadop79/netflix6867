TARGET  = report
SRCS    = $(wildcard *.tex *.bib)
# need to manually activate the bibtex module because of the premature \end{document}
FLAGS   = -m bibtex -Wrefs -Wmisc

# garbage
# FIGS    = $(wildcard figures/*.pdf)
# BIBS    = $(wildcard *.bib)

all: clean mixture.tex pdf

mixture.tex: mixture.tex.bak
	../tools/FixDia ../tools/rules.txt < mixture.tex.bak > mixture.tex

ps:
	rubber --ps $(FLAGS) $(TARGET)

pdf:
	GS_OPTIONS=-dPDFSETTINGS=/prepress rubber --pdf $(FLAGS) $(TARGET) # embeds fonts

clean:
	rubber --ps --clean $(FLAGS) $(SRCS)
