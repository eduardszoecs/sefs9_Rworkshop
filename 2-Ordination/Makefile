all: slides purl handout

slides: 	slides/IntroOrdination.Rnw
					Rscript -e "knitr::knit(\"slides/IntroOrdination.Rnw\", output = \"slides/IntroOrdination.tex\")"
					Rscript -e "tools::texi2dvi(\"slides/IntroOrdination.tex\", pdf = TRUE, clean = TRUE)"
					mv IntroOrdination.pdf slides/IntroOrdination.pdf
        
purl:   	slides/IntroOrdination.Rnw
					Rscript -e "knitr::purl(\"slides/IntroOrdination.Rnw\", output = \"code/IntroOrdination.R\", documentation = 0)"

handout:	slides/IntroOrdination.tex
					cat slides/IntroOrdination.tex | sed -e '1s/.*/\\documentclass[9pt, compress, table, handout]{beamer}/' > slides/handout_IntroOrdination.tex
					Rscript -e "tools::texi2dvi(\"slides/handout_IntroOrdination.tex\", pdf = TRUE, clean = TRUE)"
					pdfnup --frame true --a4paper --landscape --nup 1x2 --scale 0.85 --offset '-6cm 0cm' handout_IntroOrdination.pdf
					-mv handout_IntroOrdination-nup.pdf handout/handout_IntroOrdination.pdf
					-rm handout_IntroOrdination.pdf
					-rm -rf figure
					-rm -rf cache
