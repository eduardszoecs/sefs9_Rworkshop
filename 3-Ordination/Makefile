all: slides purl

slides: slides/IntroOrdination.Rnw
				Rscript -e "knitr::knit(\"slides/IntroOrdination.Rnw\", output = \"slides/IntroOrdination.tex\")"
				Rscript -e "tools::texi2dvi(\"slides/IntroOrdination.tex\", pdf = TRUE, clean = TRUE)"
				mv IntroOrdination.pdf slides/IntroOrdination.pdf
        
purl:   slides/IntroOrdination.Rnw
				Rscript -e "knitr::purl(\"slides/IntroOrdination.Rnw\", output = \"code/IntroOrdination.R\", documentation = 0)"
