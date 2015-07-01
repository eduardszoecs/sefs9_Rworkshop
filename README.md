Workshop: "Data analysis in freshwater ecology using R" @SEFS9
--------------------------------

This is the GitHub repository for the workshop  "Data analysis in freshwater ecology using R" held at the [9th Symposium for European Freshwater Sciences 2015](http://www.sefs9.ch/?page_id=1480).

* Date:   Sunday, **July 5th 2015**
* Time: **9.00 - 17.00** 
* Location:  Room **R160** on the ground floor of the Conference Center.



## Organizers

* [Jun.-Prof. Dr. Ralf B. Schäfer](https://www.uni-koblenz-landau.de/en/campus-landau/faculty7/environmental-sciences/landscape-ecology/Staff/ralf-schaefer/ralf-schaefer)
* [Avit Kumar Bhowmik](https://www.uni-koblenz-landau.de/en/campus-landau/faculty7/environmental-sciences/landscape-ecology/Staff/avitbhowmik)
* [Eduard Szöcs](https://www.uni-koblenz-landau.de/en/campus-landau/faculty7/environmental-sciences/landscape-ecology/Staff/eduardszoecs)


## Preliminaries

### Install R

If you don't already have R set up laptop, please download and install it from [here](http://cran.rstudio.com/). 
You can use the code editor of your choice. However, we recommend using [RStudio Desktop](http://www.rstudio.com/products/rstudio/download/).


### Installing packages

Packages extend the basic functionality of R and add functions or datasets.
For this course we need a few extra packages.  Please install the following packages - simply paste and run this command in your R prompt :

```{R}
install.packages(c("vegan", "rgl", "leaps", "car", "DAAG", "MuMIn", "relaimpo", 
                   "hier.part", "bootstrap", "MPV", "faraway"), 
                 dependencies = TRUE)
```


### Downloading code and data

1. [Click here](https://github.com/EDiLD/sefs9_Rworkshop/archive/master.zip) to download all the material.
2. Unzip the file - The unzipped folder contains all files and folders of this on-line repository.


## Structure of the Course

The course is structured into four parts of roughly 1.5 hours:

* [Linear and Generalized Linear Models (Ralf B. Schäfer)](https://github.com/EDiLD/sefs9_Rworkshop/tree/master/1-LinearModels) (**9.00 - 11.00**)
* [Ordinations (Eduard Szöcs)](https://github.com/EDiLD/sefs9_Rworkshop/tree/master/2-Ordination) (**11.15 - 12.15 and 13.15 - 14.45**)
* [Spatial data analysis (Avit Kumar Bhowmik)](https://github.com/EDiLD/sefs9_Rworkshop/tree/master/3-SpatialModels) (**15.00 - 16.45**)
* Course Evaluation (**16.45 - 17.00**) - see link below.

Each part has its own folder, with slides, data, code and other materials.
Each lecture is structured in: Intro - Demo - Exercise


## Feedback

Please rate our course by filling [this questionnaire](https://docs.google.com/forms/d/1wYkbajB_q5NPxzqJfS4gAwrRro7BcOd4OE_W5tmiOYo/viewform?usp=send_form). It is anonymous and will be only used to improve the course.


## License  
<a rel="license" href="http://creativecommons.org/licenses/by-sa/3.0/de/"><img alt="Creative Commons License" style="border-width:0" src="https://i.creativecommons.org/l/by-sa/3.0/de/88x31.png" /></a><br />This work is licensed under a <a rel="license" href="http://creativecommons.org/licenses/by-sa/3.0/de/">Creative Commons Attribution-ShareAlike 3.0 Germany License</a>.
