
<!-- README.md is generated from README.Rmd. Please edit that file -->

# UsefulPlants <img src="man/figures/logo.png" align="right" height="139"/>

<!-- badges: start -->

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.594920.svg)](https://doi.org/10.5281/zenodo.594920)

<!-- badges: end -->

`UsefulPlants` is an R package primarily developed to support the
results and analyses conducted in the article **“The global distribution
of plants used by humans”**. The package provides generic functions for
building species distribution models (SDMs) at large scale, thus it
helps streamline SDM analyses, improve code readability and
reproducibility.

## *Installation*

Make sure to have [*R*](https://cloud.r-project.org/ "R") or
[*Rstudio*](https://rstudio.com/products/rstudio/download/ "Rstudio")
installed on your machine. Some R packages need to be compiled from
source, so if you are on Windows, you need to install
[*Rtools*](http://cran.r-project.org/bin/windows/Rtools/) too.\\

Install *UsefulPlants* with the following instructions. If the package
`devtools` is not already installed run `install.packages("devtools")`
in your console. Setting `R_REMOTES_NO_ERRORS_FROM_WARNINGS="false` will
cause warning messages during calls to `devtools::install_github` to
become errors. So beforehand, make sure to set this environmental
variable to `true` via:

``` r
Sys.setenv(R_REMOTES_NO_ERRORS_FROM_WARNINGS="true")
```

Then, install and load the package:

``` r
devtools::install_github("IanOndo/UsefulPlants")
library(UsefulPlants)
```

## *Vignettes*

- Cleaning occurrence records
- Modelling species distribution
- Mapping species richness

## Reference
