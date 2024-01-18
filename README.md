
<!-- README.md is generated from README.Rmd. Please edit that file -->

# UsefulPlants <img src="man/figures/logo.png" align="right" height="139"/>

<!-- badges: start -->

[![DOI](https://zenodo.org/badge/217888054.svg)](https://zenodo.org/badge/latestdoi/217888054)

<!-- badges: end -->

`UsefulPlants` is an R package primarily developed as supporting
software for the results and analyses conducted in the article **“The
global distribution of plants used by humans”** (Pironon et al. 2024).
It relies upon the R package `rsdm` (Ondo 2024) for building species
distribution models (SDMs) of plants used by humans at large scale, thus
it helps to streamline global biodiversity analyses, improve code
accessibility and reproducibility.

## *Installation*

Make sure to have [*R*](https://cloud.r-project.org/ "R") or
[*Rstudio*](https://rstudio.com/products/rstudio/download/ "Rstudio")
installed on your machine. Some R packages need to be compiled from
source, so if you are on Windows, you need to install
[*Rtools*](http://cran.r-project.org/bin/windows/Rtools/) too.

Then, install and load the package:

``` r
devtools::install_github("IanOndo/UsefulPlants") 
# optionally rebuild the vignettes by setting , build_opts = c("--no-resave-data", "--no-manual")

library(UsefulPlants)
```

Please follow the vignettes from the package
[**rsdm**](https://github.com/IanOndo/rsdm) for <u>*gathering and
curating plants occurrence records*</u> and for <u>*modelling species
distributions*</u>.

The following vignettes are available, but need to be (re-)built if
necessary.

## *Vignettes*

- Mapping species distribution, richness and endemism
- Exploring the latitudinal gradient of utilised plants species richness
  and endemism
- Exploring spatial correlations between utilised plants species
  richness/endemism and human cultural diversity (incoming soon !)

## References

<div id="refs" class="references csl-bib-body hanging-indent">

<div id="ref-rsdm" class="csl-entry">

Ondo, Ian. 2024. *<span class="nocase">rsdm: A R package for
streamlining large-scale studies of species distributions and
biodiversity patterns</span>*. <https://github.com/IanOndo/rsdm>.

</div>

<div id="ref-UsefulPlants" class="csl-entry">

Pironon, S., I. Ondo, M. Diazgranados, R. Allkin, A. C. Baquero, R.
Cámara-Leret, C. Canteiro, et al. 2024. “The Global Distribution of
Plants Used by Humans.” *Science* 383 (6680): 293–97.
<https://doi.org/10.1126/science.adg8028>.

</div>

</div>
