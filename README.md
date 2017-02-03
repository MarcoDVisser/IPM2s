IPM2s (0.1) 
===========

Archived code for building two stage Integral Projection Models in R. The example code included in /inst/example shows how to fitting vital rate models and sample from the posterior of mixed models.


## Installation

There is no  release on [CRAN](http://cran.r-project.org), but to install you can download the pacakage as [zip](https://github.com/MarcoDVisser/aprof/zipball/master) 
or [tar ball](https://github.com/MarcoDVisser/aprof/tarball/master). To install decompress these and run R CMD INSTALL on the contents of the archives, or use the **devtools** package to install directly from R.


```r
## devtools is required
require(devtools)
install_github("MarcoDVisser/aprof")
```
### Dependencies
The package depends on MASS and lme4 

