
earthtide
=========

[![Build Status](https://travis-ci.org/jkennel/earthtide.svg?branch=master)](https://travis-ci.org/jkennel/earthtide) [![Coverage Status](https://img.shields.io/codecov/c/github/jkennel/earthtide/master.svg)](https://codecov.io/github/jkennel/earthtide?branch=master)

Purpose
=======

Simplify generation of earth tides and leverage the model fitting capabilities of R.

Background
==========

The R earthtide package is a port of the Fortran ETERNA 3.4 (Wenzel, 1996) predict and analyze codes with the Kudryavtsev 2004 update. The original Fortran code was rewritten in R, and C++ using Rcpp, RcppArmadillo, and RcppParallel. The package is useful for generating synthetic earth tides using highly accurate tidal catalogs for prediction and regression applications in R. Attempts were made to ensure that results were consistent with the ETERNA 3.4 results, however, there is always the possibility that a bug was introduced in during the conversion and update. For the most feature rich version and up-to-date version of ETERNA please see <http://ggp.bkg.bund.de/eterna/> that is updated by Klaus Schueller.

Wenzel, H.G. 1996: The nanogal software: Earth tide data processing package ETERNA 3.30. Bull. Inf. Marges Terrestres. 124, 9425-9439.

Kudryavtsev, S.M., 2004. Improved harmonic development of the Earth tide-generating potential. J. Geod. 77, 829–838.

Installation
============

Example
=======

``` r
library(earthtide)

et <- Earthtide$new(
  utc = as.POSIXct("2015-01-01", tz = "UTC") + 0:(24*31) * 3600,
  latitude = 52.3868,
  longitude = 9.7144,
  catalog = "ksm03",
  freq_range = data.frame(start = 0.0, end = 6.0))

et$predict(method = "gravity", astro_update = 1)
et$lod_tide()
et$pole_tide()

head(et$output)
```

    ##              datetime    gravity  lod_tide pole_tide
    ## 1 2015-01-01 00:00:00 -161.41752 0.3261749 -3.158627
    ## 2 2015-01-01 01:00:00   89.96493 0.3243483 -3.169631
    ## 3 2015-01-01 02:00:00  321.51997 0.3225816 -3.180498
    ## 4 2015-01-01 03:00:00  500.73006 0.3208721 -3.191230
    ## 5 2015-01-01 04:00:00  610.79496 0.3192170 -3.201824
    ## 6 2015-01-01 05:00:00  652.84933 0.3176136 -3.212280

![](README_files/figure-markdown_github/plot-1.png)
