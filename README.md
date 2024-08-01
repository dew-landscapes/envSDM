# `envSDM`: an R package to automate creation of species distribution models

The goal of `envSDM` is to help automate the preparation, tuning and prediction
of species distribution models. `envSDM` attempts to make decisions at each of
these steps that are robust(ish) for running SDMs for many, many taxa.

Preparation includes generating:

* optional generation of a (possibly buffered) minimum convex polygon around
presences to limit the rest of the process
* density raster of presences
* spatially thickened background points against density raster
* balanced spatial folds from the presences and background points
* environmental data for presences and background points (byo environmental
rasters)
* environmental layers that are correlated _at presences_

Tuning includes:

* three possible algorithms:
    + randomForest::randomForest()
        + always using the `randomForest()` `sampsize` argument downsample to
        the minimum number of presences
    + maxnet::maxnet()
    + predicts::envelope()
* ability to use multiple metrics for choosing a 'best' tune

Prediction includes choice of prediction to any combination of:

* full extent of environmental layers
* only within the minimum convex polygon around the presences
* only within another specified area of interest

## Installation

`envSDM` is not on [CRAN](https://CRAN.R-project.org).

You can install the development version from [GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("acanthiza/envSDM")
```
