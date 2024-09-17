The goal of `envSDM` is to help automate the preparation, tuning and prediction
of species distribution models. `envSDM` attempts to make decisions at each of
these steps that are robust(ish) for running SDMs for many, many taxa.

`envSDM` assumes that your are trying to run many, many taxa, thus there is no
option to run a single taxa in parallel. The functions are all designed around
the potential to run many taxa in parallel (assuming each taxa is run on a
single core).

Preparation includes generating:

* optional generation of a (possibly buffered) minimum convex polygon around
presences to limit the rest of the process
* density raster of presences
* spatially thickened background points against density raster
* balanced spatial folds from the presences and background points
* environmental data for presences and background points
* ensuring the environmental variables used are not correlated beyond a
threshold, per taxa

Tuning includes:

* three possible algorithms:
    + randomForest::randomForest()
        + always using the `randomForest()` `sampsize` argument downsample to
        the minimum number of presences
    + maxnet::maxnet()
    + predicts::envelope()
* ability to use multiple metrics for choosing a 'best' tune

Prediction includes:

* always predicting to the full extent of environmental layers
    + often this will include some rubbish prediction outside the predict
    boundary, but enables stacking the results later
* mask to the predict boundary (usually based on a, possibly buffered, minimum
convex polygon around presences)
    + masking back from the same full extent enables stacking the results for
    each taxa
* threshold to maximum of specificity + sensitivity

## Installation

`envSDM` is not on [CRAN](https://CRAN.R-project.org).

You can install the development version from [GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("acanthiza/envSDM")
```
