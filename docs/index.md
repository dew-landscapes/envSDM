# References

The goal of `envSDM` is to help automate the preparation, tuning and
prediction of species distribution models. `envSDM` attempts to make
decisions at each of these steps that are robust(ish) for running SDMs
for many, many taxa.

If you are looking for packages to run species distribution models (or
ecological niche models) there are plenty of better packages to choose
from:

- [tidysdm](https://cran.r-project.org/web/packages/tidysdm/index.html)
- [flexsdm](https://sjevelazco.github.io/flexsdm/)
- [ENMeval](https://cran.r-project.org/web/packages/ENMeval/index.html)

`envSDM` assumes that your are trying to run many, many taxa, thus there
is no option to run a single taxa in parallel. The functions are all
designed around the potential to run many taxa in parallel (assuming
each taxa is run on a single core). For the long running functions,
there is the option to return either: the object, or the path to an .rds
file into which the object is saved.

Preparation includes generating:

- generation of a (possibly buffered) minimum convex polygon around
  presences to limit the rest of the process (predict boundary)
- density raster of presences
- spatially thickened (Vollering et al. 2019) background points against
  density raster
- balanced spatial folds from the presences and background points
- ability to run repeated spatial cross validation
- environmental data for presences and background points
- ensuring the environmental variables used are not correlated beyond a
  threshold, per taxa

Tuning includes:

- three possible algorithms:
  - [`randomForest::randomForest()`](https://rdrr.io/pkg/randomForest/man/randomForest.html)
    - always using the `randomForest()` `sampsize` argument downsample
      to the minimum number of presences
  - [`maxnet::maxnet()`](https://rdrr.io/pkg/maxnet/man/maxnet.html)
  - [`predicts::envelope()`](https://rdrr.io/pkg/predicts/man/envelope.html)
- ability to use multiple metrics for choosing a ‘best’ tune

Prediction is limited to boundary established during the preparation.

Several values provided as options for thresholding the prediction:

- maximum of specificity + sensitivity
- equal specificity + sensitivity
- omission rate at 10%
- equal prevalence
- no omission

## Installation

`envSDM` is not on [CRAN](https://CRAN.R-project.org).

You can install the development version from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("dew-landscapes/envSDM")
```

Vollering, Julien, Rune Halvorsen, Inger Auestad, and Knut Rydgren.
2019. “Bunching up the Background Betters Bias in Species Distribution
Models.” *Ecography* 42 (10): 1717–27.
<https://doi.org/><https://doi.org/10.1111/ecog.04503>.
