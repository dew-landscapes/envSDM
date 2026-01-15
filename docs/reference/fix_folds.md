# Make sure all folds achieve min_fold_n presences

Make sure all folds achieve min_fold_n presences

## Usage

``` r
fix_folds(folds, pres, min_fold_n = 8, pres_val = 1)
```

## Arguments

- folds:

  Vector of fold ids

- pres:

  Vector of presence/absence data

- min_fold_n:

  Minimum number of presences in a fold

- pres_val:

  Value in `pres` that represent presence

## Value

Vector of adjusted fold ids ensuring each fold achieves min_fold_n
presences
