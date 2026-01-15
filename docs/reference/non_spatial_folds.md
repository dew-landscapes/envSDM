# Generate non spatial folds

Generate non spatial folds

## Usage

``` r
non_spatial_folds(use_folds, data, pa_col = "pa", pres_val = 1)
```

## Arguments

- use_folds:

  Numeric. Total number of folds

- data:

  Dataframe with `pa_col`

- pa_col:

  Name of column in data containing p/a values

- pres_val:

  Value in `pa_col` identifying presences
